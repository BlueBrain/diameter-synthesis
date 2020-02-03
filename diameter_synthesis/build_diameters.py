""" Build neurite diameters from a pre-generated model """
import json
import os
import random
import time
from collections import deque
from functools import partial
import multiprocessing
import logging

import numpy as np
from numpy.polynomial import polynomial

import neurom as nm

import diameter_synthesis.morph_functions as morph_funcs
import diameter_synthesis.plotting as plotting
import diameter_synthesis.utils as utils
from diameter_synthesis import io
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.types import STR_TO_TYPES
from diameter_synthesis.utils import get_diameters, set_diameters

TRUNK_FRAC_DECREASE = 0.1
L = logging.getLogger(__name__)

################################
# Build diameters from a model #
################################


class Worker:
    """worker for building diameters"""
    def __init__(self, model, models_params, config):
        self.model = model
        self.models_params = models_params
        self.config = config

    def __call__(self, neuron_input):

        fname = neuron_input[0]
        mtype = neuron_input[1]

        L.info("%s...", fname)
        time_0 = time.time()

        neuron = io.load_morphology(os.path.join(self.config['morph_path'], fname))

        build(neuron, self.models_params[mtype][self.model], self.config)

        io.save_neuron(neuron, self.model, self.config['new_morph_path'])

        L.info("%s... done in %.2f seconds", fname, np.round(time.time() - time_0, 2))


def build_diameters(morphologies_dict, models_params, config):
    """ Building the diameters from the generated diameter models with multiprocessing"""

    all_models = {}
    for model in config['models']:
        all_models[model] = select_model(model)

    # collect neurons paths and mtypes
    for model in config['models']:
        L.info('Generating model %s', model)
        neurons = []
        for mtype in morphologies_dict:
            for neuron in morphologies_dict[mtype]:
                _, neuron_ext = os.path.splitext(neuron)
                if (neuron_ext in {'.h5', '.asc', '.swc'} and
                        os.path.exists(os.path.join(config['morph_path'], neuron))):
                    neurons.append([neuron, mtype])

        # generate diameters in parallel
        worker = Worker(model, models_params, config)
        pool = multiprocessing.Pool(config['n_cpu'])
        pool.map(worker, neurons)


def build(neuron, models_params, config):
    """ Building the diameters from the generated diameter models of a neuron"""

    model = config['models'][0]
    neurite_types = config['neurite_types']
    extra_params = config['extra_params']
    n_samples = config['n_samples']

    if 'seed' in extra_params[model]:
        np.random.seed(extra_params[model]['seed'])

    diameter_generator = select_model(model)

    diameter_generator(neuron, models_params, neurite_types, extra_params[model])

    if n_samples > 1:
        diameters = utils.get_all_diameters(neuron)
        for _ in range(n_samples - 1):
            diameter_generator(neuron, models_params, neurite_types, extra_params[model])
            diameters_tmp = utils.get_all_diameters(neuron)
            for j, diams in enumerate(diameters_tmp):
                diameters[j] += diams

        for diams in diameters:
            diams /= n_samples

        utils.set_all_diameters(neuron, diameters)

    if 'plot' in config:
        try:
            folder = 'shapes_' + os.path.basename(config['new_morph_path'][:-1])
            plotting.plot_diameter_diff(neuron.name, config['morph_path'],
                                        neuron, model, config['neurite_types'],
                                        folder=folder, ext=config['ext'])
        except Exception as exc:  # pylint: disable=broad-except
            L.exception(neuron.name, neuron, exc)


###################
# diameter models #
###################

def select_model(model):
    """ select a model to use from available ones """
    if model == 'M1':

        params_tree = {}
        params_tree['mode_sibling'] = 'threshold'
        params_tree['mode_rall'] = 'threshold'
        params_tree['with_asymmetry'] = True
        params_tree['no_taper'] = False
        params_tree['reduction_factor_max'] = 1.0

    return partial(diametrize_model, params_tree)


def diametrize_model(params_tree, neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    with_trunk_taper = False  # experimental trunk tapering

    for neurite_type in neurite_types:

        params_tree['neurite_type'] = neurite_type
        params_tree['sibling_threshold'] = extra_params['threshold'][neurite_type]
        params_tree['rall_threshold'] = extra_params['threshold'][neurite_type]

        neurites = (neurite for neurite in neuron.neurites if neurite.type ==
                    STR_TO_TYPES[neurite_type])

        for neurite in neurites:

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.0
            n_tries_step = 1
            while wrong_tips:

                # sample a trunk diameter
                trunk_diam = trunk_diam_frac * \
                    get_trunk_diameter(neurite, params['trunk_diameters'][neurite_type])
                if with_trunk_taper:
                    trunk_taper = -get_trunk_taper(neurite, params['trunk_tapers'][neurite_type])
                else:
                    trunk_taper = None

                if trunk_diam < 0.01:
                    trunk_diam = 1.0
                    L.warning('sampled trunk diameter < 0.01, so use 1 instead')

                params_tree['trunk_diam'] = trunk_diam
                params_tree['trunk_taper'] = trunk_taper

                wrong_tips = diametrize_tree(neurite, params, params_tree)

                # if we can't get a good model, reduce the trunk diameter progressively
                n_tries += 1
                if n_tries > 2 * n_tries_step:  # if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1
                # don't try to much and keep the latest try
                if (n_tries > extra_params['trunk_max_tries'] and
                        extra_params['trunk_max_tries'] > 1):
                    L.warning('max tries attained with %s', neurite_type)
                    wrong_tips = False


def get_sibling_ratio(section, params, mode='generic', tot_length=1., threshold=0.3):
    """return a sampled sibling ratio"""
    if mode == 'generic':
        seq_value = morph_funcs.sequential_single(params['sequential'], section=section)
        sibling_ratio = sample_distribution(params, seq_value)

    elif mode == 'threshold':
        # use a threshold to make sibling ratio = 0
        seq_value = morph_funcs.sequential_single(params['sequential'], section=section)
        seq_value /= tot_length

        if seq_value > threshold:
            sibling_ratio = 0.0
        else:
            sibling_ratio = sample_distribution(params, 0)  # sample from smallest distribution
    else:
        raise Exception('type not understood')

    return sibling_ratio


def get_rall_deviation(section, params, mode='generic', tot_length=1., threshold=0.3):
    """return a sampled rall deviation"""
    if mode == 'generic':
        # sample a Rall deviation
        seq_value = morph_funcs.sequential_single(params['sequential'], section=section)
        seq_value /= tot_length
        rall_deviation = sample_distribution(params, seq_value)

    elif mode == 'exact':
        rall_deviation = 1.

    elif mode == 'threshold':
        seq_value = morph_funcs.sequential_single(params['sequential'], section=section)
        seq_value /= tot_length
        if seq_value > threshold:
            rall_deviation = 1.0
        else:
            rall_deviation = sample_distribution(params, seq_value)
    else:
        raise Exception('type not understood')

    return rall_deviation


def get_trunk_diameter(neurite, params):
    """ sample a trunk diameter """
    seq_value = morph_funcs.sequential_single(params['sequential'], neurite=neurite)
    trunk_diam = sample_distribution(params, seq_value[0])

    return trunk_diam


def get_trunk_taper(neurite, params):
    """ sample a trunk taper """
    seq_value = morph_funcs.sequential_single(params['sequential'], neurite=neurite)
    trunk_taper = sample_distribution(params, seq_value[0])

    return trunk_taper


def get_terminal_diameter(neurite, params):
    """ sample a terminal diameter """

    seq_value = morph_funcs.sequential_single(params['sequential'], neurite=neurite)
    terminal_diam = sample_distribution(params, seq_value[0])

    return terminal_diam


def get_taper(section, params, no_taper=False):
    """ sample a taper """

    if no_taper:
        return 0.0

    # sample a taper rate
    seq_value = morph_funcs.sequential_single(params, section=section)
    taper = sample_distribution(params, seq_value[0])

    return taper


def get_daughter_diameters(section, params, params_tree):
    """ return daughter diamters from parent d0 """

    reduction_factor = params_tree['reduction_factor_max'] + 1.0
    # try until we get a reduction of diameter in the branching
    while reduction_factor > params_tree['reduction_factor_max']:
        sibling_ratio = get_sibling_ratio(
            section,
            params['sibling_ratios'][params_tree['neurite_type']],
            mode=params_tree['mode_sibling'],
            tot_length=params_tree['tot_length'],
            threshold=params_tree['sibling_threshold'])
        rall_deviation = get_rall_deviation(
            section,
            params['rall_deviations'][params_tree['neurite_type']],
            mode=params_tree['mode_rall'],
            tot_length=params_tree['tot_length'],
            threshold=params_tree['rall_threshold'])

        # compute the reduction factor
        reduction_factor = morph_funcs.rall_reduction_factor(
            rall_deviation=rall_deviation, sibling_ratio=sibling_ratio)

    diam_0 = get_diameters(section)[-1]

    # if new terminal diam is too large from the previous section, reassign it
    terminal_diam = min(diam_0, params_tree['terminal_diam'])

    diam_1 = reduction_factor * diam_0
    diam_2 = sibling_ratio * diam_1

    # set minimum values if too small
    diam_1 = max(diam_1, terminal_diam)
    diam_2 = max(diam_2, terminal_diam)

    diams = [diam_1] + (len(section.children) - 1) * [diam_2]

    if params_tree['with_asymmetry']:
        # if we want to use the asymmetry to set diameters, re-oreder them
        part = []
        for child in section.children:
            part.append(sum(1 for _ in child.ipreorder()))

        # sort them by larger first and create a list of diameters
        child_sort = np.argsort(part)[::-1]
        diams = list(np.array(diams)[child_sort])

    else:
        # otherwise just shuffle them to avoid systematic bias
        random.shuffle(diams)

    return diams


def diametrize_tree(neurite, params, params_tree):
    """ diametrize a single tree """

    neurite_type = params_tree['neurite_type']

    # initialise status variables
    wrong_tips = False  # used to run until terminal diameters are thin enough

    # create a deque for breath-first
    active = deque([neurite.root_node])

    params_tree['tot_length'] = nm.get('total_length', neurite)[0]

    while active:
        section = active.popleft()

        # set trunk diam if trunk, or first diam of the section otherwise
        if section.is_root():
            init_diam = params_tree['trunk_diam']
            if params_tree['trunk_taper'] is not None:
                taper = params_tree['trunk_tapers']
            else:
                taper = get_taper(neurite, params['tapers'][neurite_type],
                                  no_taper=params_tree['no_taper'])
        else:
            init_diam = get_diameters(section)[0]
            taper = get_taper(neurite, params['tapers'][neurite_type],
                              no_taper=params_tree['no_taper'])

        # sample a terminal diameter
        params_tree['terminal_diam'] = \
            get_terminal_diameter(neurite, params['terminal_diameters'][neurite_type])

        # diametrize a section
        diametrize_section(section, init_diam, taper=taper,
                           min_diam=params_tree['terminal_diam'],
                           max_diam=params_tree['trunk_diam'])

        # if branching points has children, keep looping
        max_diam = params['terminal_diameters'][neurite_type]['params']['max']
        if len(section.children) > 0:
            diams = get_daughter_diameters(section, params, params_tree)

            # set diameters
            for i, child in enumerate(section.children):
                utils.redefine_diameter_section(child, 0, diams[i])
                active.append(child)  # add sections to queue

        # if we are at a tip, check tip diameters to restart if too large
        elif get_diameters(section)[-1] > max_diam:
            wrong_tips = True

    return wrong_tips


def diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=10.):
    '''Corrects the diameters of a section'''

    diams = [initial_diam]

    # lengths of each segments will be used for scaling of tapering
    lengths = [0] + utils.section_lengths(section)

    diams = polynomial.polyval(lengths, [initial_diam, taper])
    diams = np.clip(diams, min_diam, max_diam)

    set_diameters(section, np.array(diams, dtype=np.float32))


def replace_params(param_type, param_all_name='./model_params_all.json', model='M0'):
    """replace the param file with the all neuron fit"""

    with open(param_all_name, 'r') as filename:
        params_all = json.load(filename)

    return params_all['all_types'][model][param_type]
