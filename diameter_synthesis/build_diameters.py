import os
import glob
import shutil
import json
from tqdm import tqdm

import numpy as np
from numpy.polynomial import polynomial as polynomial
import random

from collections import deque

import neurom as nm

from diameter_synthesis import io

from diameter_synthesis.types import STR_TO_TYPES, TYPE_TO_STR
import diameter_synthesis.utils as utils
import diameter_synthesis.morph_functions as morph_funcs
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.utils import get_diameters, set_diameters
import diameter_synthesis.plotting as plotting

from multiprocessing import Pool
from functools import partial

TRUNK_FRAC_DECREASE = 0.1

##################################
## Build diameters from a model ##
##################################


def build_diameters(models, models_params, morphologies_dict, neurite_types, new_morph_path, extra_params, morph_path, plot=True, n_cpu=1, n_samples=1, ext = '.png'):
    """ Building the diameters from the generated diameter models"""

    all_models = {}
    for model in models:
        if model == 'M0':
            all_models[model] = diametrize_model_generic
        if model == 'M1':
            all_models[model] = diametrize_model_apical
        if model == 'M2':
            all_models[model] = diametrize_model_apical
        if model == 'M3':
            all_models[model] = diametrize_model_astrocyte

    # collect neurons paths and mtypes
    for model in models:
        print('Generating model', model)
        neurons = []
        for mtype in morphologies_dict:
            for neuron in morphologies_dict[mtype]:
                name, neuron_ext = os.path.splitext(neuron)
                if neuron_ext in {'.h5', '.asc', '.swc'} and os.path.exists(morph_path + '/' + neuron):
                    neurons.append([neuron, mtype])

        # set all parameters
        build_diam_poolf = partial(build_diam_pool, all_models, model, models_params,
                                   neurite_types, extra_params, morph_path, new_morph_path, plot, n_samples, ext)

        # generate diameters in parallel
        with Pool(processes=n_cpu) as p_build:  # initialise the parallel computation
            list(tqdm(p_build.imap(build_diam_poolf, neurons), total=len(neurons)))


def build_diam_pool(all_models, model, models_params, neurite_types, extra_params, morph_path, new_morph_path, plot, n_samples, ext, neuron_input):
    """ build a neuron diameters, save and plot it """

    fname = neuron_input[0]
    mtype = neuron_input[1]

    filepath = os.path.join(morph_path, fname)
    neuron = io.load_morphology(filepath)

    np.random.seed(extra_params[model]['seed'])

    all_models[model](neuron, models_params[mtype][model], neurite_types, extra_params[model])
    if n_samples > 1:
        diameters = utils.get_all_diameters(neuron)
        for i in range(n_samples - 1):
            all_models[model](neuron, models_params[mtype][model],
                              neurite_types, extra_params[model])
            diameters_tmp = utils.get_all_diameters(neuron)
            for di, diams in enumerate(diameters_tmp):
                diameters[di] += diams

        for di, diams in enumerate(diameters):
            diams /= n_samples

        utils.set_all_diameters(neuron, diameters)

    io.save_neuron(neuron, model, new_morph_path)
    if plot:
        folder = 'shapes_' + os.path.basename(new_morph_path[:-1])
        plotting.plot_diameter_diff(os.path.splitext(fname)[0], morph_path, neuron,
                                    model, neurite_types, folder=folder, ext = ext)


def diametrize_model_astrocyte(neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type ==
                    STR_TO_TYPES[neurite_type])

        for neurite in neurites:

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            n_tries_step = 1
            while wrong_tips:
                # sample a trunk diameter
                trunk_diam = trunk_diam_frac * get_trunk_diameter(neurite, params['trunk_diameter'][neurite_type])
                # try to diametrize the neurite
                wrong_tips = diametrize_tree(neurite, params, neurite_type, trunk_diam, mode_sibling='threshold', mode_rall='generic', sibling_threshold=extra_params['threshold'][neurite_type], rall_threshold=extra_params['threshold'][neurite_type], with_asymmetry=True, no_taper=False, reduction_factor_max=2.0)

                # if we can't get a good model, reduce the trunk diameter progressively
                n_tries += 1
                if n_tries > 2 * n_tries_step:  # if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1

                # don't try to much and keep the latest try
                if n_tries > extra_params['trunk_max_tries'] and extra_params['trunk_max_tries'] > 1:
                    print('max tries attained with', neurite_type)
                    wrong_tips = False



def diametrize_model_generic(neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type ==
                    STR_TO_TYPES[neurite_type])

        for neurite in neurites:

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            n_tries_step = 1
            while wrong_tips:

                # sample a trunk diameter
                trunk_diam = trunk_diam_frac * get_trunk_diameter(neurite, params['trunk_diameter'][neurite_type])
                # try to diametrize the neurite
                wrong_tips = diametrize_tree(neurite, params, neurite_type, trunk_diam, mode_sibling='generic', mode_rall='generic', sibling_threshold=extra_params['threshold'][neurite_type], rall_threshold=extra_params['threshold'][neurite_type], with_asymmetry=True, no_taper=False, reduction_factor_max=1.0)

                # if we can't get a good model, reduce the trunk diameter progressively
                n_tries += 1
                if n_tries > 2 * n_tries_step:  # if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1

                # don't try to much and keep the latest try
                if n_tries > extra_params['trunk_max_tries'] and extra_params['trunk_max_tries'] > 1:
                    print('max tries attained with', neurite_type)
                    wrong_tips = False


def diametrize_model_apical(neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type ==
                    STR_TO_TYPES[neurite_type])

        for neurite in neurites:

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            n_tries_step = 1
            while wrong_tips:

                # sample a trunk diameter
                trunk_diam = trunk_diam_frac * get_trunk_diameter(neurite, params['trunk_diameter'][neurite_type])
                # try to diametrize the neurite
                wrong_tips = diametrize_tree(neurite, params, neurite_type, trunk_diam, mode_sibling='threshold', mode_rall='threshold', sibling_threshold=extra_params['threshold'][neurite_type], rall_threshold=extra_params['threshold'][neurite_type], with_asymmetry=True, no_taper=True, reduction_factor_max=1.0)

                # if we can't get a good model, reduce the trunk diameter progressively
                n_tries += 1
                if n_tries > 2 * n_tries_step:  # if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1

                # don't try to much and keep the latest try
                if n_tries > extra_params['trunk_max_tries'] and extra_params['trunk_max_tries'] > 1:
                    print('max tries attained with', neurite_type)
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


def get_daughter_diameters(section, terminal_diam, params, neurite_type, mode_sibling='generic', mode_rall='generic', tot_length=1, sibling_threshold=0.3, rall_threshold=0.3, with_asymmetry=False, reduction_factor_max=3.0):
    """ return daughter diamters from parent d0 """

    reduction_factor = reduction_factor_max + 1.0
    while reduction_factor > reduction_factor_max:  # try until we get a reduction of diameter in the branching

        sibling_ratio = get_sibling_ratio(section, params['sibling_ratio'][neurite_type], mode=mode_sibling, tot_length=tot_length, threshold=sibling_threshold)
        rall_deviation = get_rall_deviation(
            section, params['rall_deviation'][neurite_type], mode=mode_rall, tot_length=tot_length, threshold=rall_threshold)

        # compute the reduction factor
        reduction_factor = morph_funcs.rall_reduction_factor(rall_deviation=rall_deviation, siblings_ratio=sibling_ratio)

    d0 = get_diameters(section)[-1]
    # if new terminal diam is too large from the previous section, reassign it
    terminal_diam = min(d0, terminal_diam)

    d1 = reduction_factor * d0
    d2 = sibling_ratio * d1

    # set minimum values if too small
    d1 = max(d1, terminal_diam)
    d2 = max(d2, terminal_diam)

    #if reduction_factor == 1:
    #    d2 *= 1.3

    diams = [d1] + (len(section.children) - 1) * [d2]

    if with_asymmetry:
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


def diametrize_tree(neurite, params, neurite_type, trunk_diam, mode_sibling='generic', mode_rall='generic', sibling_threshold=0.3, rall_threshold=0.3, with_asymmetry=False, no_taper=False, reduction_factor_max=1.0):
    """ diametrize a single tree """

    # initialise status variables
    wrong_tips = False  # used to run until terminal diameters are thin enough

    # create a deque for breath-first
    active = deque([neurite.root_node])

    tot_length = nm.get('total_length', neurite)[0]

    while active:
        section = active.popleft()

        # set trunk diam if trunk, or first diam of the section otherwise
        if section.is_root():
            init_diam = trunk_diam
        else:
            init_diam = get_diameters(section)[0]

        # sample a terminal diameter
        terminal_diam = get_terminal_diameter(neurite, params['terminal_diameter'][neurite_type])
        taper = get_taper(neurite, params['taper'][neurite_type], no_taper=no_taper)

        # diametrize a section
        diametrize_section(section, init_diam, taper=taper, min_diam=terminal_diam, max_diam=trunk_diam)

        # if branching points has children, keep looping
        if len(section.children) > 0:
            diams = get_daughter_diameters(section, terminal_diam, neurite_type=neurite_type, params=params, mode_sibling=mode_sibling, mode_rall=mode_rall, tot_length=tot_length, sibling_threshold=sibling_threshold, with_asymmetry=with_asymmetry, reduction_factor_max=reduction_factor_max)

            # set diameters
            for i, ch in enumerate(section.children):
                utils.redefine_diameter_section(ch, 0, diams[i])
                active.append(ch)  # add sections to queue

        # if we are at a tip, check tip diameters to restart if too large
        elif get_diameters(section)[-1] > params['terminal_diameter'][neurite_type]['params']['max']:
            wrong_tips = True

    return wrong_tips


def diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=10.):
    '''Corrects the diameters of a section'''

    #max_diam = np.clip(np.random.normal(1.7, 1.0), 1.2, 2.5)
    #min_diam = np.clip(np.random.normal(0.6, 0.5), 0.3, 0.8)

    diams = [initial_diam]

    # if the initial diameter is not in the range of the sampled terminal diameters, just reset it
    #fact = 1.0
    #fact2 = 1.0

    if initial_diam < min_diam:
        min_diam = initial_diam

    #if initial_diam > max_diam:
    #    max_diam = initial_diam

    # lengths of each segments will be used for scaling of tapering
    lengths = [0] + utils.section_lengths(section)

    diams = polynomial.polyval(lengths, [initial_diam, taper])
    diams = np.clip(diams, min_diam, max_diam) 

    set_diameters(section, np.array(diams, dtype=np.float32))


def replace_params(param_type, param_all_name='./model_params_all.json', model='M0'):
    """replace the param file with the all neuron fit"""

    with open(param_all_name, 'r') as f:
        params_all = json.load(f)

    return params_all['all_types'][model][param_type]
