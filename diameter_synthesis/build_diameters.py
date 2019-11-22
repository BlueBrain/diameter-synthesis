import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np
from numpy.polynomial import polynomial as polynomial

import neurom as nm

from neurom import COLS
from neurom.core import iter_sections
from neurom import viewer

from diameter_synthesis import io

from diameter_synthesis.types import STR_TO_TYPES, TYPE_TO_STR
import diameter_synthesis.utils as utils 
import diameter_synthesis.morph_functions as morph_funcs
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.utils import get_diameters, set_diameters 
import diameter_synthesis.plotting as plotting

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from multiprocessing import Pool
from functools import partial

##################################
## Build diameters from a model ##
##################################

def build_diameters(models, models_params, morphologies_dict, neurite_types, new_morph_path, extra_params, morph_path, plot = True, n_cpu = 1):
    """ Building the diameters from the generated diameter models"""

    all_models = {}
    for model in models:
        all_models[model] = diametrize_model_generic

    #collect neurons paths and mtypes
    for model in models:
        print('Generating model', model)
        neurons = []
        for mtype in morphologies_dict:
            for neuron in morphologies_dict[mtype]:
                name, ext = os.path.splitext(neuron)
                if ext in {'.h5', '.asc', '.swc'} and os.path.exists(morph_path + '/' + neuron):
                    neurons.append([neuron, mtype])
        
        #set all parameters
        build_diam_poolf = partial(build_diam_pool, all_models, model, models_params, neurite_types, extra_params, morph_path, new_morph_path, plot)

        #generate diameters in parallel
        with Pool(processes = n_cpu) as p_build:  #initialise the parallel computation
            list(tqdm(p_build.imap(build_diam_poolf, neurons), total = len(neurons)))

def build_diam_pool(all_models, model, models_params, neurite_types, extra_params, morph_path, new_morph_path, plot, neuron_input):
    """ build a neuron diameters, save and plot it """

    fname = neuron_input[0]
    mtype = neuron_input[1]
    name, ext = os.path.splitext(fname)

    filepath = os.path.join(morph_path, fname)
    neuron = io.load_morphology(filepath)

    all_models[model](neuron, models_params[mtype][model], neurite_types, extra_params[model])

    io.save_neuron([neuron, name], model, new_morph_path)

    if plot:
        folder = 'shapes_' + os.path.basename(new_morph_path[:-1])
        plotting.plot_diameter_diff(name, morph_path, new_morph_path, model, neurite_types, folder=folder)

def diametrize_model_generic(neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    np.random.seed(extra_params['seed'])

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type == STR_TO_TYPES[neurite_type])

        for neurite in neurites:

            max_path_dist = np.max(nm.get('section_path_distances', neurite))

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            k = 1
            while wrong_tips:
                wrong_tips = diametrize_tree(neurite, params, neurite_type, max_path_dist, trunk_diam_frac)
                n_tries += 1

                if n_tries > 10*k: #if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= 0.1
                    k+=1

                if n_tries > 90: #don't try to much, and complain
                    print('max tries attained with', neurite_type)
                    wrong_tips = False

def replace_params(param_type, param_all_name='./model_params_all.json', model='M0'): 
    """replace the param file with the all neuron fit"""

    with open(param_all_name, 'r') as f:
        params_all = json.load(f)

    return params_all['all_types'][model][param_type]


def diametrize_tree(neurite, params, neurite_type, max_path_dist, trunk_diam_frac = 1.):
        """ diametrize a single tree """

        if params['trunk_diameter'][neurite_type]['distribution'] == 'exponnorm_sequence':
            if params['trunk_diameter'][neurite_type]['params']['a'][0] == 0:
                params_tmp = replace_params('trunk_diameter')[neurite_type]
            else:
                params_tmp = params['trunk_diameter'][neurite_type]
        else:
            params_tmp = params['trunk_diameter'][neurite_type]

        tpe = morph_funcs.sequential_single(params_tmp['sequential'], neurite = neurite)
        trunk_diam = trunk_diam_frac*sample_distribution(params_tmp, tpe[0])

        wrong_tips = False

        status = {s.id: False for s in iter_sections(neurite)}
        active = [neurite.root_node]

        while active:
            for section in list(active):

                if section.is_root():
                    init_diam = trunk_diam
                else:
                    init_diam = get_diameters(section)[0]
                
                #sample a terminal diameter
                if params['terminal_diameter'][neurite_type]['params']['a'] == 0:
                    params_tmp = replace_params('terminal_diameter')[neurite_type]
                else:
                    params_tmp = params['terminal_diameter'][neurite_type]

                min_diam = params_tmp['params']['min']
                max_diam = params_tmp['params']['max']
                tpe = morph_funcs.sequential_single(params_tmp['sequential'],neurite = neurite)
                terminal_diam = sample_distribution(params_tmp, tpe[0])
                
                #diametrize a section
                if params['taper'][neurite_type]['params']['a'] == 0:
                    params_tmp = replace_params('taper')[neurite_type]
                else:
                    params_tmp = params['taper'][neurite_type]

                tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                taper = -(sample_distribution(params_tmp, tpe[0])) #prevent positive tapers
                diametrize_section(section, init_diam, taper=taper,
                                             min_diam = terminal_diam, max_diam = trunk_diam)

                status[section.id] = True  # Tapering of section complete.
                active.remove(section)
                children = np.array(section.children)
                #if branching points has children, keep looping
                if len(children) > 1:
                    
                    reduc = 2.
                    while reduc > 1.: #try until we get a reduction of diameter in the branching
                        if params['sibling_ratio'][neurite_type]['params']['scale'] == 0.0:
                            params_tmp = replace_params('sibling_ratio')[neurite_type]
                        else:
                            params_tmp = params['sibling_ratio'][neurite_type]
                        tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                        sibling_ratio = sample_distribution(params_tmp, tpe[0])

                        if params['Rall_deviation'][neurite_type]['params']['a'] == 0:
                            params_tmp = replace_params('Rall_deviation')[neurite_type]
                        else:
                            params_tmp = params['Rall_deviation'][neurite_type]

                        tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                        Rall_deviation = sample_distribution(params_tmp, tpe[0])

                        reduc = morph_funcs.Rall_reduction_factor(Rall_deviation = Rall_deviation, siblings_ratio = sibling_ratio)

                    d0 = get_diameters(section)[-1]
                    if d0 < terminal_diam:
                        terminal_diam = d0

                    d1 = reduc * d0
                    d2 = sibling_ratio * d1

                    if d1 < terminal_diam:
                        d1 = terminal_diam

                    if d2 < terminal_diam:
                        d2 = terminal_diam

                    #if the asymetry pair is opposite to d1>d2, swith diameters
                    ass_pair = morph_funcs.sequential_single('asymmetry_pair', section = section)[0]
                    if ass_pair[1]> ass_pair[0]:
                        d1, d2 = d2, d1

                    #if len(children)>2:
                    #    print(len(children), 'children for this branch')

                    for i, ch in enumerate(children):
                        new_diam = d1 if i == 0 else d2
                        utils.redefine_diameter_section(ch, 0, new_diam)

                        active.append(ch)

                #if we are at a tip, check tip diameters and stop
                else:
                    if get_diameters(section)[-1] < min_diam or get_diameters(section)[-1]  > max_diam:
                        wrong_tips = True

        return wrong_tips


def diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=100.):
    '''Corrects the diameters of a section'''

    diams = [initial_diam]

    #if the initial diameter is not in the range of the sampled terminal diameters, just reset it
    if initial_diam < min_diam:
        min_diam = initial_diam

    if initial_diam > max_diam:
        max_diam = initial_diam

    # lengths of each segments will be used for scaling of tapering
    lengths = [0] + utils.section_lengths(section)

    diams = polynomial.polyval(lengths, [initial_diam, taper])
    diams[diams<min_diam] = min_diam
    diams[diams>max_diam] = max_diam
    set_diameters(section, np.array(diams, dtype=np.float32))
