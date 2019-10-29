import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm

from neurom import COLS
from neurom.core import iter_sections
from neurom import viewer

from utils import STR_TO_TYPES, TYPE_TO_STR

import diameter_synthesis.utils as utils 
import diameter_synthesis.morph_functions as morph_funcs
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.utils import get_diameters, set_diameters 
import diameter_synthesis.plotting as plotting

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import copy
##################################
## Build diameters from a model ##
##################################

def build_diameters(models, models_params, morphologies, neurite_types, new_morph_path, extra_params, morph_path):
    """ Building the diameters from the generated diameter models"""  

    all_models = {'M0': diametrize_model_generic,
            'M1': diametrize_model_generic,
            'M2': diametrize_model_generic,
            'M3': diametrize_model_generic,
            'M4': diametrize_model_generic
    }

    tqdm_1, tqdm_2 = utils.tqdm_disable(morphologies) #to have a single progression bar

    for mtype in tqdm(morphologies, disable = tqdm_1):
        for model in models:
            for neuron in tqdm(morphologies[mtype], disable = tqdm_2):

                all_models[model](neuron[0], models_params[mtype][model], neurite_types, extra_params[model])
                utils.save_neuron(neuron, model, new_morph_path) 
                plotting.plot_diameter_diff(neuron[1], morph_path, new_morph_path, model, neurite_types, folder='new_morphologies')

def diametrize_model_generic(neuron, params, neurite_types, extra_params):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    np.random.seed(extra_params['seed'])

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type == STR_TO_TYPES[neurite_type])

        for neurite in neurites: 

            max_bo = np.max(nm.get('section_term_branch_orders', neurite)) 
            if max_bo > extra_params['max_bo_fit']: #too few such branches, so we stick to the same parameters
                max_bo = extra_params['max_bo_fit']

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            k = 1
            while wrong_tips:
                wrong_tips = diametrize_tree(neurite, params, neurite_type, max_bo, trunk_diam_frac, taper = extra_params['taper'])
                n_tries += 1

                if n_tries > 10*k: #if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= 0.01
                    k+=1

                if n_tries > 900: #don't try to much, and complain
                    print('max tries attained with', neurite_type)
                    wrong_tips = False

def diametrize_tree(neurite, params, neurite_type, max_bo, trunk_diam_frac = 1., taper = 0):
        """ diametrize a single tree """

        trunk_diam = trunk_diam_frac*sample_distribution(params['trunk_diameter'][neurite_type], max_bo)

        wrong_tips = False

        status = {s.id: False for s in iter_sections(neurite)}
        active = [neurite.root_node]


        while active:

            for section in list(active):

                if section.is_root():
                    taper = -taper
                    init_diam = trunk_diam
                else:
                    taper = -taper
                    init_diam = get_diameters(section)[0]
                
                #sample a terminal diameter
                terminal_diam = sample_distribution(params['terminal_diameter'][neurite_type])
                
                #diametrize a section
                diametrize_section(section, init_diam, taper=taper,
                                             min_diam = terminal_diam, max_diam = trunk_diam)

                status[section.id] = True  # Tapering of section complete.
                active.remove(section)
                children = np.array(section.children)
                #if branching points has children, keep looping
                if len(children) > 1:
                    
                    reduc = 2. 
                    while reduc > 1.: #try until we get a reduction of diameter in the branching
                        sibling_ratio = sample_distribution(params['sibling_ratio'][neurite_type])
                        Rall_deviation = sample_distribution(params['Rall_deviation'][neurite_type])
                        reduc = morph_funcs.Rall_reduction_factor(Rall_deviation = Rall_deviation, siblings_ratio = sibling_ratio)

                    d0 = get_diameters(section)[-1]
                    d1 = reduc * d0
                    d2 = sibling_ratio * d1  

                    if d1 < terminal_diam:
                        d1 = terminal_diam

                    if d2 < terminal_diam:
                        d2 = terminal_diam

                    if len(children)>2:
                        print(len(children), 'children for this branch')

                    for i, ch in enumerate(children):
                        new_diam = d1 if i == 0 else d2
                        
                        utils.redefine_diameter_section(ch, 0, new_diam)

                        active.append(ch)

                #if we are at a tip, check tip diameters and stop
                else:
                    if get_diameters(section)[-1] < params['terminal_diameter'][neurite_type]['params']['min'] or get_diameters(section)[-1]  > params['terminal_diameter'][neurite_type]['params']['max']:
                        wrong_tips = True

        return wrong_tips


def diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=100.):
    '''Corrects the diameters of a section'''

    if initial_diam< min_diam:
        initial_diam = min_diam
    elif initial_diam> max_diam:
        initial_diam = max_diam
    diams = [initial_diam]
    if section.is_root():
        range_ = range(1, len(section.points))
    else:
        range_ = range(0, len(section.points) - 1)

    # lengths of each segments will be used for scaling of tapering
    lengths = [0] + utils.section_lengths(section).tolist()

    for i in range_:
        # Taper should be a negative number for decreasing diameters
        new_diam = diams[-1] + taper * lengths[i]

        if new_diam >= max_diam:
            diams.append(max_diam)
        elif new_diam <= min_diam:
            diams.append(min_diam)
        else:
            diams.append(new_diam)

    set_diameters(section, np.array(diams, dtype=np.float32))


