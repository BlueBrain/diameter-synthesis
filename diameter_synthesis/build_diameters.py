import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm

from neurom import COLS
from neurom.core import iter_sections

from utils import STR_TO_TYPES, TYPE_TO_STR

import diameter_synthesis.utils as utils 
import diameter_synthesis.morph_functions as morph_funcs
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.utils import get_diameters, set_diameters 
import diameter_synthesis.plotting as plotting

##################################
## Build diameters from a model ##
##################################

def build_diameters(models, models_params, morphologies, neurite_types, new_morph_path):
    """ Building the diameters from the generated diameter models"""  

    np.random.seed(10)
    all_models = {'M0': diametrize_model_0}

    tqdm_1, tqdm_2 = utils.tqdm_disable(morphologies) #to have a single progression bar

    for mtype in tqdm(morphologies, disable = tqdm_1):
        for model in models:
            for neuron in tqdm(morphologies[mtype], disable = tqdm_2):
                plotting.plot_neuron(neuron, 'original_neurons')

                all_models[model](neuron[0], models_params[mtype][model], neurite_types)

                plotting.plot_neuron(neuron, 'new_neurons')

                utils.save_neuron(neuron, model, new_morph_path) 

def diametrize_model_0(neuron, params, neurite_types = None):
    '''Corrects the diameters of a morphio-neuron according to the model.
       Starts from the root and moves towards the tips.
    '''

    for neurite_type in neurite_types:
        neurites = (neurite for neurite in neuron.neurites if neurite.type == STR_TO_TYPES[neurite_type])

        for neurite in neurites: 

            max_bo = np.max(nm.get('section_term_branch_orders', neurite)) 
            if max_bo > 8: #too few such branches, so we stick to the same parameters
                max_bo = 8 

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            k = 1
            while wrong_tips:
                wrong_tips = diametrize_tree(neurite, params, neurite_type, max_bo, trunk_diam_frac)
                n_tries += 1

                if n_tries > 100*k: #if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= 0.1
                    k+=1

                if n_tries > 900: #don't try to much, and complain
                    print('max tries attained with', neurite_type)
                    wrong_tips = False

def diametrize_tree(neurite, params, neurite_type, max_bo, trunk_diam_frac = 1.):
        """ diametrize a single tree """

        trunk_diam = trunk_diam_frac*sample_distribution(params['trunk_diameter'][neurite_type], max_bo)

        wrong_tips = False

        status = {s.id: False for s in iter_sections(neurite)}
        active = [neurite.root_node]

        while active:
            for section in list(active):

                if section.is_root():
                    taper = -0.000
                    init_diam = trunk_diam
                else:
                    taper = -0.000
                    init_diam = get_diameters(section)[0]
                
                #sample parameters from model
                params['terminal_diameter'][neurite_type]['params']['max'] = 0.8
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

                    d1 = reduc * get_diameters(section)[-1]
                    d2 = sibling_ratio * d1  
                    for i, ch in enumerate(children):
                        new_diam = d1 if i == 0 else d2
                        
                        #make sure the new diameter is not larger that terminal diameter
                        if new_diam < terminal_diam:
                            new_diam = terminal_diam

                        utils.redefine_diameter_section(ch, 0, new_diam)

                        active.append(ch)

                #if we are at a tip, check tip diameters and stop
                else:
                    if get_diameters(section)[-1] < params['terminal_diameter'][neurite_type]['params']['min'] or get_diameters(section)[-1]  > params['terminal_diameter'][neurite_type]['params']['max']:
                        wrong_tips = True

        return wrong_tips


def diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=100.):
    '''Corrects the diameters of a section'''

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



