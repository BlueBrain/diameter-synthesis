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
    print(neuron.name)
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

            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.
            k = 1
            while wrong_tips:

                #smple a trunk diameter
                params_tmp = params['trunk_diameter'][neurite_type]
                tpe = morph_funcs.sequential_single(params_tmp['sequential'], neurite = neurite)
                trunk_diam = trunk_diam_frac*sample_distribution(params_tmp, tpe[0])

                #try to diametrize the neurite
                wrong_tips = diametrize_tree(neurite, params, neurite_type, trunk_diam)
                
                #if we can't get a good model, reduce the trunk diameter progressively
                n_tries += 1
                if n_tries > 2*k: #if we keep failing, slighly reduce the trunk diams
                    trunk_diam_frac -= 0.10
                    k+=1

                #don't try to much and keep the latest try
                if n_tries > 90: 
                    print('max tries attained with', neurite_type)
                    wrong_tips = False

            #print(n_tries, 'tries for ', neurite_type, trunk_diam_frac)
def diametrize_tree(neurite, params, neurite_type, trunk_diam):
        """ diametrize a single tree """

        #initialise status variables
        wrong_tips = False # used to run until terminal diameters are thin enough
        active = [neurite.root_node] # list of sections to diametrize

        while active:
            for section in list(active):
                
                #set trunk diam if trunk, or first diam of the section otherwise
                if section.is_root():
                    init_diam = trunk_diam
                else:
                    init_diam = get_diameters(section)[0]
                
                #sample a terminal diameter
                params_tmp = params['terminal_diameter'][neurite_type]
                tpe = morph_funcs.sequential_single(params_tmp['sequential'], neurite = neurite)
                terminal_diam = sample_distribution(params_tmp, tpe[0])
               
                #remember the min and max diam for later
                min_diam = params_tmp['params']['min']
                max_diam = params_tmp['params']['max']
               
                #sample a taper rate
                params_tmp = params['taper'][neurite_type]
                tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                taper = sample_distribution(params_tmp, tpe[0]) #prevent positive tapers

                #diametrize a section
                diametrize_section(section, init_diam, taper=taper,
                                             min_diam = terminal_diam, max_diam = trunk_diam)

                #remove section from list of section to treat
                active.remove(section)

                #if branching points has children, keep looping
                if len(section.children) > 1:
                    
                    reduc_max = 1.2 # if set to larger than 1, we allow increase of diameters 
                    reduc = reduc_max + 1
                    while reduc > reduc_max: #try until we get a reduction of diameter in the branching

                        params_tmp = params['sibling_ratio'][neurite_type]
                        tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                   
                        ########### HARDCODED !!!! #############
                        if params_tmp['sequential'] == 'asymmetry':
                            if tpe < -0.8:
                                if neurite_type =='apical':
                                    tpe = -2.0
                                else:
                                    tpe = -1.
                            else:
                                tpe = -0.1
                        ##########################################

                        sibling_ratio = sample_distribution(params_tmp, tpe)

                        #sample a Rall deviation
                        params_tmp = params['Rall_deviation'][neurite_type]
                        tpe = morph_funcs.sequential_single(params_tmp['sequential'], section = section)
                        Rall_deviation = sample_distribution(params_tmp, tpe)

                        #compute the reduction factor
                        reduc = morph_funcs.Rall_reduction_factor(Rall_deviation = Rall_deviation, siblings_ratio = sibling_ratio)

                    #set new diameters to firt points of children
                    d0 = get_diameters(section)[-1]
                    if d0 < terminal_diam:
                        terminal_diam = d0

                    d1 = reduc * d0
                    d2 = sibling_ratio * d1
                    
                    #set minimum values if too small
                    if d1 < terminal_diam:
                        d1 = terminal_diam

                    if d2 < terminal_diam:
                        d2 = terminal_diam

                    #if the asymetry of partition is opposite to d1 > d2, switch diameters
                    #first compute the partition of each child
                    part = []
                    for child in range(len(section.children)):
                        part.append(float(sum(1 for _ in section.children[child].ipreorder())))

                    #sort them by larger first and create a list of diameters
                    child_sort = np.argsort(part)[::-1]
                    ds = [d1] + (len(section.children)-1)*[d2]
                    
                    #set diameters
                    for i, ch in enumerate(section.children):

                        if len(section.children) == 1:
                            print('single child')

                        utils.redefine_diameter_section(ch, 0, ds[child_sort[i]])
                        active.append(ch)

                #if we are at a tip, check tip diameters to restart if too large
                elif get_diameters(section)[-1]  > max_diam:
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

def replace_params(param_type, param_all_name='./model_params_all.json', model='M0'): 
    """replace the param file with the all neuron fit"""

    with open(param_all_name, 'r') as f:
        params_all = json.load(f)

    return params_all['all_types'][model][param_type]



