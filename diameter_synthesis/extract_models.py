import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom 

from tns.morphio_utils import STR_TO_TYPES, TYPE_TO_STR 

import diameter_synthesis.utils as utils 
import diameter_synthesis.plotting as plotting 

def sampling_model_test(morphologies, neurite_types, tqdm_disable = False):
    """ test for sampling models """
  
    #initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    Rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        Rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []

    #loop first over all morphologies (TODO: could be parallelized)
    for neuron in tqdm(morphologies, disable = tqdm_disable):
        #for each neurite in the neuron
        for i, neurite in enumerate(neuron.neurites):
            #for each type of neurite we consider
            for neurite_type in neurite_types:
                if neurite.type == STR_TO_TYPES[neurite_type]:

                    #compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += utils.sibling_ratios(neurite)
                    Rall_deviations[neurite_type] += utils.Rall_deviations(neurite)
                    terminal_diameters[neurite_type] += utils.terminal_diameters(neurite)
                    trunk_diameters[neurite_type] += utils.trunk_diameter(neurite)
  
    #do the fits of each morphological values
    sibling_ratio_models = {}
    Rall_deviation_models= {}
    terminal_diameters_models= {}
    trunk_diameters_models= {}
    for neurite_type in neurite_types:

        #sibling ratio 
        sibling_ratio_models[neurite_type] = {} 
        sibling_ratio_models[neurite_type]['distribution'] =  'expon_rev'
        sibling_ratio_models[neurite_type]['params'] =  utils.fit_distribution(sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'])
        
        #Rall deviation
        Rall_deviation_models[neurite_type] = {} 
        Rall_deviation_models[neurite_type]['distribution'] =  'exponnorm'
        Rall_deviation_models[neurite_type]['params'] =  utils.fit_distribution(Rall_deviations[neurite_type], Rall_deviation_models[neurite_type]['distribution'])
        
        #terminal diameters
        terminal_diameters_models[neurite_type] = {} 
        terminal_diameters_models[neurite_type]['distribution'] =  'gamma'
        terminal_diameters_models[neurite_type]['params'] =  utils.fit_distribution(terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'])

        #trunk diameters
        trunk_diameters_models[neurite_type] = {} 
        trunk_diameters_models[neurite_type]['distribution'] =  'gamma_sequence'
        min_sample_num = {'basal': 10, 'apical': 5}
        trunk_diameters_models[neurite_type]['params'] =  utils.fit_distribution(trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], min_sample_num = min_sample_num[neurite_type], floc = 0)
        trunk_diameters_models[neurite_type] = utils.update_params_fit_distribution(trunk_diameters[neurite_type], trunk_diameters_models[neurite_type])

    #collect all models in one dictionary
    all_models = {
    'sibling_ratio': sibling_ratio_models,
    'Rall_deviation': Rall_deviation_models, 
    'terminal_diameter': terminal_diameters_models, 
    'trunk_diameter':  trunk_diameters_models
    }

    all_data = {
    'sibling_ratio': sibling_ratios,
    'Rall_deviation': Rall_deviations, 
    'terminal_diameter': terminal_diameters, 
    'trunk_diameter':  trunk_diameters
    }
       
    return all_models, all_data


def build(models, morphologies, neurite_types, fig_folder = 'figures', ext = '.png'):
    """ Building the models in the list of models """  

    all_models = {'M0': sampling_model_test}
    
    #to have a single progression bar
    if len(morphologies) >1:
        tqdm_1 = False
        tqdm_2 = True
    else:
        tqdm_1 = True
        tqdm_2 = False

    #extrac tthe data and the models
    models_params = {} #dictionary of model parameters for each mtype
    models_data = {} #dictionary of model parameters for each mtype
    for mtype in tqdm(morphologies, disable = tqdm_1):
        models_params[mtype] = {}
        models_data[mtype] = {}
        for model in models:
            models_params[mtype][model], models_data[mtype][model] = all_models[model](morphologies[mtype], neurite_types, tqdm_2)

    #plot the distributions and fit of the data
    print('Plot the fits...')
    shutil.rmtree(fig_folder, ignore_errors = True)
    os.mkdir(fig_folder)
    for mtype in tqdm(morphologies): # for each mtypes
        os.mkdir(fig_folder +'/' + mtype)
        for model in models: #for each diameter model
            for fit_tpe in models_data[mtype][model]: #for each fit of the data we did 
                plotting.plot_distribution_fit(models_data[mtype][model][fit_tpe], models_params[mtype][model][fit_tpe], neurite_types, fig_name = fig_folder + '/' + mtype + '/' + model +'_' + fit_tpe, ext = ext)

    return models_params
