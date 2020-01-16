import os
import glob
import shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm
from neurom import COLS
from neurom.core import iter_sections

from diameter_synthesis.types import STR_TO_TYPES, TYPE_TO_STR

import diameter_synthesis.utils as utils
from diameter_synthesis.distribution_fitting import fit_distribution, update_params_fit_distribution
import diameter_synthesis.morph_functions as morph_funcs
import diameter_synthesis.plotting as plotting

##############################################
## Build a model from a set of morphologies ##
##############################################


def sampling_model_sibling_asymmetry_trunk(morphologies, neurite_types, extra_params, tqdm_disable=False):
    """ test for sampling models """

    sibling_sequential = 'asymmetry_threshold'
    rall_deviation_sequential = 'asymmetry_threshold'
    terminal_diameters_sequential = None
    trunk_diameters_sequential = 'max_branch' 
    tapers_sequential = None

    # initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    trunk_tapers = {}
    tapers = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []
        trunk_tapers[neurite_type] = []
        tapers[neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    i = 0
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += morph_funcs.sibling_ratios(neurite, seq=sibling_sequential)
                    rall_deviations[neurite_type] += morph_funcs.rall_deviations(neurite, seq=rall_deviation_sequential)
                    terminal_diameters[neurite_type] += morph_funcs.terminal_diameters(neurite, threshold=extra_params['terminal_threshold'], seq=terminal_diameters_sequential)
                    trunk_diameters[neurite_type] += morph_funcs.trunk_diameter(neurite, seq=trunk_diameters_sequential)
                    trunk_tapers[neurite_type] += morph_funcs.taper(neurite, only_first=True)
                    tapers[neurite_type] += morph_funcs.taper(neurite, params=extra_params['taper'], seq=tapers_sequential)

    # do the fits of each morphological values
    sibling_ratio_models = {}
    rall_deviation_models = {}
    terminal_diameters_models = {}
    trunk_diameters_models = {}
    trunk_tapers_models = {}
    tapers_models = {}
    for neurite_type in neurite_types:

        # sibling ratio
        sibling_ratio_models[neurite_type] = {}
        sibling_ratio_models[neurite_type]['distribution'] = 'expon_rev'
        sibling_ratio_models[neurite_type]['sequential'] = sibling_sequential
        sibling_ratio_models[neurite_type]['params'] = fit_distribution(sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'], seq=sibling_sequential, extra_params=extra_params, name ='sibling', threshold=extra_params['threshold'][neurite_type] )

        # Rall deviation
        rall_deviation_models[neurite_type] = {}
        rall_deviation_models[neurite_type]['distribution'] = 'exponnorm'
        rall_deviation_models[neurite_type]['sequential'] = rall_deviation_sequential
        rall_deviation_models[neurite_type]['params'] = fit_distribution(rall_deviations[neurite_type], rall_deviation_models[neurite_type]['distribution'], seq=rall_deviation_sequential, extra_params=extra_params, name ='Rall', threshold=extra_params['threshold'][neurite_type] )

        # terminal diameters
        terminal_diameters_models[neurite_type] = {}
        terminal_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        terminal_diameters_models[neurite_type]['sequential'] = terminal_diameters_sequential
        terminal_diameters_models[neurite_type]['params'] = fit_distribution(terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'], seq=terminal_diameters_sequential, extra_params=extra_params)

        # trunk diameters
        trunk_diameters_models[neurite_type] = {}
        trunk_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_diameters_models[neurite_type]['sequential'] = trunk_diameters_sequential
        trunk_diameters_models[neurite_type]['params'] = fit_distribution(trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], seq=trunk_diameters_sequential, extra_params=extra_params)

        # trunk taper 
        trunk_tapers_models[neurite_type] = {}
        trunk_tapers_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_tapers_models[neurite_type]['sequential'] = None
        trunk_tapers_models[neurite_type]['params'] = fit_distribution(trunk_tapers[neurite_type], trunk_tapers_models[neurite_type]['distribution'], seq=None, extra_params=extra_params)

        # taper
        tapers_models[neurite_type] = {}
        tapers_models[neurite_type]['distribution'] = 'exponnorm'
        tapers_models[neurite_type]['sequential'] = tapers_sequential
        tapers_models[neurite_type]['params'] = fit_distribution(
            tapers[neurite_type], tapers_models[neurite_type]['distribution'], seq=tapers_sequential, extra_params=extra_params)

    # collect all models in one dictionary
    all_models = {
        'sibling_ratio': sibling_ratio_models,
        'rall_deviation': rall_deviation_models,
        'terminal_diameter': terminal_diameters_models,
        'trunk_diameter': trunk_diameters_models,
        'trunk_taper': trunk_tapers_models,
        'taper': tapers_models
    }

    all_data = {
        'sibling_ratio': sibling_ratios,
        'rall_deviation': rall_deviations,
        'terminal_diameter': terminal_diameters,
        'trunk_diameter': trunk_diameters,
        'trunk_taper': trunk_tapers,
        'taper': tapers
    }

    return all_models, all_data

def sampling_model_sibling_asymmetry(morphologies, neurite_types, extra_params, tqdm_disable=False):
    """ test for sampling models """

    sibling_sequential = 'asymmetry_threshold'
    rall_deviation_sequential = 'asymmetry_threshold'
    terminal_diameters_sequential = None
    trunk_diameters_sequential = None #'max_branch' 
    tapers_sequential = None

    # initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    trunk_tapers = {}
    tapers = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []
        trunk_tapers[neurite_type] = []
        tapers[neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    i = 0
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += morph_funcs.sibling_ratios(neurite, seq=sibling_sequential)
                    rall_deviations[neurite_type] += morph_funcs.rall_deviations(neurite, seq=rall_deviation_sequential)
                    terminal_diameters[neurite_type] += morph_funcs.terminal_diameters(neurite, threshold=extra_params['terminal_threshold'], seq=terminal_diameters_sequential)
                    trunk_diameters[neurite_type] += morph_funcs.trunk_diameter(neurite, seq=trunk_diameters_sequential)
                    trunk_tapers[neurite_type] += morph_funcs.taper(neurite, only_first=True)
                    tapers[neurite_type] += morph_funcs.taper(neurite, params=extra_params['taper'], seq=tapers_sequential)

    # do the fits of each morphological values
    sibling_ratio_models = {}
    rall_deviation_models = {}
    terminal_diameters_models = {}
    trunk_diameters_models = {}
    trunk_tapers_models = {}
    tapers_models = {}
    for neurite_type in neurite_types:

        # sibling ratio
        sibling_ratio_models[neurite_type] = {}
        sibling_ratio_models[neurite_type]['distribution'] = 'expon_rev'
        sibling_ratio_models[neurite_type]['sequential'] = sibling_sequential
        sibling_ratio_models[neurite_type]['params'] = fit_distribution(sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'], seq=sibling_sequential, extra_params=extra_params, name ='sibling', threshold=extra_params['threshold'][neurite_type] )

        # Rall deviation
        rall_deviation_models[neurite_type] = {}
        rall_deviation_models[neurite_type]['distribution'] = 'exponnorm'
        rall_deviation_models[neurite_type]['sequential'] = rall_deviation_sequential
        rall_deviation_models[neurite_type]['params'] = fit_distribution(rall_deviations[neurite_type], rall_deviation_models[neurite_type]['distribution'], seq=rall_deviation_sequential, extra_params=extra_params, name ='Rall', threshold=extra_params['threshold'][neurite_type] )

        # terminal diameters
        terminal_diameters_models[neurite_type] = {}
        terminal_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        terminal_diameters_models[neurite_type]['sequential'] = terminal_diameters_sequential
        terminal_diameters_models[neurite_type]['params'] = fit_distribution(terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'], seq=terminal_diameters_sequential, extra_params=extra_params)

        # trunk diameters
        trunk_diameters_models[neurite_type] = {}
        trunk_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_diameters_models[neurite_type]['sequential'] = trunk_diameters_sequential
        trunk_diameters_models[neurite_type]['params'] = fit_distribution(trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], seq=trunk_diameters_sequential, extra_params=extra_params)

        # trunk taper 
        trunk_tapers_models[neurite_type] = {}
        trunk_tapers_models[neurite_type]['distribution'] = 'skewnorm'
        trunk_tapers_models[neurite_type]['sequential'] = None
        trunk_tapers_models[neurite_type]['params'] = fit_distribution(trunk_tapers[neurite_type], trunk_tapers_models[neurite_type]['distribution'], seq=None, extra_params=extra_params)

        # taper
        tapers_models[neurite_type] = {}
        tapers_models[neurite_type]['distribution'] = 'exponnorm'
        tapers_models[neurite_type]['sequential'] = tapers_sequential
        tapers_models[neurite_type]['params'] = fit_distribution(
            tapers[neurite_type], tapers_models[neurite_type]['distribution'], seq=tapers_sequential, extra_params=extra_params)

    # collect all models in one dictionary
    all_models = {
        'sibling_ratio': sibling_ratio_models,
        'rall_deviation': rall_deviation_models,
        'terminal_diameter': terminal_diameters_models,
        'trunk_diameter': trunk_diameters_models,
        'trunk_taper': trunk_tapers_models,
        'taper': tapers_models
    }

    all_data = {
        'sibling_ratio': sibling_ratios,
        'rall_deviation': rall_deviations,
        'terminal_diameter': terminal_diameters,
        'trunk_diameter': trunk_diameters,
        'trunk_taper': trunk_tapers,
        'taper': tapers
    }

    return all_models, all_data


def sampling_model_trunk_path(morphologies, neurite_types, extra_params, tqdm_disable=False):
    """ test for sampling models """

    sibling_sequential = None
    rall_deviation_sequential = None
    terminal_diameters_sequential = None
    trunk_diameters_sequential = 'max_path'
    tapers_sequential = None

    # initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    tapers = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []
        tapers[neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    i = 0
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += morph_funcs.sibling_ratios(neurite, seq=sibling_sequential)
                    rall_deviations[neurite_type] += morph_funcs.rall_deviations(neurite, seq=rall_deviation_sequential)
                    terminal_diameters[neurite_type] += morph_funcs.terminal_diameters(neurite, threshold=extra_params['terminal_threshold'], seq=terminal_diameters_sequential)
                    trunk_diameters[neurite_type] += morph_funcs.trunk_diameter(neurite, seq=trunk_diameters_sequential)
                    tapers[neurite_type] += morph_funcs.taper(neurite, params=extra_params['taper'], seq=tapers_sequential)

    # do the fits of each morphological values
    sibling_ratio_models = {}
    rall_deviation_models = {}
    terminal_diameters_models = {}
    trunk_diameters_models = {}
    tapers_models = {}
    for neurite_type in neurite_types:

        # sibling ratio
        sibling_ratio_models[neurite_type] = {}
        sibling_ratio_models[neurite_type]['distribution'] = 'expon_rev'
        sibling_ratio_models[neurite_type]['sequential'] = sibling_sequential
        sibling_ratio_models[neurite_type]['params'] = fit_distribution(
            sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'], seq=sibling_sequential, extra_params=extra_params)

        # Rall deviation
        rall_deviation_models[neurite_type] = {}
        rall_deviation_models[neurite_type]['distribution'] = 'exponnorm'
        rall_deviation_models[neurite_type]['sequential'] = rall_deviation_sequential
        rall_deviation_models[neurite_type]['params'] = fit_distribution(
            rall_deviations[neurite_type], rall_deviation_models[neurite_type]['distribution'], seq=rall_deviation_sequential, extra_params=extra_params)

        # terminal diameters
        terminal_diameters_models[neurite_type] = {}
        terminal_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        terminal_diameters_models[neurite_type]['sequential'] = terminal_diameters_sequential
        terminal_diameters_models[neurite_type]['params'] = fit_distribution(
            terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'], seq=terminal_diameters_sequential, extra_params=extra_params)

        # trunk diameters
        trunk_diameters_models[neurite_type] = {}
        trunk_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_diameters_models[neurite_type]['sequential'] = trunk_diameters_sequential
        trunk_diameters_models[neurite_type]['params'] = fit_distribution(
            trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], seq=trunk_diameters_sequential, extra_params=extra_params)

        # taper
        tapers_models[neurite_type] = {}
        tapers_models[neurite_type]['distribution'] = 'exponnorm'
        tapers_models[neurite_type]['sequential'] = tapers_sequential
        tapers_models[neurite_type]['params'] = fit_distribution(
            tapers[neurite_type], tapers_models[neurite_type]['distribution'], seq=tapers_sequential, extra_params=extra_params)

    # collect all models in one dictionary
    all_models = {
        'sibling_ratio': sibling_ratio_models,
        'rall_deviation': rall_deviation_models,
        'terminal_diameter': terminal_diameters_models,
        'trunk_diameter': trunk_diameters_models,
        'taper': tapers_models
    }

    all_data = {
        'sibling_ratio': sibling_ratios,
        'rall_deviation': rall_deviations,
        'terminal_diameter': terminal_diameters,
        'trunk_diameter': trunk_diameters,
        'taper': tapers
    }

    return all_models, all_data


def sampling_model_astrocyte(morphologies, neurite_types, extra_params, tqdm_disable=False):
    """ test for sampling models """

    sibling_sequential = None
    rall_deviation_sequential = None
    terminal_diameters_sequential = None
    trunk_diameters_sequential = None
    tapers_sequential = None

    # initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    tapers = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []
        tapers[neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    i = 0
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += morph_funcs.sibling_ratios(neurite, seq=sibling_sequential)
                    rall_deviations[neurite_type] += morph_funcs.rall_deviations(neurite, seq=rall_deviation_sequential, bounds = [0, 10])
                    terminal_diameters[neurite_type] += morph_funcs.min_diameter(neurite)
                    trunk_diameters[neurite_type] += morph_funcs.max_diameter(neurite)
                    tapers[neurite_type] += morph_funcs.taper(neurite, params=extra_params['taper'], seq=tapers_sequential)

    # do the fits of each morphological values
    sibling_ratio_models = {}
    rall_deviation_models = {}
    terminal_diameters_models = {}
    trunk_diameters_models = {}
    tapers_models = {}
    for neurite_type in neurite_types:

        # sibling ratio
        sibling_ratio_models[neurite_type] = {}
        sibling_ratio_models[neurite_type]['distribution'] = 'expon_rev'
        sibling_ratio_models[neurite_type]['sequential'] = sibling_sequential
        sibling_ratio_models[neurite_type]['params'] = fit_distribution(
            sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'], seq=sibling_sequential, extra_params=extra_params)

        # Rall deviation
        rall_deviation_models[neurite_type] = {}
        rall_deviation_models[neurite_type]['distribution'] = 'exponnorm'
        rall_deviation_models[neurite_type]['sequential'] = rall_deviation_sequential
        rall_deviation_models[neurite_type]['params'] = fit_distribution(
            rall_deviations[neurite_type], rall_deviation_models[neurite_type]['distribution'], seq=rall_deviation_sequential, extra_params=extra_params)

        # terminal diameters
        terminal_diameters_models[neurite_type] = {}
        terminal_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        terminal_diameters_models[neurite_type]['sequential'] = terminal_diameters_sequential
        terminal_diameters_models[neurite_type]['params'] = fit_distribution(
            terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'], seq=terminal_diameters_sequential, extra_params=extra_params)

        # trunk diameters
        trunk_diameters_models[neurite_type] = {}
        trunk_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_diameters_models[neurite_type]['sequential'] = trunk_diameters_sequential
        trunk_diameters_models[neurite_type]['params'] = fit_distribution(
            trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], seq=trunk_diameters_sequential, extra_params=extra_params)

        # taper
        tapers_models[neurite_type] = {}
        tapers_models[neurite_type]['distribution'] = 'exponnorm'
        tapers_models[neurite_type]['sequential'] = tapers_sequential
        tapers_models[neurite_type]['params'] = fit_distribution(
            tapers[neurite_type], tapers_models[neurite_type]['distribution'], seq=tapers_sequential, extra_params=extra_params)

    # collect all models in one dictionary
    all_models = {
        'sibling_ratio': sibling_ratio_models,
        'rall_deviation': rall_deviation_models,
        'terminal_diameter': terminal_diameters_models,
        'trunk_diameter': trunk_diameters_models,
        'taper': tapers_models
    }

    all_data = {
        'sibling_ratio': sibling_ratios,
        'rall_deviation': rall_deviations,
        'terminal_diameter': terminal_diameters,
        'trunk_diameter': trunk_diameters,
        'taper': tapers
    }

    return all_models, all_data


def sampling_model_generic(morphologies, neurite_types, extra_params, tqdm_disable=False):
    """ test for sampling models """

    sibling_sequential = None
    rall_deviation_sequential = None
    terminal_diameters_sequential = None
    trunk_diameters_sequential = None
    tapers_sequential = None

    # initialise dictionaries for collecting morphological quantities
    sibling_ratios = {}
    rall_deviations = {}
    terminal_diameters = {}
    trunk_diameters = {}
    tapers = {}
    for neurite_type in neurite_types:
        sibling_ratios[neurite_type] = []
        rall_deviations[neurite_type] = []
        terminal_diameters[neurite_type] = []
        trunk_diameters[neurite_type] = []
        tapers[neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    i = 0
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    sibling_ratios[neurite_type] += morph_funcs.sibling_ratios(neurite, seq=sibling_sequential)
                    rall_deviations[neurite_type] += morph_funcs.rall_deviations(neurite, seq=rall_deviation_sequential)
                    terminal_diameters[neurite_type] += morph_funcs.terminal_diameters(neurite, threshold=extra_params['terminal_threshold'], seq=terminal_diameters_sequential)
                    trunk_diameters[neurite_type] += morph_funcs.trunk_diameter(neurite, seq=trunk_diameters_sequential)
                    tapers[neurite_type] += morph_funcs.taper(neurite, params=extra_params['taper'], seq=tapers_sequential)

    # do the fits of each morphological values
    sibling_ratio_models = {}
    rall_deviation_models = {}
    terminal_diameters_models = {}
    trunk_diameters_models = {}
    tapers_models = {}
    for neurite_type in neurite_types:

        # sibling ratio
        sibling_ratio_models[neurite_type] = {}
        sibling_ratio_models[neurite_type]['distribution'] = 'expon_rev'
        sibling_ratio_models[neurite_type]['sequential'] = sibling_sequential
        sibling_ratio_models[neurite_type]['params'] = fit_distribution(
            sibling_ratios[neurite_type], sibling_ratio_models[neurite_type]['distribution'], seq=sibling_sequential, extra_params=extra_params)

        # Rall deviation
        rall_deviation_models[neurite_type] = {}
        rall_deviation_models[neurite_type]['distribution'] = 'exponnorm'
        rall_deviation_models[neurite_type]['sequential'] = rall_deviation_sequential
        rall_deviation_models[neurite_type]['params'] = fit_distribution(
            rall_deviations[neurite_type], rall_deviation_models[neurite_type]['distribution'], seq=rall_deviation_sequential, extra_params=extra_params)

        # terminal diameters
        terminal_diameters_models[neurite_type] = {}
        terminal_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        terminal_diameters_models[neurite_type]['sequential'] = terminal_diameters_sequential
        terminal_diameters_models[neurite_type]['params'] = fit_distribution(
            terminal_diameters[neurite_type], terminal_diameters_models[neurite_type]['distribution'], seq=terminal_diameters_sequential, extra_params=extra_params)

        # trunk diameters
        trunk_diameters_models[neurite_type] = {}
        trunk_diameters_models[neurite_type]['distribution'] = 'exponnorm'
        trunk_diameters_models[neurite_type]['sequential'] = trunk_diameters_sequential
        trunk_diameters_models[neurite_type]['params'] = fit_distribution(
            trunk_diameters[neurite_type], trunk_diameters_models[neurite_type]['distribution'], seq=trunk_diameters_sequential, extra_params=extra_params)

        # taper
        tapers_models[neurite_type] = {}
        tapers_models[neurite_type]['distribution'] = 'exponnorm'
        tapers_models[neurite_type]['sequential'] = tapers_sequential
        tapers_models[neurite_type]['params'] = fit_distribution(
            tapers[neurite_type], tapers_models[neurite_type]['distribution'], seq=tapers_sequential, extra_params=extra_params)

    # collect all models in one dictionary
    all_models = {
        'sibling_ratio': sibling_ratio_models,
        'rall_deviation': rall_deviation_models,
        'terminal_diameter': terminal_diameters_models,
        'trunk_diameter': trunk_diameters_models,
        'taper': tapers_models
    }

    all_data = {
        'sibling_ratio': sibling_ratios,
        'rall_deviation': rall_deviations,
        'terminal_diameter': terminal_diameters,
        'trunk_diameter': trunk_diameters,
        'taper': tapers
    }

    return all_models, all_data


def build_models(models, morphologies, neurite_types, extra_params, fig_folder='figures', ext='.png', plot=True):
    """ Building the models in the list of models """

    all_models = {}
    for model in models:
        if model == 'M0':
            all_models[model] = sampling_model_generic
        elif model == 'M1':
            all_models[model] = sampling_model_sibling_asymmetry
        elif model == 'M2':
            all_models[model] = sampling_model_sibling_asymmetry_trunk
        elif model == 'M3':
            all_models[model] = sampling_model_astrocyte

    tqdm_1, tqdm_2 = utils.tqdm_disable(morphologies)  # to have a single progression bar

    # extract the data and the models
    models_params = {}  # dictionary of model parameters for each mtype
    models_data = {}  # dictionary of model parameters for each mtype
    for mtype in tqdm(morphologies, disable=tqdm_1):
        models_params[mtype] = {}
        models_data[mtype] = {}
        for model in models:
            models_params[mtype][model], models_data[mtype][model] = all_models[model](morphologies[mtype], neurite_types, extra_params[model], tqdm_2)

    # plot the distributions and fit of the data
    if plot:
        print('Plot the fits...')
        #shutil.rmtree(fig_folder, ignore_errors=True)
        if not os.path.isdir(fig_folder):
            os.mkdir(fig_folder)
        for mtype in tqdm(morphologies):  # for each mtypes
            if not os.path.isdir(fig_folder + '/' + mtype):
                os.mkdir(fig_folder + '/' + mtype)
            for model in models:  # for each diameter model
                for fit_tpe in models_data[mtype][model]:  # for each fit of the data we did
                    plotting.plot_distribution_fit(models_data[mtype][model][fit_tpe], models_params[mtype][model][fit_tpe], neurite_types, fig_name=fig_folder + '/' + mtype + '/' + model + '_' + fit_tpe, ext=ext)

    return models_params
