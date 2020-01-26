""" Construct diameter models from cells """
import os
from functools import partial

from tqdm import tqdm

import diameter_synthesis.morph_functions as morph_funcs
import diameter_synthesis.plotting as plotting
import diameter_synthesis.utils as utils
from diameter_synthesis.distribution_fitting import fit_distribution
from diameter_synthesis.types import STR_TO_TYPES

############################################
# Build a model from a set of morphologies #
############################################


def get_models(config):
    """ get the model parameters """

    all_models = {}
    for model in config["models"]:
        if model == "M1":
            distribution_types = {}
            distribution_types["sibling_ratios"] = ["expon_rev", "asymmetry_threshold"]
            distribution_types["rall_deviations"] = ["exponnorm", "asymmetry_threshold"]
            distribution_types["terminal_diameters"] = ["exponnorm", None]
            distribution_types["trunk_diameters"] = ["exponnorm", None]
            distribution_types["trunk_tapers"] = ["skewnorm", None]
            distribution_types["tapers"] = ["exponnorm", None]

            all_models[model] = partial(build_single_model, distribution_types)

    return all_models


def build_models(morphologies, config, single_model=False):
    """ Building the models in the list of models """

    all_models = get_models(config)

    # extract the data and the models
    models_params = {}  # dictionary of model parameters for each mtype
    models_data = {}  # dictionary of model parameters for each mtype
    for mtype in tqdm(morphologies, disable=utils.tqdm_disable(morphologies)[0]):
        models_params[mtype] = {}
        models_data[mtype] = {}
        if single_model:
            models_params[mtype], models_data[mtype] = all_models[config["models"][0]](
                morphologies[mtype],
                config["neurite_types"],
                config["extra_params"][config["models"][0]],
                utils.tqdm_disable(morphologies)[1],
            )
        else:
            for model in config["models"]:
                models_params[mtype][model], models_data[mtype][model] = all_models[
                    model
                ](
                    morphologies[mtype],
                    config["neurite_types"],
                    config["extra_params"][model],
                    utils.tqdm_disable(morphologies)[2],
                )

    # plot the distributions and fit of the data
    if config["plot"]:
        plot_models(morphologies, config, models_params, models_data, single_model=single_model)

    return models_params


def plot_models(morphologies, config, models_params, models_data, single_model=False):
    """plot the models"""
    try:
        print("Plot the fits...")

        if not os.path.isdir(config["fig_folder"]):
            os.mkdir(config["fig_folder"])

        for mtype in tqdm(morphologies):  # for each mtypes
            if not os.path.isdir(os.path.join(config["fig_folder"], mtype)):
                os.mkdir(os.path.join(config["fig_folder"], mtype))

            for model in config["models"]:  # for each diameter model
                if not single_model:
                    fit_tpes = models_data[mtype][model]
                    model_data = models_data[mtype][model]
                    model_param = models_params[mtype][model]
                    mtype = os.path.join(mtype, model + '_')
                else:
                    fit_tpes = models_data[mtype]
                    model_data = models_data[mtype]
                    model_param = models_params[mtype]

                # for each fit of the data we did
                for fit_tpe in fit_tpes:
                    fig_name = os.path.join(config["fig_folder"], mtype) + fit_tpe
                    plotting.plot_distribution_fit(
                        model_data[fit_tpe],
                        model_param[fit_tpe],
                        config["neurite_types"],
                        fig_name=fig_name,
                        ext=config["ext"],
                    )
    except Exception as exc:  # pylint: disable=broad-except
        print('Could not plot models because of', exc)


def build_single_model(
        sampling_model, morphologies, neurite_types, extra_params, tqdm_disable=False
):
    """ get diameter model from a set of dendrites """
    all_data = extract_parameters(sampling_model,
                                  morphologies,
                                  neurite_types,
                                  extra_params,
                                  tqdm_disable)

    all_models = fit_all_models(all_data,
                                sampling_model,
                                extra_params,
                                neurite_types)

    return all_models, all_data


def extract_parameters(
        sampling_model, morphologies, neurite_types, extra_params, tqdm_disable=False
):
    """ extract parameters from neurites """

    all_data = {
        "sibling_ratios": {},
        "rall_deviations": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "trunk_tapers": {},
        "tapers": {},
    }

    # initialise dictionaries for collecting morphological quantities
    for neurite_type in neurite_types:
        all_data['sibling_ratios'][neurite_type] = []
        all_data['rall_deviations'][neurite_type] = []
        all_data['terminal_diameters'][neurite_type] = []
        all_data['trunk_diameters'][neurite_type] = []
        all_data['trunk_tapers'][neurite_type] = []
        all_data['tapers'][neurite_type] = []

    # loop first over all morphologies (TODO: could be parallelized)
    for neuron in tqdm(morphologies, disable=tqdm_disable):
        # for each neurite in the neuron
        for neurite in neuron.neurites:
            # for each type of neurite we consider
            for neurite_type in neurite_types:

                if neurite.type == STR_TO_TYPES[neurite_type]:

                    # compute here all the morphological values from the neurite
                    all_data['sibling_ratios'][neurite_type] += morph_funcs.sibling_ratios(
                        neurite, seq=sampling_model["sibling_ratios"][1]
                    )
                    all_data['rall_deviations'][neurite_type] += morph_funcs.rall_deviations(
                        neurite, seq=sampling_model["rall_deviations"][1]
                    )
                    all_data['terminal_diameters'][neurite_type] += morph_funcs.terminal_diameters(
                        neurite,
                        threshold=extra_params["terminal_threshold"],
                        seq=sampling_model["terminal_diameters"][1],
                    )
                    all_data['trunk_diameters'][neurite_type] += morph_funcs.trunk_diameter(
                        neurite, seq=sampling_model["trunk_diameters"][1]
                    )
                    all_data['trunk_tapers'][neurite_type] += morph_funcs.taper(
                        neurite, only_first=True
                    )
                    all_data['tapers'][neurite_type] += morph_funcs.taper(
                        neurite,
                        params=extra_params["taper"],
                        seq=sampling_model["tapers"][1],
                    )
    return all_data


def fit_all_models(all_data, sampling_model, extra_params, neurite_types):
    """ fit the model parameters """

    all_models = {
        "sibling_ratios": {},
        "rall_deviations": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "trunk_tapers": {},
        "tapers": {},
    }

    # do the fits of each morphological values
    for neurite_type in neurite_types:

        extra_params["neurite_type"] = neurite_type

        # sibling ratio
        all_models['sibling_ratios'][neurite_type] = fit_model(
            sampling_model["sibling_ratios"],
            all_data['sibling_ratios'][neurite_type],
            extra_params,
        )

        # Rall deviation
        all_models['rall_deviations'][neurite_type] = fit_model(
            sampling_model["rall_deviations"],
            all_data['rall_deviations'][neurite_type],
            extra_params,
        )

        # terminal diameters
        all_models['terminal_diameters'][neurite_type] = fit_model(
            sampling_model["terminal_diameters"],
            all_data['terminal_diameters'][neurite_type],
            extra_params,
        )

        # trunk diameters
        all_models['trunk_diameters'][neurite_type] = fit_model(
            sampling_model["trunk_diameters"],
            all_data['trunk_diameters'][neurite_type],
            extra_params,
        )

        # trunk taper
        all_models['trunk_tapers'][neurite_type] = fit_model(
            sampling_model["trunk_tapers"],
            all_data['trunk_tapers'][neurite_type],
            extra_params,
        )

        # taper
        all_models['tapers'][neurite_type] = fit_model(
            sampling_model["tapers"],
            all_data['tapers'][neurite_type],
            extra_params,
        )

    return all_models


def fit_model(sampling_model, data, extra_params):
    """ fit a single parameter """
    output = {}
    output["distribution"] = sampling_model[0]
    output["sequential"] = sampling_model[1]
    output["params"] = fit_distribution(
        data, sampling_model[0], seq=sampling_model[1], extra_params=extra_params
    )

    return output
