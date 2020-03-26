"""Construct diameter models from cells"""
import logging
from functools import partial

from neurom import NeuriteType

import diameter_synthesis.morph_functions as morph_funcs
from diameter_synthesis.distribution_fitting import fit_distribution
from diameter_synthesis.exception import DiameterSynthesisError

L = logging.getLogger(__name__)

STR_TO_NEUROM_TYPES = {
    "apical": NeuriteType.apical_dendrite,
    "basal": NeuriteType.basal_dendrite,
}


def _get_model_builder(config):
    """get the diameter model builder"""
    all_models = {}
    for model in config["models"]:
        if model == "generic":
            distribution_types = {}
            distribution_types["sibling_ratios"] = ["expon_rev", "asymmetry_threshold"]
            distribution_types["diameter_power_relation"] = [
                "exponnorm",
                "asymmetry_threshold",
            ]
            distribution_types["terminal_diameters"] = ["exponnorm", None]
            distribution_types["trunk_diameters"] = ["exponnorm", None]
            distribution_types["tapers"] = ["expon_rev", None]
        else:
            raise DiameterSynthesisError("model not understood {}".format(model))

        all_models[model] = partial(build_single_model, distribution_types)
    return all_models


def build(morphologies, config, with_data=False):
    """ Building a model from a set of morphologies and a config file"""
    all_models = _get_model_builder(config)
    models_params, models_data = all_models[config["models"][0]](
        morphologies,
        config["neurite_types"],
        config["extra_params"][config["models"][0]],
    )
    if with_data:
        return models_params, models_data
    return models_params


def build_single_model(sampling_model, morphologies, neurite_types, extra_params):
    """get diameter model from a set of dendrites"""
    all_data = extract_parameters(
        sampling_model, morphologies, neurite_types, extra_params
    )
    all_models = fit_all_models(all_data, sampling_model, extra_params, neurite_types)
    return all_models, all_data


def extract_parameters(sampling_model, morphologies, neurite_types, extra_params):
    """extract parameters from neurites"""

    all_data = {
        "sibling_ratios": {},
        "diameter_power_relation": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "tapers": {},
    }

    for neurite_type in neurite_types:
        all_data["sibling_ratios"][neurite_type] = []
        all_data["diameter_power_relation"][neurite_type] = []
        all_data["terminal_diameters"][neurite_type] = []
        all_data["trunk_diameters"][neurite_type] = []
        all_data["tapers"][neurite_type] = []

        for neuron in morphologies:
            for neurite in neuron.neurites:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:
                    all_data["sibling_ratios"][
                        neurite_type
                    ] += morph_funcs.compute_sibling_ratios(
                        neurite, attribute_name=sampling_model["sibling_ratios"][1]
                    )
                    all_data["diameter_power_relation"][
                        neurite_type
                    ] += morph_funcs.compute_diameter_power_relation(
                        neurite,
                        attribute_name=sampling_model["diameter_power_relation"][1],
                    )
                    all_data["terminal_diameters"][
                        neurite_type
                    ] += morph_funcs.terminal_diameters(
                        neurite,
                        threshold=extra_params["terminal_threshold"],
                        attribute_name=sampling_model["terminal_diameters"][1],
                    )
                    all_data["trunk_diameters"][
                        neurite_type
                    ] += morph_funcs.trunk_diameter(
                        neurite, attribute_name=sampling_model["trunk_diameters"][1]
                    )
                    all_data["tapers"][neurite_type] += morph_funcs.taper(
                        neurite,
                        params=extra_params["taper"],
                        attribute_name=sampling_model["tapers"][1],
                    )
    return all_data


def fit_all_models(all_data, sampling_model, extra_params, neurite_types):
    """fit the model parameters"""
    all_models = {
        "sibling_ratios": {},
        "diameter_power_relation": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "tapers": {},
    }

    for neurite_type in neurite_types:

        extra_params["neurite_type"] = neurite_type
        extra_params["name"] = "sibling_ratios"
        all_models["sibling_ratios"][neurite_type] = fit_model(
            sampling_model["sibling_ratios"],
            all_data["sibling_ratios"][neurite_type],
            extra_params,
        )

        extra_params["name"] = "diameter_power_relation"
        all_models["diameter_power_relation"][neurite_type] = fit_model(
            sampling_model["diameter_power_relation"],
            all_data["diameter_power_relation"][neurite_type],
            extra_params,
        )

        extra_params["name"] = "terminal_diameter"
        all_models["terminal_diameters"][neurite_type] = fit_model(
            sampling_model["terminal_diameters"],
            all_data["terminal_diameters"][neurite_type],
            extra_params,
        )

        extra_params["name"] = "trunk_diameters"
        all_models["trunk_diameters"][neurite_type] = fit_model(
            sampling_model["trunk_diameters"],
            all_data["trunk_diameters"][neurite_type],
            extra_params,
        )

        extra_params["name"] = "tapers"
        all_models["tapers"][neurite_type] = fit_model(
            sampling_model["tapers"], all_data["tapers"][neurite_type], extra_params,
        )

    return all_models


def fit_model(sampling_model, data, extra_params):
    """ fit a single parameter """
    output = {}
    if len(data) == 0:
        return output
    output["distribution"] = sampling_model[0]
    output["sequential"] = sampling_model[1]
    output["params"] = fit_distribution(
        data,
        sampling_model[0],
        attribute_name=sampling_model[1],
        extra_params=extra_params,
    )

    return output
