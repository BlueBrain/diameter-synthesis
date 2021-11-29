"""Construct diameter models from cells."""

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
    "axon": NeuriteType.axon,
}


def _get_model_builder(config):
    """Get the diameter model builder.

    Args:
        config (dict): configuration dictionary.

    Returns:
        dict: dictionary of models.
    """
    all_models = {}
    for model in config["models"]:
        if model == "generic":
            distribution_types = {
                "sibling_ratios": ["expon_rev", None],
                "diameter_power_relation": ["exponnorm", None],
                "terminal_diameters": ["exponnorm", None],
                "trunk_diameters": ["exponnorm", None],
                "trunk_diameters_method": "last",
                "tapers": ["expon_rev", None],
            }
        elif model == "astrocyte":
            distribution_types = {
                "sibling_ratios": ["expon_rev", None],
                "diameter_power_relation": ["exponnorm", None],
                "terminal_diameters": ["exponnorm", None],
                "trunk_diameters_method": "mean",
                "trunk_diameters": ["exponnorm", None],
                "tapers": ["exponnorm", None],
            }
        else:
            raise DiameterSynthesisError(f"model not understood {model}")

        all_models[model] = partial(build_single_model, distribution_types)
    return all_models


def build(morphologies, config, with_data=False):
    """Build a model from a set of morphologies and a config file.

    Args:
        morphologies (dict): dictionary of morphologies.
        config (dict): configuration dictionary.
        with_data (bool): if set to True, the model data are also returned.

    Returns:
        dict: parameter of the models (and possibly the data extracted if `with_data==True`).
    """
    all_models = _get_model_builder(config)
    models_params, models_data = all_models[config["models"][0]](morphologies, config)
    if with_data:
        return models_params, models_data
    return models_params


def build_single_model(sampling_model, morphologies, config):
    """Build a single model.

    Args:
        sampling_model (dict): parameters for model building.
        morphologies (dict): dictionary of morphologies.
        config (dict): configuration dictionary.

    Returns:
        dict: parameter of the models data extracted.
    """
    all_data = extract_parameters(sampling_model, morphologies, config)
    all_models = fit_all_models(all_data, sampling_model, config)
    return all_models, all_data


def extract_parameters(sampling_model, morphologies, config):
    """Extract parameters from neurites to then use for model building.

    Args:
        sampling_model (dict): parameters for model building.
        morphologies (dict): dictionary of morphologies.
        config (dict): configuration dictionary.

    Returns:
        dict: data extracted.
    """
    all_data = {
        "sibling_ratios": {},
        "diameter_power_relation": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "tapers": {},
    }

    for neurite_type in config["neurite_types"]:
        all_data["sibling_ratios"][neurite_type] = []
        all_data["diameter_power_relation"][neurite_type] = []
        all_data["terminal_diameters"][neurite_type] = []
        all_data["trunk_diameters"][neurite_type] = []
        all_data["tapers"][neurite_type] = []

        for neuron in morphologies:
            for neurite in neuron.neurites:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:
                    all_data["sibling_ratios"][neurite_type] += morph_funcs.compute_sibling_ratios(
                        neurite, attribute_name=sampling_model["sibling_ratios"][1]
                    )
                    all_data["diameter_power_relation"][
                        neurite_type
                    ] += morph_funcs.compute_diameter_power_relation(
                        neurite,
                        attribute_name=sampling_model["diameter_power_relation"][1],
                    )
                    all_data["terminal_diameters"][neurite_type] += morph_funcs.terminal_diameters(
                        neurite,
                        threshold=config["terminal_threshold"],
                        attribute_name=sampling_model["terminal_diameters"][1],
                    )
                    all_data["trunk_diameters"][neurite_type] += morph_funcs.trunk_diameter(
                        neurite,
                        attribute_name=sampling_model["trunk_diameters"][1],
                        method=sampling_model["trunk_diameters_method"],
                    )
                    all_data["tapers"][neurite_type] += morph_funcs.taper(
                        neurite,
                        params=config["taper"],
                        attribute_name=sampling_model["tapers"][1],
                    )
    return all_data


def fit_all_models(all_data, sampling_model, config):
    """Fit the parameters to get models.

    Args:
        all_data (dict): parameters to fit.
        sampling_model (dict): parameters for model building.
        config (dict): configuration dictionary.

    Returns:
        dict: models.
    """
    all_models = {
        "sibling_ratios": {},
        "diameter_power_relation": {},
        "terminal_diameters": {},
        "trunk_diameters": {},
        "tapers": {},
    }

    extra_params = config.copy()
    for neurite_type in config["neurite_types"]:

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
            sampling_model["tapers"], all_data["tapers"][neurite_type], extra_params
        )

    return all_models


def fit_model(sampling_model, data, extra_params):
    """Fit a single parameter.

    Args:
        sampling_model (dict): model parameters.
        data (array): parameter to fit.
        extra_params (dict): additional parameters.

    Returns:
        dict: fit parameters.
    """
    if len(data) == 0:
        return {}
    return {
        "distribution": sampling_model[0],
        "sequential": sampling_model[1],
        "params": fit_distribution(
            data,
            sampling_model[0],
            attribute_name=sampling_model[1],
            extra_params=extra_params,
        ),
    }
