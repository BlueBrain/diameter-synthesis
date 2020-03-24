"""Build neurite diameters from a pre-generated model"""
import logging
import random
from collections import deque
from functools import partial

import numpy as np
from numpy.polynomial import polynomial

import diameter_synthesis.morph_functions as morph_funcs
import diameter_synthesis.utils as utils
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.exception import DiameterSynthesisError
from diameter_synthesis.types import STR_TO_TYPES

TRUNK_FRAC_DECREASE = 0.1
N_TRIES_BEFORE_REDUC = 5
L = logging.getLogger(__name__)


def _set_seed(params):
    """set numpy seed"""
    if "seed" in params:
        np.random.seed(params["seed"])


def _reset_caches():
    """reset the cached functions"""
    morph_funcs._sec_length.cache_clear()
    morph_funcs._partition_asymetry_length.cache_clear()
    morph_funcs._lengths_from_origin.cache_clear()
    morph_funcs._child_length.cache_clear()


def _get_neurites(neuron, neurite_type):
    """get iterator over neurites to diametrize"""
    return [
        list(neurite.iter())
        for neurite in neuron.iter()
        if neurite.is_root
        if neurite.type == STR_TO_TYPES[neurite_type]
    ]


def _get_sibling_ratio(
    params, neurite_type, asymetry_value=0.0, mode="generic", asymetry_threshold=0.3
):
    """return a sampled sibling ratio"""
    if mode == "generic":
        return sample_distribution(params["sibling_ratios"][neurite_type])
    if mode == "threshold":
        if asymetry_value > asymetry_threshold:
            return 0.0
        return sample_distribution(params["sibling_ratios"][neurite_type])
    raise DiameterSynthesisError("mode not understood {}".format(mode))


def _get_diameter_power_relation(
    params, neurite_type, asymetry_value=0.0, mode="generic", asymetry_threshold=0.3
):
    """return a sampled diameter power relation"""
    if mode == "generic":
        return sample_distribution(params["diameter_power_relation"][neurite_type])
    if mode == "threshold":
        if asymetry_value > asymetry_threshold:
            return 1.0
        return sample_distribution(params["diameter_power_relation"][neurite_type])
    if mode == "exact":
        return 1.0
    raise DiameterSynthesisError("mode not understood {}".format(mode))


def _get_trunk_diameter(params, neurite_type):
    """sample a trunk diameter"""
    return sample_distribution(params["trunk_diameters"][neurite_type])


def _get_terminal_diameter(params, neurite_type):
    """sample a terminal diameter"""
    return sample_distribution(params["terminal_diameters"][neurite_type])


def _get_taper(params, neurite_type, no_taper=False):
    """sample a taper rate"""
    if no_taper:
        return 0.0
    return sample_distribution(params["tapers"][neurite_type])


def _get_daughter_diameters(section, params, params_tree):
    """return daughter diameters from parent diameter"""

    reduction_factor = params_tree["reduction_factor_max"] + 1.0
    # try until we get a reduction of diameter in the branching
    while reduction_factor > params_tree["reduction_factor_max"]:
        asymetry_value = morph_funcs.get_additional_attribute(
            params["sibling_ratios"][params_tree["neurite_type"]]["sequential"],
            section=section,
        )

        if (
            params["sibling_ratios"][params_tree["neurite_type"]]["sequential"]
            == "asymmetry_threshold"
        ):
            asymetry_value /= params_tree["tot_length"]

        sibling_ratio = _get_sibling_ratio(
            params,
            params_tree["neurite_type"],
            asymetry_value=asymetry_value,
            mode=params_tree["mode_sibling"],
            asymetry_threshold=params_tree["asymetry_threshold"],
        )

        diameter_power_relation = _get_diameter_power_relation(
            params,
            params_tree["neurite_type"],
            asymetry_value=asymetry_value,
            mode=params_tree["mode_diameter_power_relation"],
            asymetry_threshold=params_tree["asymetry_threshold"],
        )

        reduction_factor = morph_funcs.diameter_power_relation_factor(
            diameter_power_relation, sibling_ratio
        )

    diam_0 = section.diameters[-1]
    terminal_diam = min(diam_0, params_tree["terminal_diam"])

    diam_1 = reduction_factor * diam_0
    diam_2 = sibling_ratio * diam_1

    diam_1 = max(diam_1, terminal_diam)
    diam_2 = max(diam_2, terminal_diam)

    diams = [diam_1] + (len(section.children) - 1) * [diam_2]
    if params_tree["with_asymmetry"]:
        child_sort = np.argsort(
            [morph_funcs._child_length(child) for child in section.children]
        )[::-1]
        return list(np.array(diams)[child_sort])
    return random.shuffle(diams)


def _diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=10.0):
    """Set the diameters of a section"""
    diams = polynomial.polyval(
        morph_funcs._lengths_from_origin(section), [initial_diam, taper]
    )
    section.diameters = np.clip(diams, min_diam, max_diam, out=diams)


def _diametrize_tree(neurite, params, params_tree):
    """Diametrize a single tree"""
    params_tree["tot_length"] = morph_funcs._get_total_length(neurite)
    max_diam = params["terminal_diameters"][params_tree["neurite_type"]]["params"][
        "max"
    ]
    wrong_tips = False
    active = deque([neurite[0]])
    while active:
        section = active.popleft()

        if section.is_root:
            init_diam = params_tree["trunk_diam"]
        else:
            init_diam = section.diameters[0]

        taper = _get_taper(
            params, params_tree["neurite_type"], no_taper=params_tree["no_taper"],
        )

        params_tree["terminal_diam"] = min(
            init_diam, _get_terminal_diameter(params, params_tree["neurite_type"])
        )

        _diametrize_section(
            section,
            init_diam,
            taper=taper,
            min_diam=params_tree["terminal_diam"],
            max_diam=params_tree["trunk_diam"],
        )

        if len(section.children) > 0:
            diams = _get_daughter_diameters(section, params, params_tree)

            for i, child in enumerate(section.children):
                utils._redefine_diameter_section(child, 0, diams[i])
                active.append(child)

        # if we are at a tip, check if tip diameters are small enough
        elif section.diameters[-1] > max_diam:
            wrong_tips = True

    return wrong_tips


def _diametrize_neuron(params_tree, neuron, params, neurite_types, extra_params):
    """Diametrize a morphio-neuron according to the model"""
    for neurite_type in neurite_types:
        params_tree["neurite_type"] = neurite_type
        params_tree["asymetry_threshold"] = extra_params["asymetry_threshold"][
            neurite_type
        ]

        for neurite in _get_neurites(neuron, neurite_type):
            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.0
            n_tries_step = 1
            while wrong_tips:
                trunk_diam = trunk_diam_frac * _get_trunk_diameter(params, neurite_type)

                if trunk_diam < 0.01:
                    trunk_diam = 1.0
                    L.warning("sampled trunk diameter < 0.01, so use 1 instead")

                params_tree["trunk_diam"] = trunk_diam
                wrong_tips = _diametrize_tree(neurite, params, params_tree)

                # if we can't get a good model, reduce the trunk diameter progressively
                if n_tries > N_TRIES_BEFORE_REDUC * n_tries_step:
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1
                if n_tries > extra_params["trunk_max_tries"]:
                    L.warning("max tries attained with %s", neurite_type)
                    wrong_tips = False
                n_tries += 1


def _select_model(model):
    """select a model to use from available ones"""
    if model == "generic":
        params_tree = {}
        params_tree["mode_sibling"] = "threshold"
        params_tree["mode_diameter_power_relation"] = "threshold"
        params_tree["with_asymmetry"] = True
        params_tree["no_taper"] = False
        params_tree["reduction_factor_max"] = 1.0
    else:
        raise DiameterSynthesisError("Unknown diameter model")

    return partial(_diametrize_neuron, params_tree)


def build(neuron, model, model_params, neurite_types, config):
    """Main function for building the diameters from the generated diameter models of a neuron"""
    _set_seed(config)
    _reset_caches()
    diameter_generator = _select_model(model)

    diameter_generator(neuron, model_params, neurite_types, config)
    if config["n_samples"] > 1:
        diameters = np.array(utils._get_all_diameters(neuron))
        for _ in range(config["n_samples"] - 1):
            diameter_generator(neuron, model_params, neurite_types, config)
            diameters += np.array(utils._get_all_diameters(neuron))
        utils._set_all_diameters(neuron, diameters / config["n_samples"])
