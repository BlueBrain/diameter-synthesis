"""Build neurite diameters from a pre-generated model. TO REVIEW!."""
import logging
import random
from collections import deque
from functools import partial

import numpy as np
from numpy.polynomial import polynomial

from morphio import SectionType

import diameter_synthesis.morph_functions as morph_funcs
import diameter_synthesis.utils as utils
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.exception import DiameterSynthesisError

TRUNK_FRAC_DECREASE = 0.1
N_TRIES_BEFORE_REDUC = 5
L = logging.getLogger(__name__)

STR_TO_TYPES = {
    "basal": SectionType.basal_dendrite,
    "apical": SectionType.apical_dendrite,
}


def _reset_caches():
    """Reset the cached functions."""
    morph_funcs.sec_length.cache_clear()
    morph_funcs.partition_asymetry_length.cache_clear()
    morph_funcs.lengths_from_origin.cache_clear()
    morph_funcs.n_children_downstream.cache_clear()


def _get_neurites(neuron, neurite_type):
    """Get a list of neurites to diametrize.

    Args:
        neuron (morphio.mu.Morphology): a neuron
        neurite_type (str): the neurite type to consider

    Returns:
        list: list of neurites to consider
    """
    return [
        list(neurite.iter())
        for neurite in neuron.root_sections
        if neurite.type == STR_TO_TYPES[neurite_type]
    ]


def _sample_sibling_ratio(
    params, neurite_type, asymetry_value=0.0, mode="generic", asymetry_threshold=0.3
):
    """Sample a sibling ratio from distribution.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider
        asymetry_value (float): asymetry of current branching point
        mode (str): to use or not the asymetry_threshold
        asymetry_threshold (float): asymetry threshold

    Returns:
        float: sibling ratio
    """
    if mode == "generic":
        return sample_distribution(params["sibling_ratios"][neurite_type])
    if mode == "threshold":
        if asymetry_value > asymetry_threshold:
            return 0.0
        return sample_distribution(params["sibling_ratios"][neurite_type])
    raise DiameterSynthesisError("mode not understood {}".format(mode))


def _sample_diameter_power_relation(
    params, neurite_type, asymetry_value=0.0, mode="generic", asymetry_threshold=0.3
):
    """Sample a diameter power relation from distribution.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider
        asymetry_value (float): asymetry of current branching point
        mode (str): to use or not the asymetry_threshold
        asymetry_threshold (float): asymetry threshold

    Returns:
        float: diameter power relation
    """
    if mode == "generic":
        return sample_distribution(params["diameter_power_relation"][neurite_type])
    if mode == "threshold":
        if asymetry_value > asymetry_threshold:
            return 1.0
        return sample_distribution(params["diameter_power_relation"][neurite_type])
    if mode == "exact":
        return 1.0
    raise DiameterSynthesisError("mode not understood {}".format(mode))


def _sample_trunk_diameter(params, neurite_type):
    """Sample a trunk diameter from distribution.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider

    Returns:
        float: trunk diameter
    """
    return sample_distribution(params["trunk_diameters"][neurite_type])


def _sample_terminal_diameter(params, neurite_type):
    """Sample a terminal diameter.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider

    Returns:
        float: terminal diameter
    """
    return sample_distribution(params["terminal_diameters"][neurite_type])


def _sample_taper(params, neurite_type):
    """Sample a taper rate from distributions.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider

    Returns:
        float: taper rate
    """
    return sample_distribution(params["tapers"][neurite_type])


def _sample_daughter_diameters(section, params, params_tree):
    """Compute the daughter diameters of the current section.

    Args:
        section (morphio  section): section to consider
        params (dict): model parameters
        params_tree (dict): specific parameters of the current tree

    Returns:
       list: list of daughter diameters
    """
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

        sibling_ratio = _sample_sibling_ratio(
            params,
            params_tree["neurite_type"],
            asymetry_value=asymetry_value,
            mode=params_tree["mode_sibling"],
            asymetry_threshold=params_tree["asymetry_threshold"],
        )

        diameter_power_relation = _sample_diameter_power_relation(
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
        # returns child diameters sorted by child length (major/secondary for apical tree)
        child_sort = np.argsort(
            [morph_funcs.n_children_downstream(child) for child in section.children]
        )[::-1]
        return list(np.array(diams)[child_sort])
    random.shuffle(diams)
    return diams


def _diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=10.0):
    """Diameterize a section.

    Args:
        section (morphio section): current section
        initial_diam (float): initial diameter
        taper (float): taper rate
        min_diam (flaot): minimum diameter
        max_diam (float): maximum diameter
    """
    diams = polynomial.polyval(
        morph_funcs.lengths_from_origin(section), [initial_diam, taper]
    )
    section.diameters = np.clip(diams, min_diam, max_diam)


def _diametrize_tree(neurite, params, params_tree):
    """Diametrize a tree, or neurite.

    Args:
        neurite (morphio neurite): current neurite
        params (dict): model parameters
        params_tree (dict): specific parameters of the current tree

    Returns:
        bool: True is all terminal diameters are small enough, False otherwise
    """
    params_tree["tot_length"] = morph_funcs.get_total_length(neurite)
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

        taper = _sample_taper(params, params_tree["neurite_type"])

        params_tree["terminal_diam"] = min(
            init_diam, _sample_terminal_diameter(params, params_tree["neurite_type"])
        )

        _diametrize_section(
            section,
            init_diam,
            taper=taper,
            min_diam=params_tree["terminal_diam"],
            max_diam=params_tree["trunk_diam"],
        )

        if len(section.children) > 0:
            diams = _sample_daughter_diameters(section, params, params_tree)

            for i, child in enumerate(section.children):
                utils.redefine_diameter_section(child, 0, diams[i])
                active.append(child)

        # if we are at a tip, check if tip diameters are small enough
        elif section.diameters[-1] > max_diam:
            wrong_tips = True

    return wrong_tips


def _diametrize_neuron(params_tree, neuron, params, neurite_types, config):
    """Diametrize a neuron.

    Args:
        params_tree (dict): specific parameters of the current tree
        neuron (morphio.mut.Morphology): neuron to diametrize
        params (dict): model parameters
        neurite_type (str): the neurite type to consider
        config (dict): general configuration parameters
    """
    for neurite_type in neurite_types:
        params_tree["neurite_type"] = neurite_type
        params_tree["asymetry_threshold"] = config["asymetry_threshold"][neurite_type]

        for neurite in _get_neurites(neuron, neurite_type):
            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.0
            n_tries_step = 1
            while wrong_tips:
                trunk_diam = trunk_diam_frac * _sample_trunk_diameter(
                    params, neurite_type
                )

                if trunk_diam < 0.01:
                    trunk_diam = 1.0
                    L.warning("sampled trunk diameter < 0.01, so use 1 instead")

                params_tree["trunk_diam"] = trunk_diam
                wrong_tips = _diametrize_tree(neurite, params, params_tree)

                # if we can't get a good model, reduce the trunk diameter progressively
                if n_tries > N_TRIES_BEFORE_REDUC * n_tries_step:
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1
                if n_tries > config["trunk_max_tries"]:
                    L.warning("max tries attained with %s", neurite_type)
                    wrong_tips = False
                n_tries += 1


def _select_model(model):
    """Select a diametrized model to use.

    Args:
        model (str): model name

    Returns:
        function: diamtrizer with specific params_tree
    """
    if model == "generic":
        params_tree = {}
        params_tree["mode_sibling"] = "threshold"
        params_tree["mode_diameter_power_relation"] = "threshold"
        params_tree["with_asymmetry"] = True
        params_tree["reduction_factor_max"] = 1.0
    else:
        raise DiameterSynthesisError("Unknown diameter model")

    return partial(_diametrize_neuron, params_tree)


def build(neuron, model_params, neurite_types, config):
    """Builder function for generating diameters of a neuron from the a diameter models.

    Args:
        neuron (morphio.mut.Morphology): neuron to diametrize
        model_params (dict): model parameters
        neurite_type (str): the neurite type to consider
        config (dict): general configuration parameters
    """
    if "seed" in config:
        np.random.seed(config["seed"])
    _reset_caches()

    if len(config["models"]) > 1:
        L.warning("Several models provided, we will only use the first")
    diameter_generator = _select_model(config["models"][0])

    diameter_generator(neuron, model_params, neurite_types, config)
    if config["n_samples"] > 1:
        diameters = np.array(utils.get_all_diameters(neuron))
        for _ in range(config["n_samples"] - 1):
            diameter_generator(neuron, model_params, neurite_types, config)
            diameters += np.array(utils.get_all_diameters(neuron))
        utils.set_all_diameters(neuron, diameters / config["n_samples"])
