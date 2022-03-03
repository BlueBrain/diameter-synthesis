"""Build neurite diameters from a pre-generated model."""

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
from collections import deque
from copy import copy
from copy import deepcopy
from functools import partial

import numpy as np
from morphio import IterType
from morphio import SectionType
from numpy.polynomial import polynomial

from diameter_synthesis import morph_functions
from diameter_synthesis import utils
from diameter_synthesis.distribution_fitting import sample_distribution
from diameter_synthesis.exception import DiameterSynthesisError

TRUNK_FRAC_DECREASE = 0.1
N_TRIES_BEFORE_REDUC = 5
L = logging.getLogger(__name__)

STR_TO_TYPES = {
    "basal": SectionType.basal_dendrite,
    "apical": SectionType.apical_dendrite,
    "axon": SectionType.axon,
}

TYPES_TO_STR = {
    SectionType.basal_dendrite: "basal",
    SectionType.apical_dendrite: "apical",
    SectionType.axon: "axon",
}


def _reset_caches():
    """Reset the cached functions."""
    morph_functions.sec_length.cache_clear()
    morph_functions.partition_asymmetry_length.cache_clear()
    morph_functions.lengths_from_origin.cache_clear()
    morph_functions.n_children_downstream.cache_clear()


def _get_neurites(neuron, neurite_type):
    """Get a list of neurites to diametrize.

    Args:
        neuron (morphio.mut.Morphology): a neuron.
        neurite_type (morphio.SectionType): the neurite type to consider.

    Returns:
        list: list of neurites to consider.
    """
    return [
        list(neurite.iter()) for neurite in neuron.root_sections if neurite.type == neurite_type
    ]


def _sample_sibling_ratio(
    params, neurite_type, apply_asymmetry=False, mode="generic", rng=np.random
):
    """Sample a sibling ratio from distribution.

    Args:
        params (dict): model parameters.
        neurite_type (str): the neurite type to consider.
        apply_asymmetry (bool): asymmetry of current branching point.
        mode (str): to use or not the `asymmetry_threshold`.
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        float: sibling ratio.
    """
    if mode == "generic":
        return sample_distribution(params["sibling_ratios"][neurite_type], rng=rng)
    if mode == "threshold":
        if apply_asymmetry:
            return 0.0
        return sample_distribution(params["sibling_ratios"][neurite_type], rng=rng)
    # This case should never happen since the mode is already checked in `_select_model`
    raise DiameterSynthesisError(f"mode not understood {mode}")


def _sample_diameter_power_relation(
    params, neurite_type, apply_asymmetry=False, mode="generic", rng=np.random
):
    """Sample a diameter power relation from distribution.

    Args:
        params (dict): model parameters.
        neurite_type (str): the neurite type to consider.
        apply_asymmetry (bool): asymmetry of current branching point.
        mode (str): to use or not the `asymmetry_threshold`.
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        float: diameter power relation.
    """
    if mode == "generic":
        return sample_distribution(params["diameter_power_relation"][neurite_type], rng=rng)
    if mode == "threshold":
        if apply_asymmetry:
            return 1.0
        return sample_distribution(params["diameter_power_relation"][neurite_type], rng=rng)
    if mode == "exact":
        # This case should never happen since this mode is not known by `_select_model`
        return 1.0
    # This case should never happen since the mode is already checked in `_select_model`
    raise DiameterSynthesisError(f"mode not understood {mode}")


def _sample_trunk_diameter(params, neurite_type, rng=np.random):
    """Sample a trunk diameter from distribution.

    Args:
        params (dict): model parameters.
        neurite_type (str): the neurite type to consider.
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        float: trunk diameter
    """
    return sample_distribution(params["trunk_diameters"][neurite_type], rng=rng)


def _sample_terminal_diameter(params, neurite_type, rng=np.random):
    """Sample a terminal diameter.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        float: terminal diameter
    """
    return sample_distribution(params["terminal_diameters"][neurite_type], rng=rng)


def _sample_taper(params, neurite_type, rng=np.random):
    """Sample a taper rate from distributions.

    Args:
        params (dict): model parameters
        neurite_type (str): the neurite type to consider
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        float: taper rate
    """
    return sample_distribution(params["tapers"][neurite_type], rng=rng)


def _sample_daughter_diameters(section, params, params_tree, rng=np.random):
    """Compute the daughter diameters of the current section.

    Args:
        section (morphio.Section): section to consider.
        params (dict): model parameters.
        params_tree (dict): specific parameters of the current tree.
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
       list: list of daughter diameters.
    """
    # pylint: disable=too-many-locals
    major_sections = params_tree["major_sections"]

    apply_asymmetry = False
    if params_tree.get("with_asymmetry", False):
        child_not_in_major = [child.id not in major_sections for child in section.children]
        # if children are not asymmetrical, don't apply asymmetry
        if False in child_not_in_major:
            child_sort = np.argsort(child_not_in_major)
            apply_asymmetry = True

    reduction_factor = params_tree["reduction_factor_max"] + 1.0
    # try until we get a reduction of diameter in the branching
    while reduction_factor > params_tree["reduction_factor_max"]:

        sibling_ratio = _sample_sibling_ratio(
            params,
            params_tree["neurite_type"],
            apply_asymmetry=apply_asymmetry,
            mode=params_tree["mode_sibling"],
            rng=rng,
        )

        diameter_power_relation = _sample_diameter_power_relation(
            params,
            params_tree["neurite_type"],
            apply_asymmetry=apply_asymmetry,
            mode=params_tree["mode_diameter_power_relation"],
            rng=rng,
        )

        reduction_factor = morph_functions.diameter_power_relation_factor(
            diameter_power_relation, sibling_ratio
        )

    diam_0 = section.diameters[-1]
    terminal_diam = min(diam_0, params_tree["terminal_diam"])

    diam_1 = reduction_factor * diam_0
    diam_2 = sibling_ratio * diam_1

    diam_1 = max(diam_1, terminal_diam)
    diam_2 = max(diam_2, terminal_diam)

    diams = [diam_1] + (len(section.children) - 1) * [diam_2]

    if apply_asymmetry:
        return list(np.array(diams)[child_sort])

    # At the moment we don't have enough information to do better than a random choice in this case
    rng.shuffle(diams)
    return diams


def _diametrize_section(section, initial_diam, taper, min_diam=0.07, max_diam=10.0):
    """Diameterize a section.

    Args:
        section (morphio.Section): current section.
        initial_diam (float): initial diameter.
        taper (float): taper rate.
        min_diam (float): minimum diameter.
        max_diam (float): maximum diameter.
    """
    diams = polynomial.polyval(morph_functions.lengths_from_origin(section), [initial_diam, taper])
    section.diameters = np.clip(diams, min_diam, max_diam)


def _diametrize_tree(neurite, params, params_tree, rng=np.random):
    """Diametrize a tree, or neurite.

    Args:
        neurite (morphio neurite): current neurite.
        params (dict): model parameters.
        params_tree (dict): specific parameters of the current tree.
        rng (numpy.random.Generator): the random number generator to use.

    Returns:
        bool: `True` is all terminal diameters are small enough, `False` otherwise.
    """
    params_tree["tot_length"] = morph_functions.get_total_length(neurite)
    max_diam = params["terminal_diameters"][params_tree["neurite_type"]]["params"]["max"]
    wrong_tips = False
    active = deque([neurite[0]])
    while active:
        section = active.popleft()

        if section.is_root:
            init_diam = params_tree["trunk_diam"]
        else:
            init_diam = section.diameters[0]

        taper = _sample_taper(params, params_tree["neurite_type"], rng=rng)

        params_tree["terminal_diam"] = min(
            init_diam, _sample_terminal_diameter(params, params_tree["neurite_type"], rng=rng)
        )

        _diametrize_section(
            section,
            init_diam,
            taper=taper,
            min_diam=params_tree["terminal_diam"],
            max_diam=params_tree["trunk_diam"],
        )

        if len(section.children) > 0:
            diams = _sample_daughter_diameters(section, params, params_tree, rng=rng)

            for i, child in enumerate(section.children):
                utils.redefine_diameter_section(child, 0, diams[i])
                active.append(child)

        # if we are at a tip, check if tip diameters are small enough
        elif section.diameters[-1] > max_diam:
            wrong_tips = True

    return wrong_tips


def _diametrize_neuron(params_tree, neuron, params, neurite_types, config, rng=np.random):
    """Diametrize a neuron.

    Args:
        params_tree (dict): specific parameters of the current tree.
        neuron (morphio.mut.Morphology): neuron to diametrize.
        params (dict): model parameters.
        neurite_types (str or morphio.SectionType): the neurite type to consider.
        config (dict): general configuration parameters.
        rng (numpy.random.Generator): the random number generator to use.
    """
    # pylint: disable=too-many-locals, too-many-branches
    major_sections = set()
    if params_tree["with_asymmetry"]:
        # Get sections on the major branch
        for apical_section in params.get("apical_point_sec_ids", []):
            for sec in neuron.sections[apical_section].iter(IterType.upstream):
                major_sections.add(sec.id)
    params_tree["major_sections"] = major_sections

    for neurite_type in neurite_types:
        if isinstance(neurite_type, str):
            morphio_neurite_type = STR_TO_TYPES[neurite_type]
        else:
            morphio_neurite_type = neurite_type
            neurite_type = TYPES_TO_STR[neurite_type]

        params_tree["neurite_type"] = neurite_type

        for neurite in _get_neurites(neuron, morphio_neurite_type):
            wrong_tips = True
            n_tries = 0
            trunk_diam_frac = 1.0
            n_tries_step = 1

            _params = deepcopy(params) if "taper_increase" in params_tree else params
            while wrong_tips:
                if "taper_increase" in params_tree:
                    _params["tapers"][neurite_type]["params"]["scale"] *= params_tree[
                        "taper_increase"
                    ]
                    params_tree["trunk_diam"] = neurite[0].diameters[0]
                else:
                    trunk_diam = trunk_diam_frac * _sample_trunk_diameter(
                        _params, neurite_type, rng=rng
                    )
                    if trunk_diam < 0.01:
                        trunk_diam = 1.0
                        L.warning("sampled trunk diameter < 0.01, so use 1 instead")
                    params_tree["trunk_diam"] = trunk_diam

                wrong_tips = _diametrize_tree(neurite, _params, params_tree, rng=rng)

                # if we can't get a good model, reduce the trunk diameter progressively
                if n_tries > N_TRIES_BEFORE_REDUC * n_tries_step:
                    trunk_diam_frac -= TRUNK_FRAC_DECREASE
                    n_tries_step += 1
                if n_tries > config.get("trunk_max_tries", 100):
                    L.warning("max tries attained with %s", neurite_type)
                    wrong_tips = False
                n_tries += 1


def _select_model(model):
    """Select a diametrized model to use.

    Args:
        model (str): model name.

    Returns:
        function: diametrizer with specific `params_tree`.
    """
    if model == "generic":
        params_tree = {}
        params_tree["mode_sibling"] = "threshold"
        params_tree["mode_diameter_power_relation"] = "threshold"
        params_tree["with_asymmetry"] = True
        params_tree["reduction_factor_max"] = 1.0
    elif model == "astrocyte":
        params_tree = {}
        params_tree["mode_sibling"] = "generic"
        params_tree["mode_diameter_power_relation"] = "generic"
        params_tree["with_asymmetry"] = True
        params_tree["reduction_factor_max"] = 3.0
    elif model == "neurite_based":
        params_tree = {}
        params_tree["mode_sibling"] = "threshold"
        params_tree["mode_diameter_power_relation"] = "threshold"
        params_tree["with_asymmetry"] = True
        params_tree["reduction_factor_max"] = 1.0
        params_tree["taper_increase"] = 1.02

    else:
        raise DiameterSynthesisError(f"Unknown diameter model: {model}")

    return partial(_diametrize_neuron, params_tree)


def build(neuron, neurite_types, model_params, diam_params, random_generator=np.random):
    """Builder function for generating diameters of a neuron from the a diameter models.

    Args:
        neuron (morphio.mut.Morphology): neuron to diametrize.
        neurite_types (list[str] or str or list[SectionType] or SectionType): the neurite type(s)
            to consider.
        model_params (dict): model parameters.
        diam_params (dict): general configuration parameters.
        random_generator (numpy.random.Generator): the random number generator to use.
    """
    if isinstance(neurite_types, (str, SectionType)):
        neurite_types = [neurite_types]

    if "seed" in diam_params:
        random_generator = np.random.default_rng(diam_params["seed"])

    _reset_caches()

    if len(diam_params["models"]) > 1:
        L.warning("Several models provided, we will only use the first")
    diameter_generator = _select_model(diam_params["models"][0])

    diameter_generator(neuron, model_params, neurite_types, diam_params, rng=random_generator)
    n_samples = diam_params.get("n_samples", 1)
    if n_samples > 1:
        diameters = utils.get_all_diameters(neuron)
        for _ in range(n_samples - 1):
            diameter_generator(
                neuron, model_params, neurite_types, diam_params, rng=random_generator
            )
            for i, new_diams in enumerate(utils.get_all_diameters(neuron)):
                diameters[i] += new_diams
        for i, _ in enumerate(diameters):
            diameters[i] /= n_samples
        utils.set_all_diameters(neuron, diameters)


def _save_first_diams(morphology, length):
    """Save diameters in a dict up to length."""
    diams = {}
    for root_section in morphology.root_sections:
        if root_section.type == SectionType.axon:
            dist = 0
            prev_point = root_section.points[0]
            for section in root_section.iter():
                _diams = copy(section.diameters)
                for point in section.points:
                    dist += np.linalg.norm(point - prev_point)
                    prev_point = copy(point)
                    if dist >= length:
                        break
                diams[section.id] = _diams
                if dist >= length:
                    break

    return diams


def _set_first_diams(morphology, diams, length):
    """Set diameters from a dict up to length."""
    for root_section in morphology.root_sections:
        if root_section.type == SectionType.axon:
            dist = 0
            prev_point = root_section.points[0]
            for section in root_section.iter():
                current_diams = copy(section.diameters)
                old_diams = diams[section.id]
                for i, point in enumerate(section.points):
                    dist += np.linalg.norm(point - prev_point)
                    prev_point = copy(point)
                    current_diams[i] = old_diams[i]
                    if dist >= length:
                        break
                section.diameters = current_diams
                if dist >= length:
                    break


def diametrize_axon(
    morphology,
    main_diameter=1.0,
    collateral_diameter=0.1,
    main_taper=-0.0005,
    axon_point_isec=None,
    ais_length=60,
    rng=np.random,
):
    """Diametrize axon in place without learning from reconstructed axons.

    The main axon branch (from soma to axon point) will have a tapered diameter, and the
    colaterals a constant diameter.

    If an axon point is not provided, and main_diameter > collateral_diameter, the diameters will
    decrease with taper and bifurcations, with hardcoded parameters sibling_ratio = 0.5 and
    diameter_power_relation = 0.5.

    Args:
        morphology (morphio.mut.Morphology): morphology to diametrize.
        main_diameter (float): diameter of main axon branch (from soma to `axon_point_isec`).
        collateral_diameter (float): diameter of collateral branches.
        main_taper (float): taper rate of main branch (set to 0 for no taper, should be negative).
        axon_point_isec (int): morphio section id of axon point (see :mod:`morph_tool.axon_point`
            module).
        ais_length (float): length of ais for which we keep original diameters.
        rng (numpy.random.Generator): the random number generator to use.
    """
    model_params = {
        "trunk_diameters": {
            "axon": {
                "distribution": "constant",
                "params": {"value": main_diameter},
            }
        },
        "terminal_diameters": {
            "axon": {
                "distribution": "constant",
                "params": {"value": collateral_diameter, "max": main_diameter},
            }
        },
        "tapers": {
            "axon": {
                "distribution": "constant",
                "params": {"value": main_taper},
            }
        },
        "sibling_ratios": {
            "axon": {
                "distribution": "constant",
                "params": {"value": 0.5},
            }
        },
        "diameter_power_relation": {
            "axon": {
                "distribution": "constant",
                "params": {"value": 5.0},
            }
        },
    }
    if axon_point_isec:
        model_params["apical_point_sec_ids"] = [axon_point_isec]

    config = {"models": ["generic"], "trunk_max_tries": 1, "n_samples": 1}

    diams = _save_first_diams(morphology, ais_length)
    build(morphology, "axon", model_params, diam_params=config, random_generator=rng)
    _set_first_diams(morphology, diams, ais_length)
