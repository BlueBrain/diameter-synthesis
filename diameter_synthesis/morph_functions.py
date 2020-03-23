"""morphology extraction functions"""
import logging
from functools import lru_cache

import neurom as nm
import numpy as np
from neurom.core import iter_sections

from diameter_synthesis.exception import DiameterSynthesisError
from diameter_synthesis.utils import _get_diameters, _get_mean_diameter

L = logging.getLogger(__name__)


@lru_cache(maxsize=None)
def _sec_length(section):
    """length of a section (morphio)"""
    return np.linalg.norm(np.diff(np.asarray(section.points), axis=0), axis=1).sum()


@lru_cache(maxsize=None)
def _lengths_from_origin(section):
    """path lengths from first point of section (morphio)"""
    return np.insert(
        np.cumsum(np.linalg.norm(np.diff(np.asarray(section.points), axis=0), axis=1)),
        0,
        0,
    )


@lru_cache(maxsize=None)
def _partition_asymetry_length(section):
    """compute partition asymetry with lengths (morphio).
    WARNING: Requires to divide by total length of dendrite"""
    asymetry_length = sum(
        [_sec_length(section) for section in section.children[0].iter()]
    )
    asymetry_length -= sum(
        [_sec_length(section) for section in section.children[1].iter()]
    )
    return abs(asymetry_length)


@lru_cache(maxsize=None)
def _child_length(section):
    """number of children of a section (morphio)"""
    return sum(1 for _ in section.iter())


def _get_total_length(neurite):
    """get total length or a neurite (morphio)"""
    return np.sum(
        [
            np.linalg.norm(np.diff(np.asarray(section.points), axis=0), axis=1).sum()
            for section in neurite
        ]
    )


def compute_sibling_ratios(neurite, method="mean", attribute_name=None, bounds=None):
    """compute the sibling ratios of a neurite (neurom)"""
    sibling_ratios_values = nm.get("sibling_ratio", neurite, method=method)

    if not bounds:
        bounds = [0, 1 + 1e-5]
    return add_additional_attributes(
        sibling_ratios_values, neurite, attribute_name=attribute_name, bounds=bounds
    )


def compute_diameter_power_relation(
    neurite, method="mean", attribute_name=None, bounds=None
):
    """Returns the Rall deviation the diameter of the segments of a tree (neurom)"""
    diameter_power_relation_values = nm.get(
        "diameter_power_relation", neurite, method=method
    )
    if not bounds:
        bounds = [0.0, 100.0]

    return add_additional_attributes(
        diameter_power_relation_values,
        neurite,
        attribute_name=attribute_name,
        bounds=bounds,
    )


def diameter_power_relation_factor(diameter_power_relation, sibling_ratio):
    """Returns the reduction factor for bifurcation diameter"""
    if sibling_ratio > 0:
        return ((1.0 + sibling_ratio ** -1.5) / diameter_power_relation) ** (2.0 / 3.0)
    return 1.0


def terminal_diameters(
    neurite, method="mean", threshold=1.0, attribute_name=None, bounds=None
):
    """Returns the terminal diameters of a neurite (neurom)"""
    diameters = []
    for section in iter_sections(neurite):
        diameters += list(_get_diameters(section))
    mean_diameter = np.mean(diameters)

    if method == "mean":
        term_diam = [
            np.mean(_get_diameters(t))
            for t in nm.core.Tree.ileaf(next(iter_sections(neurite)))
            if np.mean(_get_diameters(t)) < threshold * mean_diameter
        ]

    elif method == "first":
        term_diam = [
            _get_diameters(t)[-1]
            for t in nm.core.Tree.ileaf(next(iter_sections(neurite)))
            if _get_diameters(t)[-1] < threshold * mean_diameter
        ]

    else:
        raise DiameterSynthesisError("Method for singling computation not understood!")

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        term_diam, neurite, attribute_name=attribute_name, bounds=bounds
    )


def min_diameter(neurite, attribute_name=None, bounds=None):
    """get the min diameter of a neurite (neurom)"""
    min_diam = 1e5
    for section in iter_sections(neurite):
        min_diam = min(min_diam, _get_mean_diameter(section))

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        [min_diam], neurite, attribute_name=attribute_name, bounds=bounds
    )


def max_diameter(neurite, attribute_name=None, bounds=None):
    """get the max diameter of a neurite (neurom) """
    max_diam = 0
    for section in iter_sections(neurite):
        max_diam = max(max_diam, _get_mean_diameter(section))

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        [max_diam], neurite, attribute_name=attribute_name, bounds=bounds
    )


def trunk_diameter(neurite, attribute_name=None, bounds=None, method="last"):
    """get the trunc diameters (neurom)"""
    if method == "mean":
        trunk_diam = _get_mean_diameter(neurite.root_node)
    if method == "first":
        trunk_diam = _get_diameters(neurite.root_node)[0]
    if method == "last":
        trunk_diam = _get_diameters(neurite.root_node)[-1]

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        [trunk_diam], neurite, attribute_name=attribute_name, bounds=bounds
    )


def taper(neurite, params=None, attribute_name=None):
    """get the taper rates (neurom)"""

    if params is None:
        raise DiameterSynthesisError(
            "Please provide parameters for taper rate computations"
        )

    tapers = nm.get("section_taper_rates", neurite)
    tapers = tapers[abs(tapers) > params["zeros"]]
    return add_additional_attributes(
        tapers,
        neurite,
        attribute_name=attribute_name,
        bounds=[params["min"], params["max"]],
    )


def get_additional_attribute(
    attribute_name, neurite=None, section=None
):  # noqa, pylint: disable=too-many-return-statements,too-many-branches
    """return the value of an additional attribute of a parameter, gien a neurite or a section"""
    if attribute_name in ["asymmetry", "asymmetry_threshold"]:
        if neurite is not None and section is None:
            out = nm.get("partition_asymmetry_length", [neurite])
            return np.array(out)
        if section is not None and neurite is None:
            return _partition_asymetry_length(section)

    if attribute_name == "asymmetry_pair":
        if section is not None and neurite is None:
            return [nm.features.bifurcationfunc.partition_pair(section)]
        raise DiameterSynthesisError("Please provide only a section")

    if attribute_name == "tot_length":
        if neurite is not None and section is None:
            return nm.get("total_length", neurite)
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "max_path":
        if neurite is not None and section is None:
            return [np.max(nm.get("section_path_distances", neurite))]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "max_branch":
        if neurite is not None and section is None:
            return [np.max(nm.get("section_branch_orders", neurite))]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "root_strahler":
        if neurite is not None and section is None:
            return [nm.get("section_strahler_orders", neurite)[0]]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "sibling":
        if neurite is not None and section is None:
            return compute_sibling_ratios(neurite, method="mean")
        if section is not None and neurite is None:
            return nm.features.bifurcationfunc.sibling_ratio(section, method="mean")

    raise DiameterSynthesisError(
        "Please provide either a neurite or a section, not both"
    )


def add_additional_attributes(data, neurite, attribute_name=None, bounds=None):
    """from a data and a type of additional, return both within some bounds"""
    if not bounds:
        bounds = [-1000, 1000]

    if attribute_name is None:
        data = np.array(data)
        return list(data[(data >= bounds[0]) & (data < bounds[1])])

    return [
        [d, s]
        for d, s in zip(data, get_additional_attribute(attribute_name, neurite=neurite))
        if bounds[0] <= d < bounds[1]
    ]
