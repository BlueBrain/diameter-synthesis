"""Morphology extraction functions."""

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
from functools import lru_cache

import neurom as nm
import numpy as np
from neurom.core.morphology import Section
from neurom.core.morphology import iter_sections

from diameter_synthesis.exception import DiameterSynthesisError
from diameter_synthesis.utils import _get_diameters
from diameter_synthesis.utils import _get_mean_diameter

L = logging.getLogger(__name__)


def _segment_lengths(section):
    """Lengths of segments of a section (morphio only)."""
    points = section.points
    vectors = np.diff(points, axis=0)
    return np.linalg.norm(vectors, axis=1)


@lru_cache(maxsize=None)
def sec_length(section):
    """Length of a section (morphio only)."""
    return _segment_lengths(section).sum()


@lru_cache(maxsize=None)
def lengths_from_origin(section):
    """Path lengths from first point of section (morphio only)."""
    return np.insert(np.cumsum(_segment_lengths(section)), 0, 0)


@lru_cache(maxsize=None)
def partition_asymmetry_length(section):
    """Compute partition asymmetry with lengths (morphio).

    WARNING: This function is not complete, and requires to divide by total length of dendrite,
    which is done later in the code as we don't have access to whole dendrite here.
    """
    asymmetry_length = sum([sec_length(section) for section in section.children[0].iter()])
    try:
        asymmetry_length -= sum([sec_length(section) for section in section.children[1].iter()])
    except IndexError as exc:
        raise DiameterSynthesisError(
            "Bifurcation point with single child, consider using neuror sanitize"
        ) from exc
    return abs(asymmetry_length)


@lru_cache(maxsize=None)
def n_children_downstream(section):
    """Get the number of children of a section (morphio only)."""
    return sum(1 for _ in section.iter())


def get_total_length(neurite):
    """Get total length or a neurite (morphio only)."""
    return sum(map(sec_length, neurite))


def compute_sibling_ratios(neurite, method="mean", attribute_name=None, bounds=None):
    """Compute the sibling ratios of a neurite (neurom only)."""
    sibling_ratios_values = nm.get("sibling_ratios", neurite, method=method)

    if not bounds:
        bounds = [0, 1 + 1e-5]
    return add_additional_attributes(
        sibling_ratios_values, neurite, attribute_name=attribute_name, bounds=bounds
    )


def compute_diameter_power_relation(neurite, method="mean", attribute_name=None, bounds=None):
    """Compute the Rall deviation the diameter of the segments of a tree (neurom only)."""
    diameter_power_relation_values = nm.get("diameter_power_relations", neurite, method=method)
    if not bounds:
        bounds = [0.0, 100.0]

    return add_additional_attributes(
        diameter_power_relation_values,
        neurite,
        attribute_name=attribute_name,
        bounds=bounds,
    )


def diameter_power_relation_factor(diameter_power_relation, sibling_ratio):
    """Compute the reduction factor for bifurcation diameter."""
    if sibling_ratio > 0:
        return ((1.0 + sibling_ratio ** -1.5) / diameter_power_relation) ** (2.0 / 3.0)
    return 1.0


def terminal_diameters(neurite, method="mean", threshold=1.0, attribute_name=None, bounds=None):
    """Compute the terminal diameters of a neurite (neurom only)."""
    diameters = []
    for section in iter_sections(neurite):
        diameters += list(_get_diameters(section))
    mean_diameter = np.mean(diameters)

    if method == "mean":
        term_diam = [
            np.mean(_get_diameters(t))
            for t in Section.ileaf(next(iter_sections(neurite)))
            if np.mean(_get_diameters(t)) < threshold * mean_diameter
        ]

    elif method == "first":
        term_diam = [
            _get_diameters(t)[-1]
            for t in Section.ileaf(next(iter_sections(neurite)))
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
    """Get the min diameter of a neurite (neurom only)."""
    min_diam = min(map(_get_mean_diameter, iter_sections(neurite)))

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        [min_diam], neurite, attribute_name=attribute_name, bounds=bounds
    )


def max_diameter(neurite, attribute_name=None, bounds=None):
    """Get the max diameter of a neurite (neurom only)."""
    max_diam = max(map(_get_mean_diameter, iter_sections(neurite)))

    if not bounds:
        bounds = [0, 100]
    return add_additional_attributes(
        [max_diam], neurite, attribute_name=attribute_name, bounds=bounds
    )


def trunk_diameter(neurite, attribute_name=None, bounds=None, method="last"):
    """Get the trunc diameters (neurom only)."""
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


def taper(neurite, params, attribute_name=None):
    """Get the taper rates (neurom only)."""
    if params is None:
        raise DiameterSynthesisError("Please provide parameters for taper rate computations")
    tapers = nm.get("section_taper_rates", neurite)
    return add_additional_attributes(
        np.array(tapers),
        neurite,
        attribute_name=attribute_name,
        bounds=[params["min"], params["max"]],
    )


def get_additional_attribute(
    attribute_name, neurite=None, section=None
):  # noqa, pylint: disable=too-many-return-statements,too-many-branches
    """Return the value of an additional attribute of a parameter, given a neurite or a section."""
    section_only = section is not None and neurite is None
    neurite_only = neurite is not None and section is None
    if attribute_name in ["asymmetry", "asymmetry_threshold"]:
        if neurite_only:
            out = nm.get("partition_asymmetry", neurite, variant="length")
            return np.array(out)
        if section_only:
            return partition_asymmetry_length(section)

    if attribute_name == "asymmetry_pair":
        if section_only:
            return [nm.features.bifurcation.partition_pair(section)]
        raise DiameterSynthesisError("Please provide only a section")

    if attribute_name == "tot_length":
        if neurite_only:
            return [nm.get("total_length", neurite)]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "max_path":
        if neurite_only:
            return [np.max(nm.get("section_path_distances", neurite))]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "max_branch":
        if neurite_only:
            return [np.max(nm.get("section_branch_orders", neurite))]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "root_strahler":
        if neurite_only:
            return [nm.get("section_strahler_orders", neurite)[0]]
        raise DiameterSynthesisError("Please provide only a neurite")

    if attribute_name == "sibling":
        if neurite_only:
            return compute_sibling_ratios(neurite, method="mean")
        if section_only:
            return nm.features.bifurcation.sibling_ratio(section, method="mean")

    raise DiameterSynthesisError(
        "Please provide a valid attribute_name and provide either a neurite or a section, not both"
    )


def add_additional_attributes(data, neurite, attribute_name=None, bounds=None):
    """From a data and a type of additional, return both within some bounds."""
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
