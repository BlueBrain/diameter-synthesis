""" morphologe extraction functions """
import logging

import numpy as np
from numpy.polynomial import polynomial

import neurom as nm
from neurom.core import iter_sections

from diameter_synthesis.utils import (get_diameters, get_mean_diameter,
                                      section_lengths)

L = logging.getLogger(__name__)
##########################
# morphometric functions #
##########################


def sibling_ratios(neurite, method='mean', seq=None, bounds=None):
    """ compute the siblig ratios of a neurite"""
    sibling_ratios_values = nm.get('sibling_ratio', neurite, method=method)

    if not bounds:
        bounds = [0, 1 + 1e-5]

    return sequential(sibling_ratios_values, seq, neurite, bounds=bounds)


def sequential_single(seq, neurite=None, section=None):  # noqa, pylint: disable=too-many-return-statements,too-many-branches
    """ return the value for sequencial slicing"""
    if seq in ('asymmetry', 'asymmetry_threshold'):
        if neurite is not None and section is None:
            out = nm.get('partition_asymmetry_length', [neurite, ])
            return np.array(out)

        if section is not None and neurite is None:
            if len(section.children) > 1:

                var_n = float(sum([sec.length for sec in section.children[0].ipreorder()]))
                var_m = float(sum([sec.length for sec in section.children[1].ipreorder()]))

                if len(section.children) == 3:
                    var_m2 = float(np.sum([sec.length for sec in section.children[2].ipreorder()]))
                    out = abs(np.max([var_n, var_m, var_m2]) - np.min([var_n, var_m, var_m2]))
                else:
                    out = abs(var_n - var_m)
                return out
            return 1.

        raise Exception('Please provide either a neurite or a section, not both')

    if seq == 'asymmetry_pair':

        if section is not None and neurite is None:
            try:
                return [nm.features.bifurcationfunc.partition_pair(section), ]
            except BaseException:  # pylint: disable=broad-except
                # cathing the fact that some bifurcation are triple
                L.exception('triple bifurcation for partition asymetry')
                return [(0.5, 1)]
        else:
            raise Exception('Please provide either a section only')

    if seq == 'tot_length':
        if neurite is not None and section is None:
            return nm.get('total_length', neurite)
        raise Exception('Please provide only a neuritet')

    if seq == 'max_path':
        if neurite is not None and section is None:
            return [np.max(nm.get('section_path_distances', neurite)), ]
        raise Exception('Please provide only a neuritet')

    if seq == 'max_branch':
        if neurite is not None and section is None:
            return [np.max(nm.get('section_branch_orders', neurite)), ]
        raise Exception('Please provide only a neuritet')

    if seq == 'root_strahler':
        if neurite is not None and section is None:
            return [nm.get('section_strahler_orders', neurite)[0], ]
        raise Exception('Please provide only a neuritet')

    if seq == 'sibling':
        if neurite is not None and section is None:
            return sibling_ratios(neurite, method='mean', seq=None)

        if section is not None and neurite is None:
            return nm.features.bifurcationfunc.sibling_ratio(section, method='mean')

    return [0, ]  # not the best way to do that...


def sequential(data, seq, neurite, bounds=None):
    """ from a data and a type of sequential slicing, return both """

    if not bounds:
        bounds = [-1000, 1000]

    if not isinstance(seq, str):
        data = np.array(data)
        return list(data[(data > bounds[0]) & (data < bounds[1])])

    return [[d, s]
            for d, s in zip(data, sequential_single(seq, neurite)) if bounds[0] < d < bounds[1]]


def rall_deviations(neurite, method='mean', seq=None, bounds=None):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    rall_deviations_values = nm.get('diameter_power_relation', neurite, method=method)
    if not bounds:
        bounds = [0.0, 100.]

    return sequential(rall_deviations_values, seq, neurite, bounds=bounds)


def rall_reduction_factor(rall_deviation, sibling_ratio, seq=None):  # noqa, pylint: disable=unused-argument
    '''Returns the reduction factor for bifurcation diameter'''
    if sibling_ratio > 0:
        return (1. / rall_deviation * (1. + sibling_ratio**(-3. / 2.)))**(2. / 3.)
    return 1.


def terminal_diameters(neurite, method='mean', threshold=1.0, seq=None, bounds=None):
    """Returns the model for the terminations"""

    mean_diameter = np.mean(get_diameters(neurite))

    if method == 'mean':
        term_diam = [np.mean(get_diameters(t)) for t in nm.core.Tree.ileaf(
            next(iter_sections(neurite))) if np.mean(get_diameters(t)) < threshold * mean_diameter]

    elif method == 'first':
        term_diam = [get_diameters(t)[-1] for t in nm.core.Tree.ileaf(next(iter_sections(neurite)))
                     if get_diameters(t)[-1] < threshold * mean_diameter]

    else:
        raise Exception('Method for singling computation not understood!')

    if not bounds:
        bounds = [0, 100]

    return sequential(term_diam, seq, neurite, bounds=bounds)


def min_diameter(neurite, seq=None, bounds=None):
    """ get the min diameter of a neurite """
    min_diam = 1e5
    for section in iter_sections(neurite):
        min_diam = min(min_diam, get_mean_diameter(section))

    if not bounds:
        bounds = [0, 100]

    return sequential([min_diam, ], seq, neurite, bounds=bounds)


def max_diameter(neurite, seq=None, bounds=None):
    """ get the max diameter of a neurite """
    max_diam = 0
    for section in iter_sections(neurite):
        max_diam = max(max_diam, get_mean_diameter(section))

    if not bounds:
        bounds = [0, 100]

    return sequential([max_diam, ], seq, neurite, bounds=bounds)


def trunk_diameter(neurite, seq=None, method='last', bounds=None):
    """ get the trunc diameters """

    if method == 'mean':
        trunk_diam = get_mean_diameter(neurite.root_node)
    if method == 'first':
        trunk_diam = get_diameters(neurite.root_node)[0]
    if method == 'last':
        trunk_diam = get_diameters(neurite.root_node)[-1]

    if not bounds:
        bounds = [0, 100]

    return sequential([trunk_diam, ], seq, neurite, bounds=bounds)


def taper(neurite, params=None, seq=None, only_first=False):
    """ get the taper """

    fit_order = 1

    # to output the taper of the first section of the neurite only
    if only_first:
        root_section = list(iter_sections(neurite))[0]
        lengths = [0] + section_lengths(root_section)
        fit = polynomial.polyfit(lengths, get_diameters(root_section), fit_order, full=True)
        tap = fit[0][-1]
        if -.1 < tap < 0:  # only consider negative taper here
            return [-tap, ]
        return []

    tapers = nm.get('section_taper_rates', neurite)
    tapers = tapers[abs(tapers) > params['zeros']]
    return sequential(tapers, seq, neurite, bounds=[params['min'], params['max']])
