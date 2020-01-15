import os
import glob
import shutil
import json
from tqdm import tqdm

import numpy as np
from numpy.polynomial import polynomial as polynomial

import neurom as nm
from neurom.core import iter_sections

from diameter_synthesis.utils import get_diameters, set_diameters, section_lengths, get_mean_diameter

############################
## morphometric functions ##
############################

def sibling_ratios(neurite, method='mean', seq=None, bounds=[0, 1 + 1e-5]):
    """ compute the siblig ratios of a neurite"""
    s_ratios = []
    # loop over bifuraction points
    for section in iter_sections(neurite):
        if len(section.children) > 1:
            ds = []
            for sec in section.children:
                # take the average diameter of children to smooth noise out
                if method == 'mean':
                    ds.append(get_mean_diameter(sec))

                elif method == 'first':
                    # take the first diameter, but subject to noise!
                    ds.append(get_diameters(sec)[0])

                else:
                    raise Exception('Method for sibling computation not understood!')

            s_ratios.append(np.min(ds) / np.max(ds))

    return sequential(s_ratios, seq, neurite, bounds=bounds)


def sequential_single(seq, neurite=None, section=None):
    """ return the value for sequencial slicing"""
    if seq == 'asymmetry' or seq == 'asymmetry_threshold':
        if neurite is not None and section is None:
            val_min = -2
            out = []
            tot_len = nm.get('total_length', neurite)[0]
            for bif_point in iter_sections(neurite):

                if len(bif_point.children) > 1:

                    n = float(np.sum([sec.length for sec in bif_point.children[0].ipreorder()]))
                    m = float(np.sum([sec.length for sec in bif_point.children[1].ipreorder()]))

                    if len(bif_point.children) == 3:
                        m2 = float(np.sum([sec.length for sec in bif_point.children[2].ipreorder()]))
                        out.append(abs(np.max([n, m, m2]) - np.min([n, m, m2])) / tot_len)
                    else:
                        out.append(abs(n - m) / tot_len)

                elif len(bif_point.children) == 1:
                    print('single child!')

            return np.array(out)

        elif section is not None and neurite is None:
            if len(section.children) > 1:

                n = float(np.sum([sec.length for sec in section.children[0].ipreorder()]))
                m = float(np.sum([sec.length for sec in section.children[1].ipreorder()]))

                if len(section.children) == 3:
                    m2 = float(np.sum([sec.length for sec in section.children[2].ipreorder()]))
                    out = abs(np.max([n, m, m2]) - np.min([n, m, m2]))
                else:
                    out = abs(n - m)
                return out
            else:
                #print('single children, use asymetry = 1')
                return 1.

        else:
            raise Exception('Please provide either a neurite or a section, not both')

    elif seq == 'asymmetry_pair':

        if section is not None and neurite is None:
            try:
                return [nm.features.bifurcationfunc.partition_pair(section), ]
            except:
                # cathing the fact that some bifurcation are triple
                print('triple bifurcation for partition asymetry')
                return [(0.5, 1)]
        else:
            raise Exception('Please provide either a section only')

    elif seq == 'tot_length':
        if neurite is not None and section is None:
            return nm.get('total_length', neurite)
        else:
            raise Exception('Please provide only a neuritet')


    elif seq == 'max_path':
        if neurite is not None and section is None:
            return [np.max(nm.get('section_path_distances', neurite)), ]
        else:
            raise Exception('Please provide only a neuritet')

    elif seq == 'max_branch':
        if neurite is not None and section is None:
            return [np.max(nm.get('section_branch_orders', neurite)), ]
        else:
            raise Exception('Please provide only a neuritet')

    elif seq == 'root_strahler':
        if neurite is not None and section is None:
            return [list(nm.features.neuritefunc.section_strahler_orders(neurite))[0],]
        else:
            raise Exception('Please provide only a neuritet')

    elif seq == 'sibling':
        if neurite is not None and section is None:
            return sibling_ratios(neurite, method='mean', seq=None)

        elif section is not None and neurite is None:
            d1 = get_mean_diameter(section.children[0])
            d2 = get_mean_diameter(section.children[1])

            return min(d1, d2) / max(d1, d2)

    else:
        return [0, ]  # not the best way to do that...


def sequential(data, seq, neurite, bounds=[-1000, 1000]):
    """ from a data and a type of sequential slicing, return both """
    if not isinstance(seq, str):
        data = np.array(data)
        return list(data[(data > bounds[0]) & (data < bounds[1])])

    else:
        return [[d, s] for d, s in zip(data, sequential_single(seq, neurite)) if bounds[0] < d < bounds[1]]


def rall_deviations(neurite, method='mean', seq=None, bounds=[0, 1.5 - 1e-5]):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    rall_deviations = []
    for bif_point in nm.core.Tree.ibifurcation_point(next(iter_sections(neurite))):
        if len(bif_point.children) == 2:

            if method == 'mean':
                d_0 = get_mean_diameter(bif_point)
                d_1 = get_mean_diameter(bif_point.children[0])
                d_2 = get_mean_diameter(bif_point.children[1])

            elif method == 'first':
                d_0 = get_diameters(bif_point)[-1]
                d_1 = get_diameters(bif_point.children[0])[0]
                d_2 = get_diameters(bif_point.children[1])[0]

            rall_deviation = (d_1 / d_0)**(1.5) + (d_2 / d_0)**(1.5)

            rall_deviations.append(rall_deviation)

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is ' + str(len(bif_point.children)) + '!')

    return sequential(rall_deviations, seq, neurite, bounds=bounds)


def rall_reduction_factor(rall_deviation, siblings_ratio, seq=None):
    '''Returns the reduction factor for bifurcation diameter'''
    return (rall_deviation / (1. + siblings_ratio**(3. / 2.)))**(2. / 3.)


def terminal_diameters(neurite, method='mean', threshold=1.0, seq=None, bounds=[0,100]):
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

    return sequential(term_diam, seq, neurite, bounds=bounds)


def min_diameter(neurite, seq=None, bounds=[0,100]):
    """ get the min diameter of a neurite """
    min_diam = 1e5
    for section in iter_sections(neurite):
        min_diam = min(min_diam, get_mean_diameter(section))

    return sequential([min_diam, ], seq, neurite, bounds=bounds)


def max_diameter(neurite, seq=None, bounds=[0,100]):
    """ get the max diameter of a neurite """
    max_diam = 0
    for section in iter_sections(neurite):
        max_diam = max(max_diam, get_mean_diameter(section))

    return sequential([max_diam, ], seq, neurite, bounds=bounds)


def trunk_diameter(neurite, seq=None, method='first', bounds=[0,100]):
    """ get the trunc diameters """

    if method == 'mean':
        trunk_diam = get_mean_diameter(neurite.root_node)
    if method == 'first':
        trunk_diam = get_diameters(neurite.root_node)[0]

    return sequential([trunk_diam, ], seq, neurite, bounds=bounds)


def taper(neurite, min_num_points=20, fit_order=1, params=None, seq=None):
    """ get the taper """
    tapers = []
    sec_id = []
    for i, section in enumerate(iter_sections(neurite)):
        lengths = [0] + section_lengths(section)

        # do a linear fit if more than 5 points
        if len(lengths) > min_num_points:
            z = polynomial.polyfit(lengths, get_diameters(section), fit_order, full=True)

            tap = z[0][-1]
            residual = z[1][0]
            if residual < params['max_residual'] and abs(tap) > params['zeros'] and params['min'] < tap < params['max']:
                tapers.append(tap)
                sec_id.append(i)

    return sequential(tapers, seq, neurite)
