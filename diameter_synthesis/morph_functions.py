import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm
from neurom.core import iter_sections

from diameter_synthesis.utils import get_diameters, set_diameters 

############################
## morphometric functions ##
############################

def sibling_ratios(neurite, method = 'mean'):
    """ compute the siblig ratios of a neurite"""

    s_ratios = []
    #loop over bifuraction points
    for bif_point in nm.core.Tree.ibifurcation_point(iter_sections(neurite).next()):
        if len(bif_point.children) == 2:

            if method == 'mean':
                #take the average diameter of children to smooth noise out
                d1 = np.mean(get_diameters(bif_point.children[0]))
                d2 = np.mean(get_diameters(bif_point.children[1]))

            elif method == 'first':
                #take the first diameter, but subject to noise!
                d1 = get_diameters(bif_point.children[0])[0]
                d2 = get_diameters(bif_point.children[1])[0]

            else:
                raise Exception('Method for singling computation not understood!')

            s_ratios.append( np.min([d1,d2]) / np.max([d1,d2]) ) 

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return s_ratios

def Rall_deviations(neurite, method = 'mean'):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    
    Rall_deviations = []
    for bif_point in nm.core.Tree.ibifurcation_point(iter_sections(neurite).next()):
        if len(bif_point.children) == 2:

            if method == 'mean':

                d_0 = np.mean(get_diameters(bif_point))
                d_1 = np.mean(get_diameters(bif_point.children[0]))
                d_2 = np.mean(get_diameters(bif_point.children[1]))

            elif method == 'first':

                d_0 = get_diameters(bif_point)[-1]
                d_1 = get_diameters(bif_point.children[0])[-1]
                d_2 = get_diameters(bif_point.children[1])[-1]

            Rall_deviations.append( (d_1/d_0)**(3./2.) + (d_2/d_0)**(3./2.) )

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return Rall_deviations

def Rall_reduction_factor(Rall_deviation, siblings_ratio):
    '''Returns the reduction factor for bifurcation diameter'''

    return (Rall_deviation / (1. + siblings_ratio**(3./2.) ) )**(2./3.)

def terminal_diameters(neurite, method = 'mean', threshold = 0.8):
    """Returns the model for the terminations"""

    mean_diameter = np.mean(get_diameters(neurite))

    if method == 'mean':
        term_diam = [2. * np.mean(get_diameters(t)) for t in nm.core.Tree.ileaf(iter_sections(neurite).next()) if np.mean(get_diameters(t)) < threshold * mean_diameter]

    elif method == 'first':
        term_diam = [2. * get_diameters(t)[-1] for t in nm.core.Tree.ileaf(iter_sections(neurite).next()) if get_diameters(t)[-1] < threshold * mean_radii]

    else:
        raise Exception('Method for singling computation not understood!')

    return term_diam

def trunk_diameter(neurite):
    """ get the trunc diameters """

    trunk_diam =  get_diameters(neurite.root_node)[0] 
    max_bo = np.max(nm.get('section_term_branch_orders', neurite))

    return [[trunk_diam, max_bo], ]


