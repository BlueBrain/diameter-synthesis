import os, glob, shutil
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

def sibling_ratios(neurite, method = 'mean'):
    """ compute the siblig ratios of a neurite"""

    s_ratios = []
    #loop over bifuraction points
    for bif_point in nm.core.Tree.ibifurcation_point(next(iter_sections(neurite))):
        if len(bif_point.children) == 2:

            if method == 'mean':
                #take the average diameter of children to smooth noise out
                #d1 = np.mean(get_diameters(bif_point.children[0]))
                #d2 = np.mean(get_diameters(bif_point.children[1]))

                d1 = get_mean_diameter(bif_point.children[0])
                d2 = get_mean_diameter(bif_point.children[1])

            elif method == 'first':
                #take the first diameter, but subject to noise!
                d1 = get_diameters(bif_point.children[0])[0]
                d2 = get_diameters(bif_point.children[1])[0]

            else:
                raise Exception('Method for singling computation not understood!')

            
            s_ratios.append(np.min([d1,d2]) / np.max([d1,d2]))

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return s_ratios

def Rall_deviations(neurite, method = 'mean'):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    
    Rall_deviations = []
    for bif_point in nm.core.Tree.ibifurcation_point(next(iter_sections(neurite))):
        if len(bif_point.children) == 2:

            if method == 'mean':
                #d_0 = np.mean(get_diameters(bif_point))
                #d_1 = np.mean(get_diameters(bif_point.children[0]))
                #d_2 = np.mean(get_diameters(bif_point.children[1]))

                d_0 = get_mean_diameter(bif_point)
                d_1 = get_mean_diameter(bif_point.children[0])
                d_2 = get_mean_diameter(bif_point.children[1])

            elif method == 'first':
                d_0 = get_diameters(bif_point)[-1]
                d_1 = get_diameters(bif_point.children[0])[0]
                d_2 = get_diameters(bif_point.children[1])[0]

            Rall_deviation = (d_1/d_0)**(3./2.) + (d_2/d_0)**(3./2.) 

            #don't consider cases with larger diameters of daughters 
            if Rall_deviation < 1.8:
                Rall_deviations.append(Rall_deviation)

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return Rall_deviations

def Rall_reduction_factor(Rall_deviation, siblings_ratio):
    '''Returns the reduction factor for bifurcation diameter'''

    return (Rall_deviation / (1. + siblings_ratio**(3./2.) ) )**(2./3.)

def terminal_diameters(neurite, method = 'mean', threshold = 1.0):
    """Returns the model for the terminations"""

    mean_diameter = np.mean(get_diameters(neurite))

    if method == 'mean':
        term_diam = [np.mean(get_diameters(t)) for t in nm.core.Tree.ileaf(next(iter_sections(neurite))) if np.mean(get_diameters(t)) < threshold * mean_diameter]

    elif method == 'first':
        term_diam = [get_diameters(t)[-1] for t in nm.core.Tree.ileaf(next(iter_sections(neurite))) if get_diameters(t)[-1] < threshold * mean_diameter]

    else:
        raise Exception('Method for singling computation not understood!')

    return term_diam

def trunk_diameter(neurite):
    """ get the trunc diameters """

    trunk_diam =  get_diameters(neurite.root_node)[0] 
    #max_bo = np.max(nm.get('section_term_branch_orders', neurite))
    max_path_dist = np.max(nm.get('section_path_distances', neurite))

    return [[trunk_diam, max_path_dist], ]


def taper(neurite, min_num_points = 20, fit_order = 1, params = None):
    """ get the taper """

    import pylab as plt

    tapers = []
    sec_id = [] 
    for i, section in enumerate(iter_sections(neurite)):
        lengths = [0] + section_lengths(section)

        #do a linear fit if more than 5 points
        if len(lengths) > min_num_points:
            z = polynomial.polyfit(lengths, get_diameters(section), fit_order, full=True)

            """
            plt.figure()
            plt.plot(lengths, get_diameters(section),'+')
            plt.plot(lengths, np.poly1d(z[0])(lengths))
            plt.title(str(np.round(z[1][0],3)))
            plt.savefig('taper/fig_'+str(i)+'.png')
            plt.close()
            """

            tap = z[0][-1]

            #print(z)
            #if tap < 0.010 and tap > -0.01 and abs(tap)>1e-8:
            if z[1][0]<params['max_residual'] and abs(tap)>params['zeros'] and params['min'] < -tap < params['max'] :
                tapers.append(-tap)
                sec_id.append(i)

    """ 
    bos = np.array(nm.get('section_branch_orders', neurite))[sec_id]
    lens =np.array(nm.get('section_lengths', neurite))[sec_id]
    plt.figure()
    tapers = np.array(tapers)
    print(lens,tapers)
    plt.plot(lens, tapers, '+')
    plt.axis([0,np.max(lens),-np.max(abs(tapers)), np.max(abs(tapers))])
    plt.savefig('taper/fig_'+str(j)+'.png')
    plt.close()
    """ 
    
    #bos = np.array(nm.get('section_path_distances', neurite))[sec_id]
    #return [ [tap, bo ] for tap, bo in zip(tapers, bos)  ]

    return  tapers
