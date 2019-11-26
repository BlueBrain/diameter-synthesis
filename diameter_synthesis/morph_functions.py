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

def sibling_ratios(neurite, method = 'mean', seq = None):
    """ compute the siblig ratios of a neurite"""

    s_ratios = []
    #loop over bifuraction points
    for bif_point in nm.core.Tree.ibifurcation_point(next(iter_sections(neurite))):
        if len(bif_point.children) == 2:
            if method == 'mean':
                #take the average diameter of children to smooth noise out
                d1 = get_mean_diameter(bif_point.children[0])
                d2 = get_mean_diameter(bif_point.children[1])

            elif method == 'first':
                #take the first diameter, but subject to noise!
                d1 = get_diameters(bif_point.children[0])[0]
                d2 = get_diameters(bif_point.children[1])[0]

            else:
                raise Exception('Method for sibling computation not understood!')

            #print(nm.features.bifurcationfunc.partition_pair(bif_point))
            s_ratios.append(np.min([d1,d2]) / np.max([d1,d2]))

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    #print('kl', nm.get('partition_asymmetry', neurite), s_ratios)
    #return s_ratios
    return sequential(s_ratios, seq, neurite, bounds = [0, 1-1e-5])

def sequential_single(seq, neurite = None, section = None):
    """ return the value for sequencial slicing"""

    if seq == 'asymmetry': 
        if neurite is not None and section is None:
            val_min = -2
            out = []
            for bif_point in iter_sections(neurite):

                if len(bif_point.children) > 1:

                    n = float(np.sum([sec.area for sec in bif_point.children[0].ipreorder()]))
                    m = float(np.sum([sec.area for sec in bif_point.children[1].ipreorder()]))

                    if len(bif_point.children) == 3:
                        m2 = float(np.sum([sec.area for sec in bif_point.children[2].ipreorder()]))

                        value = np.log10(np.min([n,m,m2]) / np.max([n,m,m2]))
                        out.append(max(val_min, value))
                    else:
                        value = np.log10(min(n,m) / max(n,m))
                        out.append(max(val_min, value))

                elif len(bif_point.children) == 1:
                    print('single child!')

            #print(out_max, np.mean(out_mean))
            return np.array(out)#/out_max#*out_max/np.mean(out_mean))

        elif section is not None and neurite is None:
            if len(section.children) > 1:

                n = float(np.sum([sec.area for sec in section.children[0].ipreorder()]))
                m = float(np.sum([sec.area for sec in section.children[1].ipreorder()]))

                if len(section.children) == 3:
                    m2 = float(np.sum([sec.area for sec in section.children[2].ipreorder()]))
                    
                    out = np.log10(np.min([n,m,m2]) / np.max([n,m,m2]))
                else:
                    out = np.log10(min(n,m) / max(n,m))
                       
                return out

        else:
            raise Exception('Please provide either a neurite or a section, not both')

    elif seq == 'asymmetry_pair': 

        if section is not None and neurite is None:
            try:
                return [nm.features.bifurcationfunc.partition_pair(section),]
            except:
                #cathing the fact that some bifurcation are triple
                print('triple bifurcation for partition asymetry')
                return [(0.5,1)]
        else:
            raise Exception('Please provide either a section only')

    elif seq == 'max_path':
        if neurite is not None and section is None:
            return [np.max(nm.get('section_path_distances', neurite)),]
        else:
            raise Exception('Please provide only a neuritet')

    elif seq == 'sibling': 
        if neurite is not None and section is None:
            return sibling_ratios(neurite, method = 'mean', seq = None)

        elif section is not None and neurite is None:
            d1 = get_mean_diameter(section.children[0])
            d2 = get_mean_diameter(section.children[1])

            return min(d1, d2) / max(d1, d2)


    else:
        return [0,] #not the best way to do that...

def sequential(data, seq, neurite, bounds = [-1000, 1000]):
    """ from a data and a type of sequential slicing, return both """
    if not isinstance(seq, str):
        return data 

    else:
        return [[d, s] for d, s in zip(data, sequential_single(seq, neurite)) if bounds[0] < d < bounds[1] ]

def Rall_deviations(neurite, method = 'mean', seq = None):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    
    Rall_deviations = []
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

            Rall_deviation = (d_1/d_0)**(3./2.) + (d_2/d_0)**(3./2.) 

            #don't consider cases with larger diameters of daughters 
            #if 0.8 < Rall_deviation < 1.2:
            Rall_deviations.append(Rall_deviation)

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return sequential(Rall_deviations, seq, neurite, bounds = [0, 12-1e-5])

def Rall_reduction_factor(Rall_deviation, siblings_ratio, seq = None):
    '''Returns the reduction factor for bifurcation diameter'''

    return (Rall_deviation / (1. + siblings_ratio**(3./2.) ) )**(2./3.)

def terminal_diameters(neurite, method = 'mean', threshold = 1.0, seq = None):
    """Returns the model for the terminations"""

    mean_diameter = np.mean(get_diameters(neurite))

    if method == 'mean':
        term_diam = [np.mean(get_diameters(t)) for t in nm.core.Tree.ileaf(next(iter_sections(neurite))) if np.mean(get_diameters(t)) < threshold * mean_diameter]

    elif method == 'first':
        term_diam = [get_diameters(t)[-1] for t in nm.core.Tree.ileaf(next(iter_sections(neurite))) if get_diameters(t)[-1] < threshold * mean_diameter]

    else:
        raise Exception('Method for singling computation not understood!')

    return sequential(term_diam, seq, neurite)

def trunk_diameter(neurite, seq = None):
    """ get the trunc diameters """

    trunk_diam =  get_diameters(neurite.root_node)[0] 

    return sequential([trunk_diam,], seq, neurite)

def taper(neurite, min_num_points = 20, fit_order = 1, params = None, seq = None):
    """ get the taper """

    #import pylab as plt

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

    return sequential(tapers, seq, neurite)
