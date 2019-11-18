import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm
from neurom import COLS
from neurom.core import iter_sections

from neurom import NeuriteType
from morphio import SectionType

import logging 
L = logging.getLogger(__name__)

STR_TO_TYPES = {'apical': SectionType.apical_dendrite,
                'basal': SectionType.basal_dendrite,
                'axon': SectionType.axon}

TYPE_TO_STR = {SectionType.apical_dendrite: 'apical',
               SectionType.basal_dendrite: 'basal',
               SectionType.axon: 'axon',
               SectionType.soma: 'soma'}


NEUROM_TYPE_TO_STR = {NeuriteType.apical_dendrite: 'apical',
                      NeuriteType.basal_dendrite: 'basal',
                      NeuriteType.soma: 'soma',
                      NeuriteType.axon: 'axon'}

STR_TO_NEUROM_TYPES = {'apical': NeuriteType.apical_dendrite,
                       'basal': NeuriteType.basal_dendrite,
                       'soma': NeuriteType.soma,
                       'axon': NeuriteType.axon}

ROUND = 4 #number of digits for the fitted parameters
MIN_DATA_POINTS = 20 #minimum number of points to fit a distribution
A_MAX = 4 #maximum value for the a (shape) parameter of fits (can get really large when low number of points)

FORBIDDEN_MTYPES = ['L4_NGC', 'L4_CHC', 'L6_CHC', 'L6_DBC', ]

##############################
## loading/saving functions ##
##############################

def load_morphologies_from_folder(morph_path):
    """ Load the morphologies from a list of files in a folder """


    morphologies = {}
    morphologies['generic_type'] = []
    for fname in tqdm(os.listdir(morph_path)):
        filepath = os.path.join(morph_path, fname)
        if fname.endswith(('.h5', '.asc', '.swc')) and os.path.exists(filepath):
            neuron = nm.load_neuron(filepath)
            morphologies['generic_type'].append([neuron, os.path.splitext(fname)[0]])

    return morphologies 


def load_morphologies_from_dict(morph_path, name_dict):
    """ Load the morphologies from a list of files """

    tqdm_1, tqdm_2 = tqdm_disable(name_dict) #to have a single progression bar
    # just to get single progression bar 
    morphologies = {}
    for mtype in tqdm(name_dict, disable = tqdm_1):
        morphologies[mtype] = []
        for fname in tqdm(name_dict[mtype], disable = tqdm_2):
            filepath = os.path.join(morph_path, fname)
            if fname.endswith(('.h5', '.asc', '.swc')) and os.path.exists(filepath):
                neuron = nm.load_neuron(filepath)
                morphologies[mtype].append([neuron, os.path.splitext(fname)[0]])

    return morphologies 

def create_morphologies_dict(morph_path, mtypes_sort = 'all', n_morphs_max = None, n_mtypes_max = None, xml_file = 'neuronDB.xml', ext = '.asc', prefix = "", super_mtypes_path='../scripts/diameter_types/'):
    """ Create dict to load the morphologies from a directory, by mtypes or all at once """

    #first load the neuronDB.xml file 
    from xml.etree import ElementTree as ET
    FileDB = ET.parse(morph_path + xml_file)
    root = FileDB.findall('listing')[0]
    morphs = root.findall('morphology')

    #if not max number of mtypes, take all
    if not n_mtypes_max:
        n_mtypes_max =  1e10
    if not n_morphs_max:
        n_morphs_max =  1e10

    name_dict = {}
    if mtypes_sort == 'all':
        name_dict['all_types'] = []
    
    if len(morphs)>0:
        for m in morphs:
            try:
                # Define mtypes
                mtype = m.find('mtype').text
                # Define subtypes (if they exist)
                if m.find('msubtype').text:
                    mtype = mtype + ':' + m.find('msubtype').text

                if mtype not in FORBIDDEN_MTYPES: #hack to not consider some bad mtypes

                    #load morphologies by mtypes, one list of morphologies for each type
                    if mtypes_sort == 'mtypes':
                        #if it is a new mtype, add en entry to name_dict
                        if mtype not in name_dict and len(name_dict) < n_mtypes_max:
                            name_dict[mtype] = [prefix + m.find('name').text + ext]
                        elif mtype in name_dict:
                            name_dict[mtype] += [prefix + m.find('name').text + ext]

                    #collect all of the mtypes in a single dict entry
                    elif mtypes_sort == 'all' and len(name_dict['all_types']) < n_morphs_max:
                        name_dict['all_types'] += [prefix + m.find('name').text + ext]
                    
                    #use super_mtypes if precomputed
                    elif mtypes_sort == 'super_mtypes':
                        super_mtypes_file = super_mtypes_path + 'super_mtypes.json'
                        with open(super_mtypes_file, 'r') as f:
                            super_mtypes_dict = json.load(f)
                        super_mtype = super_mtypes_dict[mtype]

                        if super_mtype not in name_dict:
                            name_dict[super_mtype] = [prefix + m.find('name').text + ext]
                        elif super_mtype in name_dict:
                            name_dict[super_mtype] += [prefix + m.find('name').text + ext]
            except Exception as e:
                print('Failed to process', e)

    return name_dict

def load_morphologies(morph_path, mtypes_sort = 'all', n_morphs_max = None, n_mtypes_max = None, xml_file = './neuronDB.xml', ext = '.asc', prefix = ""):
    """ Load the morphologies from a directory, by mtypes or all at once """

    # try to load the file, if there is no file, set a constant mtype
    try: 
        name_dict = create_morphologies_dict(morph_path, mtypes_sort = mtypes_sort, n_morphs_max = n_morphs_max, n_mtypes_max = n_mtypes_max, xml_file = xml_file, ext = ext, prefix = prefix)
        return load_morphologies_from_dict(morph_path, name_dict)
    except:
        L.warning('No neuronDB.xml file found, use same mtype for all')
        return load_morphologies_from_folder(morph_path)



 
def save_neuron(neuron, model, folder):
        """ save the neuron morphology for later analysis """

        if not os.path.isdir(folder):
            os.mkdir(folder)

        neuron[0].write(folder + '/' + model + '_' + neuron[1] + '.asc')
  
def load_neuron(fname, model, folder):
        """ load the neuron morphology for later analysis """
        if model:
            return nm.load_neuron(folder + '/' + model + '_' + fname + '.asc')
        else:
            return nm.load_neuron(folder + '/' + fname + '.asc')
           

 
#################################
## diameter handling functions ##
#################################

def set_diameters(section, diameters):
    """hack to set diameters with neurom"""

    new_points = section.points
    new_points[:, COLS.R] = diameters/2.
    section.points = new_points 

def get_mean_diameter(section):
    """hack to get diameters with neurom"""

    vecs = np.diff(section.points, axis=0)[:, COLS.XYZ]
    lengths = [np.sqrt(np.dot(p, p)) for p in vecs]

    segment_mean_diams = section.points[1:, COLS.R] + section.points[:-1, COLS.R]
    mean_diam = np.sum(segment_mean_diams*lengths)/np.sum(lengths)
     
    return mean_diam


def get_diameters(section):
    """hack to get diameters with neurom"""

    return section.points[:, COLS.R]*2 

def redefine_diameter_section(section, diam_ind, diam_new):
    """Hack to replace one diameter at index diam_ind with value diam_new"""

    diameters = get_diameters(section)
    diameters[diam_ind] = diam_new
    set_diameters(section, diameters)      

##########################
## additional functions ##
##########################

def section_lengths(section):
    """Computes all segment lengths within section"""

    vecs = np.diff(section.points, axis=0)[:, COLS.XYZ]
    d2 = [np.dot(p, p) for p in vecs]
    return list(np.cumsum(np.sqrt(d2)))

def tqdm_disable(morphologies):
    """ to have a single progression bar """

    if len(morphologies) >1:
        tqdm_1 = False
        tqdm_2 = True
    else:
        tqdm_1 = True
        tqdm_2 = False

    return tqdm_1, tqdm_2
 
def set_bins(data, n_bins, n_min = 20):
    """ find a good set of bins to avoid undersampling """

    #try to bin uniformly
    max_data = np.max(data)
    min_data = np.min(data)
    diff_max = max_data - min_data

    values, bins = np.histogram(data, bins = n_bins, range = (min_data, max_data) )

    def reduce_bounds(values, bins, data):
        max_data = np.max(data)
        min_data = np.min(data)

        #if the last bins have to few points, reduce the window
        while values[-1] < n_min and max_data-min_data>0.1*diff_max: #second condition is to prevent shrinking for ever
            max_data = bins[-2]
            values, bins = np.histogram(data, bins = n_bins, range = (min_data, max_data) )

        #if the first bins have to few points, reduce the window
        while values[0] < n_min and max_data-min_data>0.1*diff_max:
            min_data = bins[1]
            values, bins = np.histogram(data, bins = n_bins, range = (min_data, max_data) )

        return values, bins

    #find new bounds
    values, bins = reduce_bounds(values, bins, data)

    #if bins have to few elements, reduce the number of bins and readjust bounds
    while len(values[values<n_min])>0 and n_bins>2:
        n_bins -= 1
        max_data = np.max(data)
        min_data = np.min(data)
        values, bins = np.histogram(data, bins = n_bins, range = (min_data, max_data) )

        values, bins = reduce_bounds(values, bins, data)

    return bins, values

