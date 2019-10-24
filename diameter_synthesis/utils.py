import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm
from neurom import COLS
from neurom.core import iter_sections

from neurom import NeuriteType
from morphio import SectionType

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


##############################
## loading/saving functions ##
##############################

def load_morphologies_from_dict(morph_path, name_dict):
    """ Load the morphologies from a list of files """

    tqdm_1, tqdm_2 = tqdm_disable(name_dict) #to have a single progression bar
    # just to get single progression bar 
    morphologies = {}
    for mtype in tqdm(name_dict, disable = tqdm_1):
        morphologies[mtype] = []
        for fname in tqdm(name_dict[mtype], disable = tqdm_2):
            name, ext = os.path.splitext(fname)
            if ext in {'.h5', '.asc', '.swc'} and os.path.exists(morph_path + '/' + fname):
                neuron = nm.load_neuron(morph_path + '/' + fname)
                morphologies[mtype].append([neuron, name])

    return morphologies 

def load_morphologies(morph_path, by_mtypes = True, n_morphs_max = None, n_mtypes_max = None, xml_file = './neuronDB.xml', ext = '.asc'):
    """ Load the morphologies from a directory, by mtypes or all at once """
    
    #load morphologies by mtypes, one listt of morphologies for each type
    if by_mtypes:

        #if not max number of mtypes, take all
        if not n_mtypes_max:
            n_mtypes_max =  1e10


        from xml.etree import ElementTree as ET
        FileDB = ET.parse(morph_path + xml_file)
        root = FileDB.findall('listing')[0]
        morphs = root.findall('morphology')
    
        name_dict = {}
        for m in morphs:
            try:
                # Define mtypes
                mtype = m.find('mtype').text
                # Define suntypes (if they exist)
                if m.find('msubtype').text:
                    mtype = mtype + ':' + m.find('msubtype').text

                #if it is a new mtype, add en entry to name_dict
                if mtype not in name_dict.keys() and len(name_dict) < n_mtypes_max:
                    name_dict[mtype] = [m.find('name').text + ext]
                elif mtype in name_dict.keys():
                    name_dict[mtype] += [m.find('name').text + ext]

            except:
                print('Failed to process', m)

    #load all the morphologies together
    else:
        name_dict = {}
        if n_morphs_max is not None:
            name_dict['all_types'] = os.listdir(morph_path)[:n_morphs_max]
        else:
            name_dict['all_types'] = os.listdir(morph_path)
    
    return load_morphologies_from_dict(morph_path, name_dict)
 
def save_neuron(neuron, model, folder):
        """ save the neuron morphology for later analysis """

        if not os.path.isdir(folder):
                os.mkdir(folder)

        neuron[0].write(folder + '/' + model + '_' + neuron[1] + '.asc')

 
#################################
## diameter handling functions ##
#################################

def set_diameters(section, diameters):
    """hack to set diameters with neurom"""

    new_points = section.points
    new_points[:, COLS.R] = diameters/2.
    section.points = new_points 

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
    return np.cumsum(np.sqrt(d2))

def tqdm_disable(morphologies):
    """ to have a single progression bar """

    if len(morphologies) >1:
        tqdm_1 = False
        tqdm_2 = True
    else:
        tqdm_1 = True
        tqdm_2 = False

    return tqdm_1, tqdm_2

