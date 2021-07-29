from xml.etree import ElementTree as ET
import os
import shutil
from tqdm import tqdm

#diable warnings
import morphio
morphio.set_maximum_warnings(0)


import neurom as nm

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from neurom.view import matplotlib_impl

def plot_neurons_per_mtype_from_xml(folder = '.', filename='neuronDB.xml', output_dir='Subtypes/', ext = '.png'):
    '''
    Generates directories of mtypes that plot all cells
    of this type, based on the data in filename
    Filename should be of neuronDB.xml type
    '''
    FileDB = ET.parse(folder+filename)
    root = FileDB.findall('listing')[0]
    morphs = root.findall('morphology')

    mtypes = set()
    for m in morphs:
        try:
            # Define mtypes
            mtype = m.find('mtype').text
            # Define suntypes (if they exist)
            if m.find('msubtype').text:
                mtype = mtype + ':' + m.find('msubtype').text
            mtypes.add(mtype)
        except:
            print('Failed to process', m)

    # Create a directory for each m-type, subtype
    for m in mtypes:
        os.mkdir(output_dir + m)

    for m in tqdm(morphs):
        mtype = output_dir + m.find('mtype').text
        if m.find('msubtype').text:
            mtype = mtype + ':' + m.find('msubtype').text
        # Copy all morphologies to the corresponding directory according to mtypes

        morph = nm.load_morphology(folder + m.find('name').text + '.asc')
        matplotlib_impl.plot_morph(morph)
        plt.savefig(mtype + '/' + m.find('name').text + ext, dpi = 500)
        plt.close()

if __name__ == '__main__':

    if not os.path.isdir('neuron_plots'):
        os.mkdir('neuron_plots')
    else:
        shutil.rmtree('neuron_plots/')
        os.mkdir('neuron_plots')

    plot_neurons_per_mtype_from_xml(folder='../../examples/05_RepairUnravel-asc/', output_dir='neuron_plots/', ext = '.png')
