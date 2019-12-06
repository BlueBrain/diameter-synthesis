import neurom as nm
import os, glob, shutil

# diable warnings
import morphio
morphio.set_maximum_warnings(0)

import diameter_synthesis.types as tpes 

"""
script to split cells in two types: pyramidal cells, and interneuron cells
"""

#input_folder = '../../examples/05_RepairUnravel-asc/'
input_folder = '../extract_morphologies/selected_morphologies/'
PC_folder = 'PC_neurons_selected'
Inter_folder = 'Inter_neurons_selected'
if not os.path.isdir(PC_folder):
    os.mkdir(PC_folder)
if not os.path.isdir(Inter_folder):
    os.mkdir(Inter_folder)

for fname in os.listdir(input_folder):

    if fname.endswith(('.h5', '.asc', '.swc')):
        neuron = nm.load_neuron(input_folder+fname)
        IN = True
        for neurite in neuron.neurites:
            if neurite.type == tpes.STR_TO_TYPES['apical']:
                shutil.copyfile(input_folder + fname, PC_folder +'/' + fname)
                print(fname, '==> PC cell')
                IN = False
                break
        if IN:
            print(fname, '==> IN cell')
            shutil.copyfile(input_folder + fname, Inter_folder + '/' + fname)
