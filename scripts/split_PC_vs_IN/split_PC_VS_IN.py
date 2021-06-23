"""Script to split cells in two types: pyramidal cells, and interneuron cells."""

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

import neurom as nm
import os, glob, shutil

# diable warnings
import morphio
morphio.set_maximum_warnings(0)

import diameter_synthesis.types as tpes


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
        neuron = nm.load_morphology(input_folder+fname)
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
