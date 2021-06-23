"""Script find GIDs."""

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

import numpy as np
import matplotlib.pyplot as plt

import bglibpy
import neurom as nrm
from bglibpy import bluepy
import os, re, json

from utils import *

bglibpy.set_verbose(100)

morph_dir_bio = '../../examples/05_RepairUnravel-asc/'
morph_dir_diam = 'new_neurons/'

mtypes=['L5_TPC:A',]

list_cells_full = os.listdir(morph_dir_diam)
list_cells = []

for cell in list_cells_full:
    list_cells.append( os.path.splitext(cell)[0])

ssim = initialize_sim(morph_dir_bio)
gids = get_cells_gid(ssim, list_cells, mtypes=mtypes)

print('Found gids for ', len(gids), 'cells')

with open('gids.json', 'w') as json_file:
    json.dump(gids, json_file, sort_keys=True, indent=4)





