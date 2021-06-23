"""Create the JSON file use as example.

Copyright (C) 2021  Blue Brain Project, EPFL

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import json
import numpy as np

np.random.seed(1)

extract_models_params = {
    'morph_path': './morphologies/',
    'mtypes_file':  None,
    #'morph_path': '/gpfs/bbp.cscs.ch/project/proj82/simulations/cell_synthesis_distributions/2017_repaired_cells/',
    #'mtypes_file': '/gpfs/bbp.cscs.ch/project/proj82/simulations/cell_synthesis_distributions/dat_files_correction/neurondb_05_RepairUnravel.dat',
    'new_morph_path': './diametrized_morphologies/',
    'models': ['astrocyte'],
    'neurite_types': ['basal', 'axon'],
    'models_params_file': 'model_params_mtypes.json',
    'fig_folder': 'model_figures',
    'terminal_threshold': 2.,
    'taper': {'max': 1.1, 'min': -1.1},
    'asymmetry_threshold': {'axon': 0.05, 'basal': 0.05},
    'n_samples': 1,
    'seed': 0,
    'trunk_max_tries': 100,
    'n_cpu': 10,
}

with open("diametrizer_params.json", "w") as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)
