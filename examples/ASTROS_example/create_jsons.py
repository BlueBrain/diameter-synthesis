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
    'morph_path': 'astros',
    'mtypes_sort': 'all',
    'n_morphs_max': None,
    'n_mtypes_max': 60,
    'models': ['M1',],
    'neurite_types': ['basal', 'axon'],
    'models_params_file': 'model_params.json',
    'plot': True,
    'fig_folder': 'figures',
    'ext': '.png',
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 2.,
                'taper': {'max_residual': 1000, 'zeros':-1e-8, 'max': 100.005, 'min': -100.0025},
                'threshold': {'basal': 1., 'axon':1.}
                }

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': 'astros',
    'mtypes_sort': 'all',
    'n_morphs_max': None,
    'n_mtypes_max': 60,
    'models': ['M1'],
    'neurite_types': ['basal', 'axon'],
    'models_params_file': 'model_params.json',
    'new_morph_path': 'new_morphologies/',
    #'new_morph_path': '../scripts/extract_morphometrics/new_morphologies_super_mtypes/',
    'plot': True,
    'n_cpu': 3,
    'n_samples': 1,
    'ext': '.png',
    'extra_params': {}
}

for i, model in enumerate(generate_diameters_params['models']):
    generate_diameters_params['extra_params'][model] = {
                'seed': i ,
                'trunk_max_tries': 1,
                'threshold': {'basal': 1., 'axon':1.}
               }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)
