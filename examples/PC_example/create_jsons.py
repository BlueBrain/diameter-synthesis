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

m_sort = 'mtypes'

extract_models_params = {
    'morph_path': 'PC_neurons/',
    #'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/',
    'mtypes_sort': m_sort,
    'models': ['generic',],
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_' + m_sort + '.json',
    'plot': False,
    'fig_folder': 'figures_' + m_sort,
    'ext': '.svg',
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 2.,
                'taper': {'max_residual': 100, 'zeros':1e-8, 'max': 0.005, 'min': -0.010},
                'threshold': {'apical': 0.2, 'basal': 1.}
                }

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': 'PC_neurons/',
    'mtypes_sort': m_sort,
    'models': ['generic2'],
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_' + m_sort + '.json',
    'new_morph_path': 'new_morphologies_' + m_sort + '/',
    'plot': True,
    'n_cpu': 5,
    'n_samples': 10,
    'ext': '.png',
    'extra_params': {}
}

for i, model in enumerate(generate_diameters_params['models']):
    generate_diameters_params['extra_params'][model] = {
                'seed': 1,
                'trunk_max_tries': 100,
                'threshold': {'apical': 0.2, 'basal': 1.}
               }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)
