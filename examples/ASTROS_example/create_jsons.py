import json 
import numpy as np
np.random.seed(1)

extract_models_params = {
    'morph_path': '/home/arnaudon/Dropbox/BBP/data_tests/astros',
    'mtypes_sort': 'all', 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0',], 
    'neurite_types': ['basal', 'axon',],
    'models_params_file': 'model_params.json',
    'plot': True,
    'fig_folder': 'figures', 
    'ext': '.png',
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 1., 
                'trunk_min_sample_num': {'basal': 10, 'apical': 10},
                'trunk_floc': None, 
                'orders': {'a': 0, 'loc':3, 'scale':3, 'min':3, 'max':2},
                'taper': {'max_residual': 1000, 'zeros':-1e-8, 'max': 100.005, 'min': -100.0025}
                }

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '/home/arnaudon/Dropbox/BBP/data_tests/astros',
    'mtypes_sort': 'all',
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0'],
    'neurite_types': ['basal', 'axon',],
    'models_params_file': 'model_params.json',
    'new_morph_path': 'new_morphologies/', 
    #'new_morph_path': '../scripts/extract_morphometrics/new_morphologies_super_mtypes/', 
    'plot': True, 
    'n_cpu': 10, 
    'extra_params': {}
}

for i, model in enumerate(generate_diameters_params['models']):
    generate_diameters_params['extra_params'][model] = {
                'seed': i 
               }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)