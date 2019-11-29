import json 
import numpy as np
np.random.seed(1)

extract_models_params = {
    'morph_path': '/home/arnaudon/Dropbox/BBP/data_tests/L5_TPC:A',
    'mtypes_sort': 'all', 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M2',], 
    'neurite_types': ['apical',],
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
                'orders': {'a': 4, 'loc':4, 'scale':4, 'min':4, 'max':4},
                'taper': {'max_residual': 10, 'zeros':1e-8, 'max': 0.002, 'min': -0.005}
                }

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '/home/arnaudon/Dropbox/BBP/data_tests/L5_TPC:A',
    'mtypes_sort': 'all',
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M2'],
    'neurite_types': ['apical',],
    'models_params_file': 'model_params.json',
    'new_morph_path': 'new_morphologies/', 
    #'new_morph_path': '../scripts/extract_morphometrics/new_morphologies_super_mtypes/', 
    'plot': True, 
    'n_cpu': 1, 
    'n_samples': 10,
    'extra_params': {}
}

for i, model in enumerate(generate_diameters_params['models']):
    generate_diameters_params['extra_params'][model] = {
                'seed': 10,
                'sibling_threshold': 0.2,
                'rall_threshold': 0.2
               }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)
