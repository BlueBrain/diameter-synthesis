import json 
import numpy as np
np.random.seed(1)

m_sort = 'mtypes'

extract_models_params = {
    'morph_path': './morphologies/',
    'mtypes_file': '/gpfs/bbp.cscs.ch/project/proj82/issues/BBPP82-318/neurondb_05_RepairUnravel.dat',
    'mtypes_sort': m_sort, 
    'models': ['generic',], 
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_mtypes.json',
    'fig_folder': 'figures_' + m_sort, 
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 2., 
                #'taper': {'max_residual': 100, 'zeros': 1e-8, 'max': 0.005, 'min': -0.010},
                'taper': {'max': 1e-6, 'min': -0.01},
                #'taper': {'max_residual': 100, 'zeros': 1e-8, 'max': -0.002, 'min': -0.002},
                'threshold': {'apical': 0.15, 'basal': 1.}
                }


with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

generate_diameters_params = {}
for i, model in enumerate(extract_models_params['models']):
    generate_diameters_params[model] = {
        #'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-h5/',
        'morph_path': './morphologies/',
        'mtypes_file': '/gpfs/bbp.cscs.ch/project/proj82/issues/BBPP82-318/neurondb_05_RepairUnravel.dat',
        'mtypes_sort': m_sort,
        'neurite_types': ['basal', 'apical'],
        'models_params_file': 'model_params_' + m_sort + '.json',
        'new_morph_path': './diametrized_morphologies/', 
        #'plot': False, 
        'n_cpu': 10, 
        'n_samples': 2,
        'seed': 0,
        'trunk_max_tries': 100,
        'asymetry_threshold': {'apical': 0.15, 'basal': 1.}
    }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)
