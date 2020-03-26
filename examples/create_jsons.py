import json 
import numpy as np
np.random.seed(1)

extract_models_params = {
    #'morph_path': './morphologies/',
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-h5-sanitized/',
    'mtypes_file':  '/gpfs/bbp.cscs.ch/project/proj82/issues/BBPP82-318/neurondb_05_RepairUnravel.dat',
    'models': ['generic',], 
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_mtypes.json',
    'fig_folder': 'model_figures', 
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 2., 
                'taper': {'max': 1e-6, 'min': -0.01},
                'threshold': {'apical': 0.2, 'basal': 1.}
                }


with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

generate_diameters_params = {}
for i, model in enumerate(extract_models_params['models']):
    generate_diameters_params[model] = {
        'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-h5-sanitized/',
        #'morph_path': './morphologies/',
        'mtypes_file': '/gpfs/bbp.cscs.ch/project/proj82/issues/BBPP82-318/neurondb_05_RepairUnravel.dat',
        'neurite_types': ['basal', 'apical'],
        'models_params_file': 'model_params.json',
        'new_morph_path': './diametrized_morphologies/', 
        'n_cpu': 1, 
        'n_samples': 2,
        'seed': 0,
        'trunk_max_tries': 100,
        'asymetry_threshold': {'apical': 0.2, 'basal': 1.}
    }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)
