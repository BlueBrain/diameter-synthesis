import json 
import numpy as np
np.random.seed(1)

extract_models_params = {
    'morph_path': '05_RepairUnravel-asc/', 
    #'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'mtypes_sort': 'mtypes', 
    'n_morphs_max': None, 
    'n_mtypes_max': 10,
    'models': ['M2',], # 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_mtypes.json',
    'plot': True,
    'fig_folder': 'figures_mtypes', 
    'ext': '.png',
    'extra_params': {}

}

for model in extract_models_params['models']:
    extract_models_params['extra_params'][model] = {
                'terminal_threshold': 1., 
                'trunk_min_sample_num': {'basal': 5, 'apical': 5},
                'trunk_floc': None, 
                'orders': {'a': 4, 'loc':4, 'scale':4, 'min':4, 'max':4},
                'taper': {'max_residual': 10, 'zeros':1e-8, 'max': 0.002, 'min': -0.005},
                'threshold': 0.2,
                }

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    #'morph_path': '../scripts/extract_morphologies/selected_morphologies/',
    'morph_path': '05_RepairUnravel-asc/', 
    'mtypes_sort': 'all',
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M2'], #, 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_all.json',
    'new_morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    #'new_morph_path': '../scripts/extract_morphometrics/new_morphologies_super_mtypes/', 
    'plot': True, 
    'n_cpu': 1, 
    'n_samples': 2,
    'extra_params': {}
}

for i, model in enumerate(generate_diameters_params['models']):
    generate_diameters_params['extra_params'][model] = {
                'seed': 10,
                'threshold': 0.2,
               }

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


reextract_models_params = {
    #'morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    'morph_path': '../scripts/extract_morphometrics/new_morphologies_all/', 
    'mtypes_sort': 'all', 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0',],# 'M1', 'M2',], # 'M3','M4','M5', 'M6', 'M6', 'M8', 'M9'],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params_remodel_all.json',
    'fig_folder': 'figures_remodel_all', 
    'ext': '.png',
    'plot': True,
    'extra_params': {
            'M0': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None,
                'orders': {'a': 1, 'loc':2, 'scale':2, 'min':2, 'max':2}
                },
            'M1': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M2': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M3': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 9, 'apical': 5},
                'trunk_floc': None},
            'M4': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M5': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M6': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M7': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M8': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
            'M9': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': None},
    }
} 


with open('reextract_models_params.json', 'w') as json_file:
    json.dump(reextract_models_params, json_file, sort_keys=True, indent=4)


