import json 


extract_models_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'by_mtypes': False, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0','M1', 'M2', 'M3', 'M4'],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params_all.json',
    'fig_folder': 'figures_all', 
    'ext': '.png',
    'extra_params': {
            'M0': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M1': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M2': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M3': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M4': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
    }
} 


with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '../scripts/extract_morphologies/selected_morphologies/',
    'by_mtypes': False,
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0', 'M1', 'M2', 'M3', 'M4'],
    'neurite_types': ['basal'],
    'models_params_file': 'model_params_all.json',
    'new_morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    'extra_params': {
            'M0': {
                'seed': 1, 
                'max_bo_fit': 10,
                'trunk_floc': 0, 
                'taper': 0},
            'M1': {
                'seed': 10, 
                'max_bo_fit': 10,
                'trunk_floc': 0, 
                'taper': 0},
            'M2': {
                'seed': 100, 
                'max_bo_fit': 10,
                'trunk_floc': 0, 
                'taper': 0},
            'M3': {
                'seed': 1000, 
                'max_bo_fit': 10,
                'trunk_floc': 0, 
                'taper': 0},
            'M4': {
                'seed': 10000, 
                'max_bo_fit': 10,
                'trunk_floc': 0, 
                'taper': 0}
    }
} 

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


reextract_models_params = {
    'morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    'by_mtypes': False, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0','M1', 'M2', 'M3'],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params_remodel.json',
    'fig_folder': 'figures_remodel', 
    'ext': '.png',
    'extra_params': {
            'M0': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M1': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M2': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M3': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M4': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
    }
} 


with open('reextract_models_params.json', 'w') as json_file:
    json.dump(reextract_models_params, json_file, sort_keys=True, indent=4)


