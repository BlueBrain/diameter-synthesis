import json 

taper = -0.0005
extract_models_params = {
    'morph_path': '05_RepairUnravel-asc/', 
    #'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'by_mtypes': False, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
    'neurite_types': ['basal'],
    'models_params_file': 'model_params_all.json',
    'fig_folder': 'figures_all', 
    'ext': '.png',
    'extra_params': {
            'M0': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M1': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M2': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M3': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M4': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M5': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M6': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M7': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M8': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
            'M9': {
                'terminal_threshold': 1.2, 
                'trunk_min_sample_num': {'basal': 20, 'apical': 5},
                'trunk_floc': 0},
    }
} 


with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    #'morph_path': '../scripts/extract_morphologies/selected_morphologies/',
    'morph_path': '05_RepairUnravel-asc/', 
    'by_mtypes': False,
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0', 'M1',], # 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9'],
    'neurite_types': ['basal'],
    'models_params_file': 'model_params_all.json',
    #'new_morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    'new_morph_path': '../scripts/extract_morphometrics/new_morphologies_all/', 
    'plot': False, 
    'n_cpu': 10, 
    'extra_params': {
            'M0': {
                'seed': 1, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M1': {
                'seed': 2, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M2': {
                'seed': 3, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M3': {
                'seed': 4, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M4': {
                'seed': 5, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M5': {
                'seed': 6, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M6': {
                'seed': 7, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M7': {
                'seed': 8, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M8': {
                'seed': 9, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
            'M9': {
                'seed': 10, 
                'max_bo_fit': 15,
                'trunk_floc': 0, 
                'taper': taper},
    }
} 

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


reextract_models_params = {
    #'morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
    'morph_path': '../scripts/extract_morphometrics/new_morphologies_all/', 
    'by_mtypes': True, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0'], #,'M1', 'M2', 'M3','M4','M5', 'M6', 'M6', 'M8', 'M9'],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params_remodel_all.json',
    'fig_folder': 'figures_remodel_all', 
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
                'trunk_min_sample_num': {'basal': 9, 'apical': 5},
                'trunk_floc': 0},
            'M4': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M5': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M6': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M7': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M8': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
            'M9': {
                'terminal_threshold': 0.8, 
                'trunk_min_sample_num': {'basal': 8, 'apical': 5},
                'trunk_floc': 0},
    }
} 


with open('reextract_models_params.json', 'w') as json_file:
    json.dump(reextract_models_params, json_file, sort_keys=True, indent=4)


