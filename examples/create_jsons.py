import json 


extract_models_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'by_mtypes': False, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0', ],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params_all.json',
    'fig_folder': 'figures_all', 
    'ext': '.png'
} 

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '../scripts/extract_morphologies/selected_morphologies/',
    'by_mtypes': False,
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0',],
    'neurite_types': ['basal'],
    'models_params_file': 'model_params_all.json',
    'new_morph_path': '../scripts/diameter-checks/new_morphologies_all/', 
} 

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


