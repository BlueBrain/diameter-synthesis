import json 


extract_models_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'by_mtypes': False, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models': ['M0', ],
    'neurite_types': ['basal',],
    'models_params_file': 'model_params.json',
    'fig_folder': 'figures', 
    'ext': '.png'
} 

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'by_mtypes': False,
    'n_morphs_max': 50, 
    'n_mtypes_max': 2,
    'models': ['M0',],
    'neurite_types': ['basal'],
    'models_params_file': 'model_params.json',
    'new_morph_path': 'new_diameters/', 
} 

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


