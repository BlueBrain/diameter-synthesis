import json 


extract_models_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'neurite_types': ['basal', 'apical'],
    'models': ['M0', ],
    'by_mtypes': True, 
    'n_morphs_max': None, 
    'n_mtypes_max': 60,
    'models_params_file': 'model_params.json',
    'fig_folder': 'figures', 
    'ext': '.png'
} 

with open('extract_models_params.json', 'w') as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)


generate_diameters_params = {
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj81/InputData/2017Release/OriginalData/05_RepairUnravel-asc/', 
    'models': ['M0',],
    'models_params_file': 'model_params.json',
} 

with open('generate_diameters_params.json', 'w') as json_file:
    json.dump(generate_diameters_params, json_file, sort_keys=True, indent=4)


