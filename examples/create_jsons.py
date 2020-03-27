import json
import numpy as np

np.random.seed(1)

extract_models_params = {
    #'morph_path': './morphologies/',
    #'mtypes_file':  None, 
    'morph_path': '/gpfs/bbp.cscs.ch/project/proj82/simulations/cell_synthesis_distributions/2017_repaired_cells/',
    'mtypes_file': '/gpfs/bbp.cscs.ch/project/proj82/simulations/cell_synthesis_distributions/dat_files_correction/neurondb_05_RepairUnravel.dat',
    'new_morph_path': './diametrized_morphologies/', 
    'models': ['generic'],
    'neurite_types': ['basal', 'apical'],
    'models_params_file': 'model_params_mtypes.json',
    'fig_folder': 'model_figures', 
    'terminal_threshold': 2., 
    'taper': {'max': 1e-6, 'min': -0.01},
    'asymetry_threshold': {'apical': 0.2, 'basal': 1.},
    'n_samples': 2,
    'seed': 0,
    'trunk_max_tries': 100,
    'n_cpu': 1,
}

with open("diametrizer_params.json", "w") as json_file:
    json.dump(extract_models_params, json_file, sort_keys=True, indent=4)
