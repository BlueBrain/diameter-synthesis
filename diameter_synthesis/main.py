import diameter_synthesis.utils as utils
from diameter_synthesis.build_diameters import build_diameters
from diameter_synthesis.build_models import build_models
import os
import glob
import shutil
import json
from tqdm import tqdm
import numpy as np

# diable warnings
import morphio
morphio.set_maximum_warnings(0)


# hack to convert numpy types to python types for json compatibility


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

def run_models(config_file):
    """ Run the model extraction from config file"""

    # load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)

    print('Loading morphologies...')
    # load all the morphologies
    morphologies_tmp = utils.load_morphologies(
        config['morph_path'], n_morphs_max=config['n_morphs_max'], mtypes_sort=config['mtypes_sort'], n_mtypes_max=config['n_mtypes_max'])

    #remove empty mtypes if missing data (not clean...)
    morphologies = morphologies_tmp.copy() 
    for m in morphologies_tmp:
        if len(morphologies_tmp[m]) == 0:
            del morphologies[m]

    print('Extracting model parameters...')
    # compute the model
    models_params = build_models(morphologies, config)

    # save the models parameters in a json file
    with open(config['models_params_file'], 'w') as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NpEncoder)


def run_diameters(config_file):
    """ Build new diameters from config file and diameter model"""

    # load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)

    print('Loading morphologies...')
    # load all the morphologies
    morphologies_dict = utils.create_morphologies_dict(
        config['morph_path'], n_morphs_max=config['n_morphs_max'], mtypes_sort=config['mtypes_sort'], n_mtypes_max=config['n_mtypes_max'])

    print('Generate diameters...')
    with open(config['models_params_file'], 'r') as f:
        models_params = json.load(f)

    # generate diamet[ers
    build_diameters(morphologies_dict, models_params, config)
            
