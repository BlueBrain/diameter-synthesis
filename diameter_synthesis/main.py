import os, glob, shutil
import click
import json
from tqdm import tqdm
import numpy as np

#diable warnings
import morphio
morphio.set_maximum_warnings(0)

from diameter_synthesis.build_models import build_models
from diameter_synthesis.build_diameters import build_diameters
import diameter_synthesis.utils as utils 

#hack to convert numpy types to python types for json compatibility
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

@click.group()
def cli():
    pass

@cli.command('run_models')
@click.argument('config_file', type=click.Path(exists=True))
def run_models(config_file):
    """ Run the model extraction from config file"""
    
    #load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)

    print('Loading morphologies...')
    #load all the morphologies 
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'])

    print('Extracting models parameters...')
    #compute the model
    models_params = build_models(config['models'], morphologies, config['neurite_types'], config['extra_params'], fig_folder = config['fig_folder'], ext = config['ext'], plot = config['plot'])

    #save the models parameters in a json file
    with open(config['models_params_file'], 'w') as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NpEncoder)

@cli.command('run_diameters')
@click.argument('config_file', type=click.Path(exists=True))
def run_diameters(config_file):
    """ Build new diameters from config file and diameter model"""
    
    #load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)

    print('Loading morphologies...')
    #load all the morphologies 
    morphologies_dict = utils.create_morphologies_dict(config['morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'])

    print('Generate diameters...')
    with open(config['models_params_file'], 'r') as f:
        models_params = json.load(f)

    #generate diameters
    models_params = build_diameters(config['models'], models_params , morphologies_dict, config['neurite_types'], config['new_morph_path'], config['extra_params'], config['morph_path'], plot = config['plot'], n_cpu = config['n_cpu'])
