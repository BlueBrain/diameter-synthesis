import os, glob, shutil
import click
import json
from tqdm import tqdm

#diable warnings
import morphio
morphio.set_maximum_warnings(0)

import diameter_synthesis.extract_models as extract_models
import diameter_synthesis.generate_diameters as generate_diameters
import diameter_synthesis.utils as utils 

@click.group()
def cli():
    pass

@cli.command('build_models')
@click.argument('config_file', type=click.Path(exists=True))
def build_models(config_file):
    """ Run the model extraction from config file"""
    
    #load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)

    print('Loading morphologies...')
    #load all the morphologies 
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], by_mtypes = config['by_mtypes'], n_mtypes_max = config['n_mtypes_max'])

    print('Extracting models parameters...')
    #compute the model
    models_params = extract_models.build_models(config['models'], morphologies, config['neurite_types'], fig_folder = config['fig_folder'], ext = config['ext'])

    #save the models parameters in a json file
    with open(config['models_params_file'], 'w') as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4)

@cli.command('build_diameters')
@click.argument('config_file', type=click.Path(exists=True))
def build_diameters(config_file):
    """ Build new diameters from config file and diameter model"""
    
    #load the config file
    with open(config_file, 'r') as f:
        config = json.load(f)


    print('Loading morphologies...')
    #load all the morphologies 
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], by_mtypes = config['by_mtypes'], n_mtypes_max = config['n_mtypes_max'])

    print('Generate diameters...')
    with open(config['models_params_file'], 'r') as f:
        models_params = json.load(f)

    #generate diameters
    models_params = generate_diameters.build_diameters(config['models'], models_params , morphologies, config['neurite_types'], config['new_morph_path'])
