""" main functions to learn and generate diameters """
import json
import logging
import os

import morphio
import numpy as np

import diameter_synthesis.utils as utils
from diameter_synthesis.build_diameters import build_diameters
from diameter_synthesis.build_models import build_models

morphio.set_maximum_warnings(0)
L = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


class NumpyEncoder(json.JSONEncoder):
    '''To encode numpy arrays'''
    def default(self, o):  # pylint: disable=method-hidden
        '''encoder'''
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.integer):
            return int(o)
        return json.JSONEncoder.default(self, o)


def run_models(config_file):
    """ Run the model extraction from config file"""

    # load the config file
    with open(config_file, 'r') as filename:
        config = json.load(filename)

    L.info('Loading morphologies...')
    # load all the morphologies
    morphologies = utils.load_morphologies(
        config['morph_path'],
        mtypes_sort=config['mtypes_sort'],
        mtypes_file=config['mtypes_file'])

    L.info('Extracting model parameters...')
    # compute the model
    models_params = build_models(morphologies, config)

    # save the models parameters in a json file
    with open(config['models_params_file'], 'w') as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NumpyEncoder)


def run_diameters(config_file):
    """ Build new diameters from config file and diameter model"""

    # load the config file
    with open(config_file, 'r') as filename:
        config = json.load(filename)

    L.info('Loading morphologies...')
    # load all the morphologies
    morphologies_dict = utils.create_morphologies_dict(
        config['morph_path'],
        mtypes_sort=config['mtypes_sort'],
        mtypes_file=config['mtypes_file'])

    L.info('Generate diameters...')
    with open(config['models_params_file'], 'r') as filename:
        models_params = json.load(filename)

    # generate diameters
    build_diameters(morphologies_dict, models_params, config)
