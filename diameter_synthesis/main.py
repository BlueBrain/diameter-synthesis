"""Main functions to learn and generate diameters."""

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import json
import logging
import multiprocessing
import os
from pathlib import Path

import neurom as nm
import numpy as np
from morphio.mut import Morphology
from tqdm import tqdm

from diameter_synthesis import utils
from diameter_synthesis.build_diameters import build as build_diameters
from diameter_synthesis.build_models import build as build_model
from diameter_synthesis.plotting import plot_distribution_fit

L = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


class NumpyEncoder(json.JSONEncoder):
    """Class to encode numpy arrays."""

    def default(self, o):  # pylint: disable=method-hidden
        """Numpy encoder."""
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.integer):
            return int(o)
        return json.JSONEncoder.default(self, o)


def plot_models(morphologies, config, models_params, models_data, ext=".png"):
    """Plot the models.

    Args:
        morphologies (dict): a dict with mtype->[morphologies].
        config (dict): the config to use.
        models_params (dict): the models parameters.
        models_data (dict): the models data.
        ext (str): the file extension used to export the figures.
    """
    L.info("Plot the fits...")

    if not ext.startswith("."):
        ext = f".{ext}"

    if not Path(config["fig_folder"]).exists():
        os.mkdir(config["fig_folder"])

    for mtype in tqdm(morphologies):
        if not (Path(config["fig_folder"]) / mtype).exists():
            os.mkdir(Path(config["fig_folder"]) / mtype)
        for model in config["models"]:
            fit_tpes = models_data[model][mtype]
            model_data = models_data[model][mtype]
            model_param = models_params[model][mtype]

            for fit_tpe in fit_tpes:
                fig_name = Path(config["fig_folder"]) / mtype / fit_tpe
                plot_distribution_fit(
                    model_data[fit_tpe],
                    model_param[fit_tpe],
                    config["neurite_types"],
                    fig_name=fig_name,
                    ext=ext,
                )


def _build_all_models(morphologies, config, plot=False, ext=".png"):
    """Build all the models in the list of models."""
    models_params = {}
    models_data = {}
    for model in config["models"]:
        models_params[model] = {}
        models_data[model] = {}
        for mtype in tqdm(morphologies):
            config_model = config.copy()
            config_model["models"] = [model]
            models_params[model][mtype], models_data[model][mtype] = build_model(
                morphologies[mtype], config_model, with_data=True
            )

    if plot:
        plot_models(morphologies, config, models_params, models_data, ext=ext)

    return models_params


def run_models(config_file, plot, ext=".png"):
    """Run the model extraction from config file.

    Args:
        config_file (str): the path to the configuration file.
        plot (bool): plot the models once they are built.
        ext (str): the file extension used to export the plots.
    """
    with open(config_file, "r", encoding="utf-8") as filename:
        config = json.load(filename)

    L.info("Loading morphologies...")
    morphologies_dict = utils.create_morphologies_dict(
        config["morph_path"],
        mtypes_file=config.get("mtypes_file", None),
    )

    morphologies = {
        mtype: [nm.load_morphology(i) for i in morphologies_dict[mtype]]
        for mtype in tqdm(morphologies_dict)
    }

    L.info("Extracting model parameters...")
    models_params = _build_all_models(morphologies, config, plot=plot, ext=ext)

    with open(config["models_params_file"], "w", encoding="utf-8") as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NumpyEncoder)


class DiameterWorker:
    """Worker for building diameters.

    Args:
        model (str): The model to use.
        models_params (dict): The parameters of the models containing the model name as key.
        config (dict): The configuration to use.
    """

    def __init__(self, model, models_params, config):
        """Init function."""
        self.model = model
        self.models_params = models_params[model]
        self.config = config

    def __call__(self, neuron_input):
        """Call function."""
        fname = neuron_input[0]
        mtype = neuron_input[1]

        neuron = Morphology(fname)

        build_diameters(
            neuron,
            self.config["neurite_types"],
            self.models_params[mtype],
            self.config,
        )

        if not Path(self.config["new_morph_path"]).exists():
            os.mkdir(self.config["new_morph_path"])

        if fname.parts[-2] == mtype:
            if not (Path(self.config["new_morph_path"]) / fname.parts[-2]).exists():
                os.mkdir(Path(self.config["new_morph_path"]) / fname.parts[-2])
            neuron.write(Path(self.config["new_morph_path"]) / fname.parts[-2] / fname.name)
        neuron.write(Path(self.config["new_morph_path"]) / fname.name)


def run_diameters(config_file, models_params_file):
    """Build new diameters from config file and diameter model.

    Args:
        config_file (str): the path to the configuration file.
        models_params_file (str): the path to the file containing the model parameters.
    """
    with open(config_file, "r", encoding="utf-8") as filename:
        config = json.load(filename)

    with open(models_params_file, "r", encoding="utf-8") as filename:
        models_params = json.load(filename)

    for model in config["models"]:
        L.info("Generating diameter with model %s", model)

        morphologies_dict = utils.create_morphologies_dict(
            config["morph_path"],
            mtypes_file=config.get("mtypes_file", None),
        )

        worker = DiameterWorker(model, models_params, config)

        all_neurons = [
            [neuron, mtype] for mtype in morphologies_dict for neuron in morphologies_dict[mtype]
        ]

        with multiprocessing.Pool(config.get("n_cpu", 1)) as pool:
            list(tqdm(pool.imap(worker, all_neurons), total=len(all_neurons)))


def diametrize_single_neuron(neuron, config=None, apical_point_sec_ids=None):
    """Diametrize single neuron by learning diameter model from it.

    Args:
        neuron (mophio.mut.Morphology): neuron to consider.
        config (dict): dict with entry 'model' and 'diameters' with corresponding dicts, if None,
            default dict will be used.
        apical_point_sec_ids (list): list of apical points if any.
    """
    if config is None:
        config = {
            "model": {
                "taper": {"max": 1e-06, "min": -0.1},
                "terminal_threshold": 2.0,
                "models": ["generic"],
                "neurite_types": ["basal", "apical"],
            },
            "diameters": {
                "models": ["generic"],
                "n_samples": 1,
                "neurite_types": ["basal", "apical"],
                "seed": 0,
                "trunk_max_tries": 100,
            },
        }

    model_params = build_model([nm.load_neuron(neuron)], config["model"])
    if apical_point_sec_ids is not None:
        model_params["apical_point_sec_ids"] = apical_point_sec_ids
    build_diameters(neuron, ["basal", "apical"], model_params, config["diameters"])
