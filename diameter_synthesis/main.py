"""Main functions to learn and generate diameters."""
import json
import logging
import multiprocessing
import os
from pathlib import Path
from tqdm import tqdm
import numpy as np

from morphio.mut import Morphology

import neurom as nm

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


def plot_models(morphologies, config, models_params, models_data, ext="png"):
    """Plot the models."""
    L.info("Plot the fits...")

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


def _build_all_models(morphologies, config, plot=False, ext="png"):
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


def run_models(config_file, plot, ext="png"):
    """Run the model extraction from config file."""
    with open(config_file, "r") as filename:
        config = json.load(filename)

    L.info("Loading morphologies...")
    morphologies_dict = utils.create_morphologies_dict(
        config["morph_path"],
        mtypes_file=config["mtypes_file"],
    )

    morphologies = {
        mtype: [nm.load_morphology(i) for i in morphologies_dict[mtype]]
        for mtype in tqdm(morphologies_dict)
    }

    L.info("Extracting model parameters...")
    models_params = _build_all_models(morphologies, config, plot=plot, ext=ext)

    with open(config["models_params_file"], "w") as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NumpyEncoder)


class DiameterWorker:
    """Worker for building diameters."""

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
            self.models_params[mtype],
            self.config["neurite_types"],
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
    """Build new diameters from config file and diameter model."""
    with open(config_file, "r") as filename:
        config = json.load(filename)

    with open(models_params_file, "r") as filename:
        models_params = json.load(filename)

    for model in config["models"]:
        L.info("Generating diameter with model %s", model)

        morphologies_dict = utils.create_morphologies_dict(
            config["morph_path"],
            mtypes_file=config["mtypes_file"],
        )

        worker = DiameterWorker(model, models_params, config)

        all_neurons = [
            [neuron, mtype] for mtype in morphologies_dict for neuron in morphologies_dict[mtype]
        ]

        pool = multiprocessing.Pool(config["n_cpu"])  # pylint: disable=consider-using-with
        list(tqdm(pool.imap(worker, all_neurons), total=len(all_neurons)))
        pool.close()
        pool.join()
