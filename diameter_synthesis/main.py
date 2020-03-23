""" main functions to learn and generate diameters """
import json
import logging
import multiprocessing
import os

from morphio import RawDataError, UnknownFileType
from morphio.mut import Morphology
from tqdm import tqdm

import diameter_synthesis.utils as utils
from diameter_synthesis import io
from diameter_synthesis.build_diameters import build as build_diameters
from diameter_synthesis.build_models import build as build_model
from diameter_synthesis.types import NumpyEncoder

L = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


def plot_models(morphologies, config, models_params, models_data, ext="png"):
    """plot the models"""
    import diameter_synthesis.plotting as plot  # pylint: disable=import-outside-toplevel

    L.info("Plot the fits...")

    if not os.path.isdir(config["fig_folder"]):
        os.mkdir(config["fig_folder"])

    for mtype in tqdm(morphologies):
        if not os.path.isdir(os.path.join(config["fig_folder"], mtype)):
            os.mkdir(os.path.join(config["fig_folder"], mtype))
        for model in config["models"]:
            fit_tpes = models_data[model][mtype]
            model_data = models_data[model][mtype]
            model_param = models_params[model][mtype]

            for fit_tpe in fit_tpes:
                fig_name = os.path.join(config["fig_folder"], mtype, fit_tpe)
                plot.plot_distribution_fit(
                    model_data[fit_tpe],
                    model_param[fit_tpe],
                    config["neurite_types"],
                    fig_name=fig_name,
                    ext=ext,
                )


def _build_all_models(morphologies, config, plot=False, ext="png"):
    """Building all the models in the list of models"""
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
    """ Run the model extraction from config file"""
    with open(config_file, "r") as filename:
        config = json.load(filename)

    L.info("Loading morphologies...")
    morphologies_dict = utils.create_morphologies_dict(
        config["morph_path"],
        mtypes_sort=config["mtypes_sort"],
        mtypes_file=config["mtypes_file"],
    )
    morphologies = io.load_morphologies_from_dict(
        config["morph_path"], morphologies_dict
    )

    L.info("Extracting model parameters...")
    models_params = _build_all_models(morphologies, config, plot=plot, ext=ext)

    with open(config["models_params_file"], "w") as json_file:
        json.dump(models_params, json_file, sort_keys=True, indent=4, cls=NumpyEncoder)


class DiameterWorker:
    """worker for building diameters"""

    def __init__(self, model, models_params, config, with_subfolders=False):
        self.model = model
        self.models_params = models_params[model]
        self.config = config[model]
        self.with_subfolders = with_subfolders

    def __call__(self, neuron_input):
        fname = neuron_input[0]
        mtype = neuron_input[1]
        neuron_name = os.path.splitext(fname)[0]

        try:
            file_format = ".h5"
            neuron = Morphology(
                os.path.join(self.config["morph_path"], neuron_name + file_format)
            )
        except (RawDataError, UnknownFileType):
            file_format = ".asc"
            neuron = Morphology(
                os.path.join(self.config["morph_path"], neuron_name + file_format)
            )

        build_diameters(
            neuron,
            self.model,
            self.models_params[mtype],
            self.config["neurite_types"],
            self.config,
        )

        if self.with_subfolders:
            save_path = os.path.join(self.config["new_morph_path"], mtype)
            if not os.path.exists(save_path):
                os.mkdir(save_path)
        else:
            save_path = self.config["new_morph_path"]

        neuron.write(os.path.join(save_path, neuron_name + file_format))


def run_diameters(config_file, models_params_file):
    """Build new diameters from config file and diameter model"""
    with open(config_file, "r") as filename:
        config = json.load(filename)

    with open(models_params_file, "r") as filename:
        models_params = json.load(filename)

    for model in config:
        L.info("Generating diameter with model %s", model)

        morphologies_dict = utils.create_morphologies_dict(
            config[model]["morph_path"],
            mtypes_sort=config[model]["mtypes_sort"],
            mtypes_file=config[model]["mtypes_file"],
        )

        # dirty hack to create subfolder as the original file structure
        if (
            len(
                os.path.dirname(morphologies_dict[list(morphologies_dict.keys())[0]][0])
            )
            > 1
        ):
            if not os.path.exists(config[model]["new_morph_path"]):
                os.mkdir(config[model]["new_morph_path"])
            with_subfolders = True
        else:
            with_subfolders = False

        worker = DiameterWorker(
            model, models_params, config, with_subfolders=with_subfolders
        )
        pool = multiprocessing.Pool(config[model]["n_cpu"])

        all_neurons = [
            [neuron, mtype]
            for mtype in morphologies_dict
            for neuron in morphologies_dict[mtype]
        ]
        list(tqdm(pool.imap(worker, all_neurons), total=len(all_neurons)))
