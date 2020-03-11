""" Helper functions for loading and writing morphologies """

import logging
import os

import neurom as nm
from morphio import RawDataError, UnknownFileType

L = logging.getLogger(__name__)


def _is_valid_spec(filename):
    return filename.endswith((".h5", ".asc", ".swc"))


def iter_morphology_filenames(directory):
    """ Returns an iterator on morphology filenames """
    return filter(_is_valid_spec, os.listdir(directory))


def iter_morphology_filepaths(directory, filenames=None):
    """ Returns a generator of morphology filepaths by joining directory path and filenames """
    filenames = filenames or iter_morphology_filenames(directory)
    return (os.path.join(directory, filename) for filename in filenames)


def load_neuron(name, directory, model_name=""):
    """ load the neuron morphology for later analysis """
    prefix = "{}_".format(model_name) if model_name else ""
    try:
        try:
            filepath = os.path.join(directory, "{}{}.h5".format(prefix, name))
            return nm.load_neuron(filepath)
        except (RawDataError, UnknownFileType):
            filepath = os.path.join(directory, "{}{}.asc".format(prefix, name))
            return nm.load_neuron(filepath)
    except (RawDataError, UnknownFileType):
        L.exception("file not found")


def load_morphologies(filepaths):
    """ Returns a morphology generator. If the path is not valid the generator
        skips the path.
    """
    for filepath in filepaths:
        try:
            cell = load_neuron(
                os.path.splitext(os.path.basename(filepath))[0],
                os.path.dirname(filepath),
            )
        except (RawDataError, UnknownFileType):
            continue
        yield cell


def load_morphologies_from_folder(directory, filenames):
    """ Load the morphologies from a list of files in a folder """
    filepaths_it = iter_morphology_filepaths(directory, filenames=filenames)
    return list(load_morphologies(filepaths_it))


def load_morphologies_from_dict(directory, filenames_per_mtype):
    """ Load the morphologies from a list of files

    Args:
        morphology_directory: string
            The directory containing all morphologies
        filenames_per_mtype: dict
            Dictionary with mtype names as keys and cell filenames as values
    Returns:
        Dictionary with mtypes as keys and cells as values
    """
    return {
        mtype: load_morphologies_from_folder(directory, filenames)
        for mtype, filenames in filenames_per_mtype.items()
    }


def save_neuron(neuron, model, folder):
    """ save the neuron morphology for later analysis """

    if not os.path.exists(folder):
        L.warning("Directory %s is created.", folder)
        os.mkdir(folder)

    if model == "generic":
        filepath = os.path.join(folder, "{}.asc".format(neuron.name))
    else:
        filepath = os.path.join(folder, "{}_{}.asc".format(model, neuron.name))
    neuron.write(filepath)
