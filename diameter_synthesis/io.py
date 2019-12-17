""" Helper functions for loading and writing morphologies """

import os
import logging

import neurom as nm
from neurom.exceptions import RawDataError
#from neurom.exceptions import UnknownFileType


L = logging.getLogger(__name__)


def _is_valid_spec(filename):
    return filename.endswith(('.h5', '.asc', '.swc'))


def iter_morphology_filenames(directory):
    """ Returns an iterator on morphology filenames """
    return filter(_is_valid_spec, os.listdir(directory))


def iter_morphology_filepaths(directory, filenames=None):
    """ Returns a generator of morphology filepaths by joining directory path and filenames """
    filenames = filenames or iter_morphology_filenames(directory)
    return (os.path.join(directory, filename) for filename in filenames)


def load_morphology(filepath):
    """ Returns a morphology object using NeuroM """
    print(type(nm.load_neuron(filepath)))
    return nm.load_neuron(filepath)


def load_neuron(name, model_name, directory):
    """ load the neuron morphology for later analysis """
    prefix = '{}_'.format(model_name) if model_name else ''
    try:
        filepath = os.path.join(directory, '{}{}.h5'.format(prefix, name))
    except:
        filepath = os.path.join(directory, '{}{}.asc'.format(prefix, name))
    return load_morphology(filepath)


def load_morphologies(filepaths):
    """ Returns a morphology generator. If the path is not valid the generator
        skips the path.
    """
    for filepath in filepaths:
        try:
            cell = load_morphology(filepath)
        except (RawDataError, UnknownFileType) as exc:
            L.warning('%s failed to load: %s', filepath, exc)
            continue
        yield cell


def load_morphologies_from_folder(directory, filenames, n_morphs_max=None):
    """ Load the morphologies from a list of files in a folder """

    filepaths_it = iter_morphology_filepaths(directory, filenames=filenames)
    morphs_it = load_morphologies(filepaths_it)

    if n_morphs_max is None:
        return list(morphs_it)

    return [cell for n, cell in enumerate(morphs_it) if n < n_morphs_max]


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
    return {mtype: load_morphologies_from_folder(directory, filenames) for
            mtype, filenames in filenames_per_mtype.items()}


def save_neuron(neuron, model, folder):
    """ save the neuron morphology for later analysis """

    if not os.path.isdir(folder):
        L.warning('Directory %s is created.', folder)
        os.mkdir(folder)

    filepath = os.path.join(folder, '{}_{}.asc'.format(model, neuron.name))
    neuron.write(filepath)
