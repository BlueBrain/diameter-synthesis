""" utils functions """
import json
import logging
import os
from xml.etree import ElementTree as ET

import numpy as np
from neurom import COLS
import pandas as pd

L = logging.getLogger(__name__)


def _create_morphologies_dict_dat(  # pylint: disable=unused-argument
    morph_path, mtypes_file="neurondb.dat", prefix="", ext=".asc"
):
    """ Create dict to load the morphologies from a directory, with dat file """
    morph_name = pd.read_csv(mtypes_file, sep=" ")
    name_dict = {}
    for morph in morph_name.values:
        if morph[2] in name_dict:
            name_dict[morph[2]] += [prefix + morph[0] + ext]
        else:
            name_dict[morph[2]] = [prefix + morph[0] + ext]
    return name_dict


def _create_morphologies_dict_json(
    morph_path, mtypes_file="neuronDB.xml", prefix="",
):
    """ Create dict to load the morphologies from a directory, with json """

    with open(mtypes_file, "r") as filename:
        morph_name = json.load(filename)

    name_dict = {}
    for fname in os.listdir(morph_path):
        filepath = os.path.join(morph_path, fname)
        if fname.endswith((".h5", ".asc", ".swc")) and os.path.exists(filepath):
            # get the mtype
            mtype = morph_name[os.path.splitext(fname)[0]]
            if mtype in name_dict:
                name_dict[mtype] += [prefix + fname]
            else:
                name_dict[mtype] = [prefix + fname]

    return name_dict


def _create_morphologies_dict_folder(
    morph_path, prefix="",
):
    """ Create dict to load the morphologies from a directory, from folders """

    n_morphs = 0
    name_dict = {}
    for fold_name in os.listdir(morph_path):
        for fname in os.listdir(os.path.join(morph_path, fold_name)):
            if fname.endswith((".h5", ".asc", ".swc")):
                # get the mtype
                mtype = fold_name
                # if it is an old mtype, convert it
                if mtype in name_dict:
                    name_dict[mtype] += [os.path.join(fold_name, prefix + fname)]
                else:
                    name_dict[mtype] = [os.path.join(fold_name, prefix + fname)]
                n_morphs += 1
    return name_dict


def _create_morphologies_dict_xml(
    morph_path, mtypes_sort="all", mtypes_file="neuronDB.xml", ext=".asc", prefix="",
):
    """ Create dict to load the morphologies from a directory, from xml """
    # first load the neuronDB.xml file

    filedb = ET.parse(morph_path + mtypes_file)
    root = filedb.findall("listing")[0]
    morphs = root.findall("morphology")

    name_dict = {}
    if mtypes_sort == "all":
        name_dict["all_types"] = []

    if len(morphs) > 0:
        for morph in morphs:
            try:
                # Define mtypes
                mtype = morph.find("mtype").text
                # Define subtypes (if they exist)
                if morph.find("msubtype").text:
                    mtype = mtype + ":" + morph.find("msubtype").text

                if mtype not in name_dict:
                    name_dict[mtype] = [prefix + morph.find("name").text + ext]
                elif mtype in name_dict:
                    name_dict[mtype] += [prefix + morph.find("name").text + ext]

            except Exception as exc:  # pylint: disable=broad-except
                L.exception("Failed to process %s", exc)

    return name_dict


def _create_morphologies_dict_all(
    morph_path, prefix="",
):
    """ Create dict to load the morphologies from a directory, all together """

    n_morphs = 0
    name_dict = {}
    name_dict["generic_type"] = []
    for fname in os.listdir(morph_path):
        filepath = os.path.join(morph_path, fname)
        if fname.endswith((".h5", ".asc", ".swc")) and os.path.exists(filepath):
            name_dict["generic_type"] += [prefix + fname]
            n_morphs += 1

    return name_dict


def create_morphologies_dict(
    morph_path, mtypes_sort="all", mtypes_file="neuronDB.xml", ext=".asc", prefix="",
):
    """ Create dict to load the morphologies from a directory, by mtypes or all at once """
    try:
        name_dict = _create_morphologies_dict_dat(
            morph_path, mtypes_file=mtypes_file, prefix=prefix
        )
        L.info("found dat file")
        return name_dict
    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        name_dict = _create_morphologies_dict_json(
            morph_path, mtypes_file=mtypes_file, prefix=prefix
        )
        L.info("found morph.json file")
        return name_dict
    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        name_dict = _create_morphologies_dict_folder(morph_path, prefix=prefix)
        L.info("found folder structure per mtype")
        return name_dict
    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        name_dict = _create_morphologies_dict_xml(
            morph_path,
            mtypes_sort=mtypes_sort,
            mtypes_file=mtypes_file,
            ext=ext,
            prefix=prefix,
        )
        L.info("found neuronDB.xml")
        return name_dict
    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        L.info("use all files as single mtype")
        name_dict = _create_morphologies_dict_all(morph_path, prefix=prefix)
        return name_dict
    except Exception as exc:  # pylint: disable=broad-except
        L.info("Could not load any files with exception: %s", exc)


def _set_diameters(section, diameters):
    """set diameters (neurom)"""
    section.morphio_section.diameters = diameters


def _get_mean_diameter(section):
    """Section mean diameter by averaging the segment truncated cone
    diameters and weighting them by their length. (morphio)"""
    segment_lengths = np.linalg.norm(section.points[1:] - section.points[:-1], axis=1)
    segment_mean_diams = (section.diameters[1:] + section.diameters[:-1]) / 2.0
    return np.sum(segment_mean_diams * segment_lengths) / segment_lengths.sum()


def _get_all_diameters(neuron):
    """get all neuron diameters (morphio)"""
    return [section.diameters for section in neuron.iter()]


def _set_all_diameters(neuron, diameters):
    """set all neuron diameters (morphio)"""
    for diameter, section in zip(diameters, neuron.iter()):
        section.diameters = diameter


def _get_diameters(section):
    """get diameters (neurom)"""
    return section.points[:, COLS.R] * 2.0


def _redefine_diameter_section(section, diam_ind, diam_new):
    """Hack to replace one diameter at index diam_ind with value diam_new (morphio)"""
    diameters = section.diameters
    diameters[diam_ind] = diam_new
    section.diameters = diameters


def tqdm_disable(morphologies):
    """ to have a single progression bar """

    if len(morphologies) > 1:
        tqdm_1 = False
        tqdm_2 = True
    else:
        tqdm_1 = True
        tqdm_2 = False

    return tqdm_1, tqdm_2


def set_bins(data, n_bins, n_min=20):
    """ find a good set of bins to avoid undersampling """

    # try to bin uniformly
    max_data = np.max(data)
    min_data = np.min(data)
    diff_max = max_data - min_data

    values, bins = np.histogram(data, bins=n_bins, range=(min_data, max_data))

    def reduce_bounds(values, bins, data):
        max_data = np.max(data)
        min_data = np.min(data)

        # if the last bins have to few points, reduce the window
        # second condition is to prevent shrinking for ever
        while values[-1] < n_min and max_data - min_data > 0.1 * diff_max:
            max_data = bins[-2]
            values, bins = np.histogram(data, bins=n_bins, range=(min_data, max_data))

        # if the first bins have to few points, reduce the window
        while values[0] < n_min and max_data - min_data > 0.1 * diff_max:
            min_data = bins[1]
            values, bins = np.histogram(data, bins=n_bins, range=(min_data, max_data))

        return values, bins

    # find new bounds
    values, bins = reduce_bounds(values, bins, data)

    # if bins have to few elements, reduce the number of bins and readjust bounds
    while len(values[values < n_min]) > 0 and n_bins > 4:
        n_bins -= 1
        max_data = np.max(data)
        min_data = np.min(data)
        values, bins = np.histogram(data, bins=n_bins, range=(min_data, max_data))

        values, bins = reduce_bounds(values, bins, data)
    return bins, values
