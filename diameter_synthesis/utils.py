""" utils functions """
import json
import logging
import os
from xml.etree import ElementTree as ET

import numpy as np

from neurom import COLS
from neurom.core import iter_sections

from .io import load_morphologies_from_dict

L = logging.getLogger(__name__)

ROUND = 4  # number of digits for the fitted parameters
MIN_DATA_POINTS = 1  # minimum number of points to fit a distribution
# maximum value for the a (shape) parameter of fits (can get really large when low number of points)
A_MAX = 4
A_MIN = 0.3
SPLINE_SMOOTH = 0.0001

# taken from jira ticket BBPP82-127
FORBIDDEN_MTYPES = [
    "L4_NGC",
    "L4_CHC",
    "L5_NGC",
    "L6_BP",
    "L6_CHC",
    "L6_DBC",
]


def create_morphologies_dict_json(
        morph_path,
        mtype_file="neuronDB.xml",
        prefix="",
):
    """ Create dict to load the morphologies from a directory, with json """

    with open(mtype_file, "r") as filename:
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


def create_morphologies_dict_folder(
        morph_path,
        prefix="",
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


def create_morphologies_dict_xml(
        morph_path,
        mtypes_sort="all",
        mtype_file="neuronDB.xml",
        ext=".asc",
        prefix="",
):
    """ Create dict to load the morphologies from a directory, from xml """
    # first load the neuronDB.xml file

    filedb = ET.parse(morph_path + mtype_file)
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

                # hack to not consider some bad mtypes
                if mtype not in FORBIDDEN_MTYPES:
                    # if it is a new mtype, add en entry to name_dict
                    if mtype not in name_dict:
                        name_dict[mtype] = [prefix + morph.find("name").text + ext]
                    elif mtype in name_dict:
                        name_dict[mtype] += [prefix + morph.find("name").text + ext]

            except Exception as exc:  # pylint: disable=broad-except
                L.exception("Failed to process %s", exc)

    return name_dict


def create_morphologies_dict_all(
        morph_path,
        prefix="",
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
        morph_path,
        mtypes_sort="all",
        mtype_file="neuronDB.xml",
        ext=".asc",
        prefix="",
):
    """ Create dict to load the morphologies from a directory, by mtypes or all at once """

    try:

        name_dict = create_morphologies_dict_json(
            morph_path,
            mtype_file=mtype_file,
            prefix=prefix)

        L.info("found morph.json file")
        return name_dict

    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        name_dict = create_morphologies_dict_folder(
            morph_path,
            prefix=prefix)

        L.info("found folder structure per mtype")
        return name_dict

    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        name_dict = create_morphologies_dict_xml(
            morph_path,
            mtypes_sort=mtypes_sort,
            mtype_file=mtype_file,
            ext=ext,
            prefix=prefix)

        L.info("found neuronDB.xml")
        return name_dict

    except BaseException:  # pylint: disable=broad-except
        pass

    try:
        L.info("use all files as single mtype")

        name_dict = create_morphologies_dict_all(
            morph_path,
            prefix=prefix)

        return name_dict

    except Exception as exc:  # pylint: disable=broad-except
        L.info("Could not load any files with exception: %s", exc)


def load_morphologies(
        morph_path,
        mtypes_sort="all",
        mtype_file="./neuronDB.xml",
        ext=".asc",
        prefix="",
):
    """ Load the morphologies from a directory, by mtypes or all at once """
    name_dict = create_morphologies_dict(
        morph_path,
        mtypes_sort=mtypes_sort,
        mtype_file=mtype_file,
        ext=ext,
        prefix=prefix,
    )

    return load_morphologies_from_dict(morph_path, name_dict)


###############################
# diameter handling functions #
###############################


def set_diameters(section, diameters):
    """hack to set diameters with neurom"""
    section.morphio_section.diameters = diameters


def get_mean_diameter(section):
    """ Section mean diameter by averaging the segment truncated cone
    diameters and weighting them by their length.
    """
    points = section.morphio_section.points[:, COLS.XYZ]
    radii = section.morphio_section.points[:, COLS.R]

    segment_lengths = np.linalg.norm(points[1:] - points[:-1], axis=1)

    segment_mean_diams = radii[1:] + radii[:-1]

    return np.sum(segment_mean_diams * segment_lengths) / segment_lengths.sum()


def get_all_diameters(neuron):
    """get all neuron diameters"""
    return list(map(get_diameters, iter_sections(neuron)))


def set_all_diameters(neuron, diameters):
    """get all neuron diameters"""

    i = 0
    for neurite in neuron.neurites:
        for section in iter_sections(neurite):
            set_diameters(section, diameters[i])
            i += 1


def get_diameters(section):
    """hack to get diameters with neurom (faster to access morphio directly)"""
    return section.morphio_section.diameters


def redefine_diameter_section(section, diam_ind, diam_new):
    """Hack to replace one diameter at index diam_ind with value diam_new"""

    diameters = get_diameters(section)
    diameters[diam_ind] = diam_new
    set_diameters(section, diameters)


########################
# additional functions #
########################


def section_lengths(section):
    """Computes all segment lengths within section"""

    vecs = np.diff(section.points, axis=0)[:, COLS.XYZ]
    len_square = [np.dot(p, p) for p in vecs]
    return list(np.cumsum(np.sqrt(len_square)))


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
        while (values[-1] < n_min and max_data - min_data > 0.1 * diff_max):
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
