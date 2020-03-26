""" utils functions """
import json
import logging
import os
from pathlib import Path
from collections import defaultdict
from xml.etree import ElementTree as ET

import numpy as np
from neurom import COLS
import pandas as pd

from diameter_synthesis.exception import DiameterSynthesisError

L = logging.getLogger(__name__)


def _create_morphologies_dict_dat(morph_path, mtypes_file="neurondb.dat"):
    """Create dict to load the morphologies from a directory, with dat file"""
    morph_name = pd.read_csv(mtypes_file, sep=" ")
    name_dict = defaultdict(list)
    ext = next(Path(morph_path).iterdir()).suffix
    for morph in morph_name.values:
        name_dict[morph[2]] += [Path(morph_path) / (morph[0] + ext)]
    return name_dict


def _create_morphologies_dict_folder(morph_path):
    """Create dict to load the morphologies from a directory, from folders"""
    name_dict = defaultdict(list)
    for mtype in Path(morph_path).iterdir():
        for fname in mtype.iterdir():
            if fname.suffix in [".h5", ".asc", ".swc"]:
                name_dict[mtype.parts[-1]] += [fname]
    return name_dict


def _create_morphologies_dict_all(morph_path):
    """Create dict to load the morphologies from a directory, all together"""
    name_dict = {"generic_type": []}
    for fname in Path(morph_path).iterdir():
        if fname.suffix in [".h5", ".asc", ".swc"] and fname.exists():
            name_dict["generic_type"] += [fname]

    return name_dict


def create_morphologies_dict(morph_path, mtypes_file=None):
    """ Create dict to load the morphologies from a directory, by mtypes or all at once """
    if mtypes_file is None and next(Path(morph_path).iterdir()).is_dir():
        L.info("found folder structure per mtype")
        return _create_morphologies_dict_folder(morph_path)
    if mtypes_file is not None:
        mtype_path = Path(morph_path) / mtypes_file
        if mtype_path.exists():
            if mtype_path.suffix == ".dat":
                L.info("found dat file")
                return _create_morphologies_dict_dat(
                    morph_path, mtypes_file=mtypes_file
                )
            raise DiameterSynthesisError(
                "neurondb file format {} not implemented".format(mtype_path.suffix)
            )

    L.info("use all files as single mtype")
    return _create_morphologies_dict_all(morph_path)


def _set_diameters(section, diameters):
    """set diameters (neurom v2 only)"""
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


#########################
# old parser functions
#########################
def _create_morphologies_dict_json(
    morph_path, mtypes_file="neuronDB.xml", prefix="",
):
    """ Create dict to load the morphologies from a directory, with json """

    with open(mtypes_file, "r") as filename:
        morph_name = json.load(filename)

    name_dict = {}
    for fname in os.listdir(morph_path):
        filepath = Path(morph_path) / fname
        if fname.endswith((".h5", ".asc", ".swc")) and filepath.exists():
            mtype = morph_name[Path(fname).stem]
            if mtype in name_dict:
                name_dict[mtype] += [prefix + fname]
            else:
                name_dict[mtype] = [prefix + fname]

    return name_dict


def _create_morphologies_dict_xml(
    morph_path, mtypes_file="neuronDB.xml", ext=".asc", prefix="",
):
    """ Create dict to load the morphologies from a directory, from xml """
    # first load the neuronDB.xml file

    filedb = ET.parse(morph_path + mtypes_file)
    root = filedb.findall("listing")[0]
    morphs = root.findall("morphology")

    name_dict = {}
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
