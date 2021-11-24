"""Utils functions."""

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
import os
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from neurom import COLS

from diameter_synthesis.exception import DiameterSynthesisError

L = logging.getLogger(__name__)


def _create_morphologies_dict_dat(morph_path, mtypes_file="neurondb.dat"):
    """Create dict to load the morphologies from a directory, with dat file.

    Args:
        morph_path (str): path to morphologies.
        mtype_file (str): path to dat file.

    Returns:
        dict: dictionary of morphologies keyed by mtypes.
    """
    # pylint wrongly determines `morph_name` as TextFileReader
    # pylint: disable=no-member
    morph_name = pd.read_csv(mtypes_file, sep=r"\s+", header=None)
    name_dict = defaultdict(list)
    if not morph_name.empty:
        first_name = morph_name.loc[0, 0]  # pylint: disable=no-member
        file_list = Path(morph_path).glob(first_name + "*")
        try:
            ext = next(file_list).suffix
        except StopIteration as e:
            raise DiameterSynthesisError(f"Could not find a file starting with {first_name}") from e
    for morph in morph_name.values:  # pylint: disable=no-member
        name_dict[morph[2]] += [Path(morph_path) / (morph[0] + ext)]
    return name_dict


def _create_morphologies_dict_folder(morph_path):
    """Create dict to load the morphologies from a directory, from folders.

    Args:
        morph_path (str): path to morphologies.

    Returns:
        dict: dictionary of morphologies keyed by mtypes.
    """
    name_dict = defaultdict(list)
    for mtype in Path(morph_path).iterdir():
        for fname in mtype.iterdir():
            if fname.suffix in [".h5", ".asc", ".swc"]:
                name_dict[mtype.parts[-1]] += [fname]
    return name_dict


def _create_morphologies_dict_all(morph_path):
    """Create dict to load the morphologies from a directory, all together.

    Args:
        morph_path (str): path to morphologies.

    Returns:
        dict: dictionary of morphologies with single key.
    """
    name_dict = {"generic_type": []}
    for fname in Path(morph_path).iterdir():
        if fname.suffix in [".h5", ".asc", ".swc"] and fname.exists():
            name_dict["generic_type"] += [fname]

    return name_dict


def create_morphologies_dict(morph_path, mtypes_file=None):
    """Create dict to load the morphologies from a directory, by mtype.

     Args:
        morph_path (str): path to morphologies.
        mtype_file (str): path to dat file.

    Returns:
        dict: dictionary of morphologies with single key.
    """
    if mtypes_file is None and next(Path(morph_path).iterdir()).is_dir():
        L.info("found folder structure per mtype")
        return _create_morphologies_dict_folder(morph_path)
    if mtypes_file is not None:
        # mtype_path = Path(morph_path) / mtypes_file
        mtype_path = Path(mtypes_file)
        if mtype_path.exists():
            if mtype_path.suffix == ".dat":
                L.info("found dat file")
                return _create_morphologies_dict_dat(morph_path, mtypes_file=mtypes_file)
            raise DiameterSynthesisError(
                f"neurondb file format {mtype_path.suffix} not implemented"
            )

    L.info("use all files as single mtype")
    return _create_morphologies_dict_all(morph_path)


def _set_diameters(section, diameters):
    """Set diameters.

    Args:
        section (neurom.Section): section to diametrize.
        diameters (list/ndarray): diameters.
    """
    section.morphio_section.diameters = diameters


def _get_mean_diameter(section):
    """Section mean diameter.

    It is obtained by averaging the segment truncated cone diameters and weighting them by their
    length. (neurom only)

    Args:
        section (neurom.Section): section to consider.

    Returns:
        float: mean diameter of the section.
    """
    segment_lengths = np.linalg.norm(
        section.points[1:, COLS.XYZ] - section.points[:-1, COLS.XYZ], axis=1
    )
    segment_mean_diams = section.points[1:, COLS.R] + section.points[:-1, COLS.R]
    return np.sum(segment_mean_diams * segment_lengths) / segment_lengths.sum()


def get_all_diameters(neuron):
    """Get all neuron diameters (morphio only).

    Args:
        neuron (morphio.mut.Morphology): neuron to consider

    Returns:
        list: all diameters
    """
    return [np.array(section.diameters) for section in neuron.iter()]


def set_all_diameters(neuron, all_diameters):
    """Set all neuron diameters (morphio only).

    Args:
        neuron (morphio.mut.Morophology): neuron to consider.
        all_diameters (list): list of section diameters.
    """
    for diameters, section in zip(all_diameters, neuron.iter()):
        section.diameters = diameters


def _get_diameters(section):
    """Get diameters (neurom only).

    Args:
        section (neurom.Section): section to consider.

    Return:
        list: diameters of section
    """
    return section.points[:, COLS.R] * 2.0


def redefine_diameter_section(section, diam_ind, diam_new):
    """Replace given diameters at indices diam_ind with values diam_new (morphio only).

    Args:
        section (neurom.Section): section to consider.
        diam_ind (list): indices of diameters, or points on the section.
        diam_new (list): corresponding diameters.
    """
    diameters = section.diameters
    diameters[diam_ind] = diam_new
    section.diameters = diameters


#######################################################################
# old parser functions (will be removed or integrated at some points) #
#######################################################################
def _create_morphologies_dict_json(
    morph_path,
    mtypes_file="neuronDB.xml",
    prefix="",
):
    """Create dict to load the morphologies from a directory, with json."""
    with open(mtypes_file, "r", encoding="utf-8") as filename:
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
    morph_path,
    mtypes_file="neuronDB.xml",
    ext=".asc",
    prefix="",
):
    """Create dict to load the morphologies from a directory, from xml."""
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
