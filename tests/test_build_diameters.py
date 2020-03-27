from pathlib import Path
import os
import json
import pytest

import morphio.mut
from neurom.core._neuron import Neuron
import neurom as nm

from diameter_synthesis import build_diameters
from diameter_synthesis.exception import DiameterSynthesisError

_path = Path(__file__).absolute().parent


def _copy_diameters(neuron_a, neuron_b):
    for section_a, section_b in zip(
        nm.iter_sections(neuron_a), nm.iter_sections(neuron_b)
    ):
        section_a.diameters = section_b.diameters


def test_build():
    """test the main build function to diametrize a neuron"""

    neurite_types = ["basal", "apical"]

    mtype = "L5_TPC:A"
    model = "generic"

    with open(_path / "data" / "config.json", "r") as filename:
        config = json.load(filename)

    with open(_path / "data" / "model_params.json", "r") as filename:
        model_params = json.load(filename)

    neuron = morphio.mut.Morphology(_path / "data" / "C030796A-P3.h5")
    model = "generic"
    build_diameters.build(
        neuron, model_params[model][mtype], neurite_types, config[model]
    )

    # neuron.write(_path / "data" / "C030796A-P3_diametrized.h5")

    neuron_diametrized = morphio.Morphology(
        _path / "data" / "C030796A-P3_diametrized.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def _compare_diameters(neuron_a, neuron_b):
    """assert of all diameters are the same"""
    for section_a, section_b in zip(neuron_a.iter(), neuron_b.iter()):
        assert all(section_a.diameters == section_b.diameters)


if __name__ == "__main__":
    test_build()
