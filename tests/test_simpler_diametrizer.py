"""Test the build_models module."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import morphio.mut
import neurom as nm
from morphio import SectionType

from diameter_synthesis import build_diameters
from diameter_synthesis import build_models

from .testing_tools import _compare_diameters
from .testing_tools import compare_dicts


def test_build_model(single_pop, simpler_model_params, simpler_model_data):
    """Test the build function."""
    # Test with generic model
    config = {
        "models": ["simpler"],
        "neurite_types": ["basal_dendrite", "apical_dendrite"],
    }
    res = build_models.build(single_pop, config, with_data=True)
    res_models_params = build_models.build(single_pop, config, with_data=False)

    assert len(res) == 2
    assert compare_dicts(res, [simpler_model_params, simpler_model_data], precision=3)
    assert compare_dicts(res_models_params, simpler_model_params, precision=3)


def test_build_model_missing_neurite_type(tmpdir, neuron, simpler_model_params, simpler_model_data):
    """Test the build function with missing neurites of given types."""
    # Remove axon from neuron
    for sec in neuron.root_sections:
        if sec.type == SectionType.axon:
            neuron.delete_section(sec)
    neuron_path = tmpdir / "neuron.swc"
    neuron.write(neuron_path)
    single_pop = [nm.load_morphology(neuron_path)]

    # Test with generic model
    config = {
        "models": ["simpler"],
        "neurite_types": ["basal_dendrite", "axon"],
    }
    res = build_models.build(single_pop, config, with_data=True)
    res_models_params = build_models.build(single_pop, config, with_data=False)

    # Build expected results
    simpler_model_params["axon"] = []
    del simpler_model_params["apical_dendrite"]
    for i in simpler_model_data:
        i["axon"] = []
        del i["apical_dendrite"]
    simpler_model_data[2]["axon"] = None

    # Check results
    assert len(res) == 2
    assert compare_dicts(res, [simpler_model_params, simpler_model_data], precision=3)
    assert compare_dicts(res_models_params, simpler_model_params, precision=3)


def test_build_diameters(simpler_config, simpler_model_params, small_morph, test_data_path):
    """Test the main build function to diametrize a neuron."""
    neurite_types = ["basal_dendrite", "apical_dendrite"]

    build_diameters.build(small_morph, neurite_types, simpler_model_params, simpler_config)
    expected = morphio.mut.Morphology(test_data_path / "simpler_morph_diametrized.asc")
    _compare_diameters(expected, small_morph)
