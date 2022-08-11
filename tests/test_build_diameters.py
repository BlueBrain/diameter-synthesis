"""Test the build_diameters module."""

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

import logging

import morphio.mut
import numpy as np
import pytest
from morphio import PointLevel
from morphio import SectionType

from diameter_synthesis import build_diameters
from diameter_synthesis.exception import DiameterSynthesisError

from .testing_tools import _compare_diameters

NEURITE_TYPES = ["basal_dendrite", "apical_dendrite"]


def test_build_simple(config, model_params, small_morph, test_data_path):
    """Test the main build function to diametrize a neuron."""
    build_diameters.build(small_morph, NEURITE_TYPES, model_params, config)
    expected = morphio.mut.Morphology(test_data_path / "small_morph_diametrized.asc")
    _compare_diameters(expected, small_morph)


def test_build_simple_with_apical_sections(config, model_params, small_morph, test_data_path):
    """Test the main build function to diametrize a neuron with a known apical section."""
    model_params["apical_point_sec_ids"] = [6]

    build_diameters.build(small_morph, NEURITE_TYPES, model_params, config)
    expected = morphio.mut.Morphology(test_data_path / "small_morph_apical_diametrized.asc")
    _compare_diameters(expected, small_morph)


def test_build_simple_with_several_apical_sections(
    config, model_params, small_morph, test_data_path
):
    """Test the main build function to diametrize a neuron with a known apical section."""
    # Duplicate the apical
    for section in small_morph.sections[4].iter():
        if section.id != 4:
            append_fun = small_morph.sections[section.parent.id + 13].append_section
        else:
            append_fun = small_morph.append_root_section
        pts = section.points.copy()
        pts[:, 0] += 1000
        append_fun(PointLevel(pts, section.diameters.copy()), SectionType(3))

    model_params["apical_point_sec_ids"] = [6, 19]

    build_diameters.build(small_morph, NEURITE_TYPES, model_params, config)

    expected = morphio.mut.Morphology(test_data_path / "small_morph_several_apical_diametrized.asc")
    _compare_diameters(expected, small_morph)


def test_build(config, model_params, neuron, neuron_diametrized):
    """Test the main build function to diametrize a neuron."""
    build_diameters.build(neuron, NEURITE_TYPES, model_params, config)
    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_seed(config, model_params, neuron, neuron_diametrized):
    """Test main build without seed in config."""
    rng = np.random.default_rng(config.pop("seed"))

    build_diameters.build(neuron, NEURITE_TYPES, model_params, config, rng)
    _compare_diameters(neuron_diametrized, neuron)


def test_build_multiple_models_warning(config, model_params, neuron, neuron_diametrized, caplog):
    """Test main build with multiple models in config."""
    config["models"].append("SKIPPED MODEL")

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, NEURITE_TYPES, model_params, config)

    module, level, entry = caplog.record_tuples[0]
    assert module == "diameter_synthesis.build_diameters"
    assert level == 30
    assert entry == "Several models provided, we will only use the first"

    _compare_diameters(neuron_diametrized, neuron)


def test_build_one_sample(test_data_path, config, model_params, neuron):
    """Test main build with 1 sample in config."""
    config["n_samples"] = 1

    build_diameters.build(neuron, NEURITE_TYPES, model_params, config)

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_1_sample.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_sequential(test_data_path, config, model_params, neuron):
    """Test main build with no sequential in config."""
    config["models"] = ["astrocyte"]
    config["n_samples"] = 1
    model_params["sibling_ratios"]["apical_dendrite"]["sequential"] = None

    build_diameters.build(neuron, NEURITE_TYPES, model_params, config)

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_astrocyte.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_small_n_tries_warning(config, model_params, neuron, caplog):
    """Test main build with trunk_max_tries = 1."""
    config["trunk_max_tries"] = 1

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, NEURITE_TYPES, model_params, config)
    assert len(caplog.record_tuples) == 5
    neurite_types = [
        "apical_dendrite",
        "apical_dendrite",
        "apical_dendrite",
        "apical_dendrite",
        "apical_dendrite",
    ]
    for i in zip(neurite_types, caplog.record_tuples):
        neurite_type, (module, level, entry) = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == f"max tries attained with {neurite_type}"


def test_build_small_trunk_diam_warning(config, model_params, neuron, caplog):
    """Test trunk diam warning in main build."""
    TRUNK_FRAC_DECREASE = build_diameters.TRUNK_FRAC_DECREASE
    try:
        build_diameters.TRUNK_FRAC_DECREASE = 999
        caplog.clear()
        caplog.set_level(logging.WARNING)
        build_diameters.build(neuron, NEURITE_TYPES, model_params, config)
    finally:
        build_diameters.TRUNK_FRAC_DECREASE = TRUNK_FRAC_DECREASE

    assert len(caplog.record_tuples) == 1
    for i in caplog.record_tuples:
        module, level, entry = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == "sampled trunk diameter < 0.01, so use 1 instead"


def test_select_model():
    """Test the _select_model function."""
    # pylint: disable=protected-access
    # pylint: disable=no-member
    f_generic = build_diameters._select_model("generic")
    assert f_generic.args == (
        {
            "mode_sibling": "threshold",
            "mode_diameter_power_relation": "threshold",
            "with_asymmetry": True,
            "reduction_factor_max": 1.0,
        },
    )

    f_astrocyte = build_diameters._select_model("astrocyte")
    assert f_astrocyte.args == (
        {
            "mode_sibling": "generic",
            "mode_diameter_power_relation": "generic",
            "with_asymmetry": True,
            "reduction_factor_max": 3.0,
        },
    )

    with pytest.raises(DiameterSynthesisError):
        build_diameters._select_model("UNKNOWN")


def test_sample_sibling_ratio(model_params):
    """Test the _sample_sibling_ratio function."""
    # pylint: disable=protected-access
    d_generic = build_diameters._sample_sibling_ratio(
        model_params, "basal_dendrite", False, mode="generic"
    )
    assert d_generic == pytest.approx(0.6925189949520985)

    d_threshold = build_diameters._sample_sibling_ratio(
        model_params, "basal_dendrite", False, mode="threshold"
    )
    assert d_threshold == pytest.approx(0.5222098881400881)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_sibling_ratio(model_params, "basal_dendrite", False, mode="UNKNOWN")


def test_sample_diameter_power_relation(model_params):
    """Test the _sample_diameter_power_relation function."""
    # pylint: disable=protected-access
    d_generic = build_diameters._sample_diameter_power_relation(
        model_params, "basal_dendrite", False, mode="generic"
    )
    assert d_generic == pytest.approx(5.281970058326886)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params, "basal_dendrite", False, mode="threshold"
    )
    assert d_threshold == pytest.approx(5.915270977501953)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params, "basal_dendrite", False, mode="exact"
    )
    assert d_threshold == pytest.approx(1)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_diameter_power_relation(
            model_params, "basal_dendrite", False, mode="UNKNOWN"
        )


def test_sample_daughter_diameters(neuron, model_params):
    """Test the _sample_diameter_power_relation function."""
    # pylint: disable=protected-access
    param_tree = {
        "mode_sibling": "threshold",
        "mode_diameter_power_relation": "threshold",
        "with_asymmetry": True,
        "reduction_factor_max": 1.0,
        "neurite_type": "basal_dendrite",
        "asymmetry_threshold": 1.0,
        "trunk_diam": 2.131444808766456,
        "tot_length": 1085.1178037524223,
        "terminal_diam": 0.35602967635703164,
        "major_sections": [],
    }
    d = build_diameters._sample_daughter_diameters(neuron.sections[0], model_params, param_tree)
    assert d == pytest.approx([0.9856948, 0.68261237])

    param_tree["with_asymmetry"] = False
    d_no_asymmetry = build_diameters._sample_daughter_diameters(
        neuron.sections[0], model_params, param_tree
    )
    assert d_no_asymmetry == pytest.approx([1.23559368, 0.77282112])


def test_diametrize_axon(small_morph):
    """Test the axon diametrizer."""
    build_diameters.diametrize_axon(small_morph)

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [0.2, 0.2, 0.2],  # basal
        1: [0.2, 0.2, 0.2],  # basal
        2: [0.2, 0.2, 0.2],  # basal
        3: [0.2, 0.2, 0.2],  # basal
        4: [0.2, 0.2, 0.2],  # apical
        5: [0.2, 0.2, 0.2],  # apical
        6: [0.2, 0.2, 0.2],  # apical
        7: [0.2, 0.2, 0.2],  # apical
        8: [0.2, 0.2, 0.2],  # apical
        9: [0.2, 0.2, 0.2],  # apical
        10: [0.2, 0.2, 0.2],  # apical
        11: [0.2, 0.2, 0.2],  # apical
        12: [0.2, 0.2, 0.2],  # apical
        13: [0.2, 0.2, 0.9],  # axon
        14: [0.2, 0.2, 0.9],  # axon
        15: [0.376628, 0.326628, 0.276628],  # basal
        16: [0.753256, 0.703256, 0.653256],  # basal
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph, rtol=1e-6)
