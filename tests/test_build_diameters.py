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


def test_build_simple(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    build_diameters.build(small_morph, neurite_types, model_params[model][mtype], config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.2997824, 0.9000233, 0.72673666],  # basal
        1: [1.6028887, 1.2359781, 0.91419667],  # basal
        2: [0.7737082, 0.6514245, 0.61029184],  # basal
        3: [0.727182, 0.5553919, 0.540321],  # basal
        4: [3.2547958, 2.7254605, 2.1961255],  # apical
        5: [1.2261434, 0.5947401, 0.44064492],  # apical
        6: [1.2281191, 1.0033262, 0.85901225],  # apical
        7: [0.5838221, 0.40629458, 0.38336557],  # apical
        8: [0.6178318, 0.6467553, 0.70707476],  # apical
        9: [0.58992374, 0.4161563, 0.4003501],  # apical
        10: [0.5538684, 0.47851545, 0.4900443],  # apical
        11: [0.44196147, 0.42748648, 0.42748648],  # apical
        12: [0.43044534, 0.38146895, 0.39314994],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build_simple_with_apical_sections(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron with a known apical section."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    model_params[model][mtype]["apical_point_sec_ids"] = [6]

    build_diameters.build(small_morph, neurite_types, model_params[model][mtype], config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.6146091, 1.1505961, 0.74188745],  # basal
        1: [1.8352976, 1.5289228, 1.2671413],  # basal
        2: [0.96662533, 0.72055197, 0.6243407],  # basal
        3: [0.894905, 0.653613, 0.56571424],  # basal
        4: [3.5404124, 2.9960587, 2.4517055],  # apical
        5: [0.4869547, 0.45853767, 0.4576969],  # apical
        6: [2.4517055, 2.073075, 1.769174],  # apical
        7: [0.8419878, 0.60131085, 0.44629398],  # apical
        8: [1.3291866, 1.1325629, 0.979848],  # apical
        9: [0.59224963, 0.48624295, 0.4726326],  # apical
        10: [0.6087809, 0.511539, 0.4896373],  # apical
        11: [0.41421217, 0.3673572, 0.39478624],  # apical
        12: [0.43502575, 0.43636972, 0.4498969],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build_simple_with_several_apical_sections(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron with a known apical section."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    # Duplicate the apical
    for section in small_morph.sections[4].iter():
        if section.id != 4:
            append_fun = small_morph.sections[section.parent.id + 13].append_section
        else:
            append_fun = small_morph.append_root_section
        pts = section.points.copy()
        pts[:, 0] += 1000
        append_fun(PointLevel(pts, section.diameters.copy()), SectionType(3))

    model_params[model][mtype]["apical_point_sec_ids"] = [6, 19]

    build_diameters.build(small_morph, neurite_types, model_params[model][mtype], config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.6581707, 1.0526989, 0.6265956],  # basal
        1: [1.8037863, 1.4524208, 1.1289493],  # basal
        2: [0.7909685, 0.64985496, 0.56315005],  # basal
        3: [0.86445045, 0.6507469, 0.56271935],  # basal
        4: [3.3207638, 2.783594, 2.294971],  # apical
        5: [0.54798895, 0.4854589, 0.48255223],  # apical
        6: [2.294971, 1.9488478, 1.6358488],  # apical
        7: [0.47745314, 0.44113955, 0.48057836],  # apical
        8: [1.6358488, 1.5305059, 1.4305559],  # apical
        9: [0.86752427, 0.6409085, 0.47136444],  # apical
        10: [0.9211162, 0.71516794, 0.56216586],  # apical
        11: [0.43250427, 0.41028976, 0.4410184],  # apical
        12: [0.45779485, 0.42024216, 0.39077753],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
        17: [1.2681566, 0.8854807, 0.66164625],  # apical
        18: [0.5188144, 0.45719615, 0.45719615],  # apical
        19: [0.66164625, 0.59810007, 0.6117419],  # apical
        20: [0.6117419, 0.46071988, 0.46071988],  # apical
        21: [0.4938531, 0.46393174, 0.49159345],  # apical
        22: [0.42158556, 0.41690236, 0.41896287],  # apical
        23: [0.41675425, 0.4613071, 0.49716002],  # apical
        24: [0.41612378, 0.4461469, 0.46464235],  # apical
        25: [0.4088277, 0.4126462, 0.42915756],  # apical
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build(config, model_params, neuron, neuron_diametrized):
    """Test the main build function to diametrize a neuron."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])

    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_seed(config, model_params, neuron, neuron_diametrized):
    """Test main build without seed in config."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    seed = config[model].pop("seed")
    rng = np.random.default_rng(seed)

    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model], rng)

    _compare_diameters(neuron_diametrized, neuron)


def test_build_multiple_models_warning(config, model_params, neuron, neuron_diametrized, caplog):
    """Test main build with multiple models in config."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["models"].append("SKIPPED MODEL")

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])

    module, level, entry = caplog.record_tuples[0]
    assert module == "diameter_synthesis.build_diameters"
    assert level == 30
    assert entry == "Several models provided, we will only use the first"

    _compare_diameters(neuron_diametrized, neuron)


def test_build_one_sample(test_data_path, config, model_params, neuron):
    """Test main build with 1 sample in config."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["n_samples"] = 1

    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_1_sample.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_sequential(test_data_path, config, model_params, neuron):
    """Test main build with no sequential in config."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "astrocyte"

    config[model] = config["generic"]
    config[model]["models"] = [model]
    config[model]["n_samples"] = 1
    model_params[model] = model_params["generic"]
    model_params[model][mtype]["sibling_ratios"]["apical"]["sequential"] = None

    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_astrocyte.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_small_n_tries_warning(config, model_params, neuron, caplog):
    """Test main build with trunk_max_tries = 1."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["trunk_max_tries"] = 1

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])

    assert len(caplog.record_tuples) == 8
    neurite_types = [
        "basal",
        "apical",
        "apical",
        "basal",
        "apical",
        "apical",
        "apical",
        "apical",
    ]
    for i in zip(neurite_types, caplog.record_tuples):
        neurite_type, (module, level, entry) = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == f"max tries attained with {neurite_type}"


def test_build_small_trunk_diam_warning(config, model_params, neuron, caplog):
    """Test trunk diam warning in main build."""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    TRUNK_FRAC_DECREASE = build_diameters.TRUNK_FRAC_DECREASE
    try:
        build_diameters.TRUNK_FRAC_DECREASE = 999
        caplog.clear()
        caplog.set_level(logging.WARNING)
        build_diameters.build(neuron, neurite_types, model_params[model][mtype], config[model])
    finally:
        build_diameters.TRUNK_FRAC_DECREASE = TRUNK_FRAC_DECREASE

    assert len(caplog.record_tuples) == 8
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
    mtype = "L5_TPC:A"
    model = "generic"
    d_generic = build_diameters._sample_sibling_ratio(
        model_params[model][mtype], "basal", False, mode="generic"
    )
    assert d_generic == pytest.approx(0.8253851329090136)

    d_threshold = build_diameters._sample_sibling_ratio(
        model_params[model][mtype], "basal", False, mode="threshold"
    )
    assert d_threshold == pytest.approx(0.7244487906052952)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_sibling_ratio(
            model_params[model][mtype], "basal", False, mode="UNKNOWN"
        )


def test_sample_diameter_power_relation(model_params):
    """Test the _sample_diameter_power_relation function."""
    # pylint: disable=protected-access
    mtype = "L5_TPC:A"
    model = "generic"
    d_generic = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", False, mode="generic"
    )
    assert d_generic == pytest.approx(4.942498757112429)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", False, mode="threshold"
    )
    assert d_threshold == pytest.approx(5.494060358805385)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", False, mode="exact"
    )
    assert d_threshold == pytest.approx(1)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_diameter_power_relation(
            model_params[model][mtype], "basal", False, mode="UNKNOWN"
        )


def test_sample_daughter_diameters(neuron, model_params):
    """Test the _sample_diameter_power_relation function."""
    # pylint: disable=protected-access
    mtype = "L5_TPC:A"
    model = "generic"

    param_tree = {
        "mode_sibling": "threshold",
        "mode_diameter_power_relation": "threshold",
        "with_asymmetry": True,
        "reduction_factor_max": 1.0,
        "neurite_type": "basal",
        "asymmetry_threshold": 1.0,
        "trunk_diam": 2.131444808766456,
        "tot_length": 1085.1178037524223,
        "terminal_diam": 0.35602967635703164,
        "major_sections": [],
    }
    d = build_diameters._sample_daughter_diameters(
        neuron.sections[0], model_params[model][mtype], param_tree
    )
    assert d == pytest.approx([0.9373572556376409, 0.7736807430277025])

    param_tree["with_asymmetry"] = False
    d_no_asymmetry = build_diameters._sample_daughter_diameters(
        neuron.sections[0], model_params[model][mtype], param_tree
    )
    assert d_no_asymmetry == pytest.approx([1.120733059071045, 0.8804980886414063])


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
