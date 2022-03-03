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
        0: [1.3947613, 0.8297218, 0.55475914],  # basal
        1: [1.4258235, 1.1060853, 0.9514928],  # basal
        2: [0.67992866, 0.5844326, 0.5498939],  # basal
        3: [0.68993413, 0.55794257, 0.5157046],  # basal
        4: [4.201871, 3.6724274, 3.142984],  # apical
        5: [0.5332143, 0.40053114, 0.36231595],  # apical
        6: [3.142984, 2.6174855, 2.1424375],  # apical
        7: [1.2872376, 0.79716194, 0.5337995],  # apical
        8: [1.3557708, 1.0472242, 0.78575134],  # apical
        9: [0.6548395, 0.45775166, 0.46969777],  # apical
        10: [0.64654696, 0.71682185, 0.8162918],  # apical
        11: [0.59667104, 0.48239225, 0.48059732],  # apical
        12: [0.6458615, 0.5328817, 0.5050491],  # apical
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
        0: [1.5356907, 1.0041976, 0.66725886],  # basal
        1: [1.5833796, 1.1596994, 0.8448763],  # basal
        2: [0.59417355, 0.48759967, 0.481851],  # basal
        3: [0.5885381, 0.53938335, 0.51282257],  # basal
        4: [3.1666424, 2.6880143, 2.2228446],  # apical
        5: [0.51484984, 0.39293796, 0.35675287],  # apical
        6: [2.2228446, 1.7585274, 1.3411181],  # apical
        7: [0.91158104, 0.66549647, 0.4635315],  # apical
        8: [1.0329021, 0.8817488, 0.75879776],  # apical
        9: [0.665402, 0.48054734, 0.4159872],  # apical
        10: [0.63206977, 0.5530502, 0.544706],  # apical
        11: [0.48506194, 0.4511692, 0.50656253],  # apical
        12: [0.48876578, 0.43318385, 0.42985684],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
        17: [1.9270298, 1.5110503, 1.1639264],  # apical
        18: [0.6205703, 0.56107974, 0.5355403],  # apical
        19: [1.1639264, 1.0494831, 1.0048978],  # apical
        20: [0.85187876, 0.61111706, 0.57416606],  # apical
        21: [0.83893585, 0.7002033, 0.65131676],  # apical
        22: [0.571281, 0.57256776, 0.6110884],  # apical
        23: [0.59457844, 0.54023635, 0.53511035],  # apical
        24: [0.52277744, 0.55350095, 0.58257633],  # apical
        25: [0.5206192, 0.49911612, 0.50282],  # apical
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
