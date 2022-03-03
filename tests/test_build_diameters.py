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
        0: [1.3711128234863281, 1.0080387592315674, 0.704459547996521],  # basal
        1: [1.8643510341644287, 1.2897943258285522, 0.9409240484237671],  # basal
        2: [0.6652331352233887, 0.6105318665504456, 0.6113356351852417],  # basal
        3: [0.7050515413284302, 0.6147821545600891, 0.6037149429321289],  # basal
        4: [3.007146120071411, 2.592045307159424, 2.1769447326660156],  # apical
        5: [0.5481259822845459, 0.45678243041038513, 0.44404134154319763],  # apical
        6: [2.1249232292175293, 1.820948839187622, 1.5981926918029785],  # apical
        7: [1.0208232402801514, 0.6448897123336792, 0.48835086822509766],  # apical
        8: [1.011188268661499, 0.7493581771850586, 0.579250693321228],  # apical
        9: [0.4641014039516449, 0.47053343057632446, 0.5290390253067017],  # apical
        10: [0.4879864752292633, 0.5238271355628967, 0.6018319725990295],  # apical
        11: [0.43771225214004517, 0.39942067861557007, 0.3820285201072693],  # apical
        12: [0.4310743808746338, 0.3950076997280121, 0.415347158908844],  # apical
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
        0: [1.4650635719299316, 1.0481635332107544, 0.7526174187660217],  # basal
        1: [1.848629355430603, 1.266541600227356, 0.8283082246780396],  # basal
        2: [0.6392568349838257, 0.6155299544334412, 0.6065273880958557],  # basal
        3: [0.6419328451156616, 0.5688144564628601, 0.5582500100135803],  # basal
        4: [3.5188217163085938, 3.0272274017333984, 2.535632610321045],  # apical
        5: [0.9256387948989868, 0.5633050203323364, 0.4497975707054138],  # apical
        6: [2.048628091812134, 1.9843870401382446, 1.9606387615203857],  # apical
        7: [1.0638917684555054, 0.6236304640769958, 0.4124460816383362],  # apical
        8: [1.2614320516586304, 1.0193065404891968, 0.8004169464111328],  # apical
        9: [0.5992140769958496, 0.49981755018234253, 0.4857638478279114],  # apical
        10: [0.575838565826416, 0.6442711353302002, 0.7361840009689331],  # apical
        11: [0.5330768823623657, 0.437779039144516, 0.41082221269607544],  # apical
        12: [0.5104016661643982, 0.434490442276001, 0.391508013010025],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
        17: [1.926880121231079, 1.6383421421051025, 1.3690083026885986],  # apical
        18: [0.7419866323471069, 0.606209397315979, 0.5698315501213074],  # apical
        19: [1.3392184972763062, 1.0754210948944092, 0.9233089685440063],  # apical
        20: [0.7070451974868774, 0.517382800579071, 0.5084635019302368],  # apical
        21: [0.7122340798377991, 0.6030470132827759, 0.5860720276832581],  # apical
        22: [0.544852614402771, 0.5220872163772583, 0.5411176681518555],  # apical
        23: [0.5539101362228394, 0.5045477151870728, 0.5091089606285095],  # apical
        24: [0.47311487793922424, 0.48335814476013184, 0.5041438937187195],  # apical
        25: [0.4828456938266754, 0.48316454887390137, 0.5292574167251587],  # apical
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
