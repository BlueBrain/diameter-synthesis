import logging

import numpy as np
import morphio.mut
import pytest
from morphio import PointLevel
from morphio import SectionType

from diameter_synthesis import build_diameters
from diameter_synthesis.exception import DiameterSynthesisError

from .testing_tools import _compare_diameters
from .testing_tools import reset_random_seed


def test_build_simple(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    build_diameters.build(small_morph, model_params[model][mtype], neurite_types, config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.12723756, 0.8501913, 0.67725074],  # basal
        1: [1.88172221, 1.60697436, 1.33222628],  # basal
        2: [0.7807585, 0.5290243, 0.52312386],  # basal
        3: [0.80801165, 0.62176675, 0.59818721],  # basal
        4: [3.19154882, 2.68855309, 2.1855576],  # apical
        5: [1.23321462, 0.85168982, 0.63268346],  # apical
        6: [1.1493628, 0.83900893, 0.72536111],  # apical
        7: [0.52783465, 0.44990396, 0.42342067],  # apical
        8: [0.51564646, 0.40626049, 0.37496084],  # apical
        9: [0.36634988, 0.40507823, 0.45199862],  # apical
        10: [0.36634988, 0.40620297, 0.47370234],  # apical
        11: [0.41420197, 0.42738733, 0.46979269],  # apical
        12: [0.41356379, 0.38290364, 0.389568],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build_simple_with_apical_sections(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron with a known apical section"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    model_params[model][mtype]["apical_point_sec_ids"] = [6]

    build_diameters.build(small_morph, model_params[model][mtype], neurite_types, config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.57616556, 1.10349214, 0.71787846],  # basal
        1: [1.9246546, 1.59596884, 1.29954779],  # basal
        2: [0.8029623, 0.5610494, 0.53921282],  # basal
        3: [0.92826194, 0.69051898, 0.61900127],  # basal
        4: [3.54046178, 3.02864218, 2.51682234],  # apical
        5: [0.49307173, 0.38693029, 0.38693029],  # apical
        6: [2.51682234, 2.15803671, 1.8276546],  # apical
        7: [0.90557307, 0.65467209, 0.48454937],  # apical
        8: [1.32304311, 1.06992745, 0.8688401],  # apical
        9: [0.71307242, 0.5796622, 0.47997007],  # apical
        10: [0.68834245, 0.61472285, 0.61060423],  # apical
        11: [0.51838523, 0.38276276, 0.38276276],  # apical
        12: [0.46621996, 0.42117757, 0.44140521],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build_simple_with_several_apical_sections(config, model_params, small_morph):
    """Test the main build function to diametrize a neuron with a known apical section"""
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

    build_diameters.build(small_morph, model_params[model][mtype], neurite_types, config[model])

    expected = morphio.mut.Morphology(small_morph)
    diameters = {
        0: [1.39890289, 0.94757283, 0.60740244],  # basal
        1: [1.97971463, 1.43837976, 1.02925932],  # basal
        2: [0.81800926, 0.58568162, 0.54873765],  # basal
        3: [0.73983711, 0.59171486, 0.59403551],  # basal
        4: [3.71459246, 3.34208155, 2.96956968],  # apical
        5: [0.44482154, 0.3843236, 0.37343067],  # apical
        6: [2.96956968, 2.40652895, 1.8434881],  # apical
        7: [1.08870435, 0.77451909, 0.49428192],  # apical
        8: [1.25662541, 0.99705613, 0.84291184],  # apical
        9: [0.63015556, 0.48776197, 0.44297838],  # apical
        10: [0.62148929, 0.47264618, 0.45764357],  # apical
        11: [0.41042924, 0.38227296, 0.41449428],  # apical
        12: [0.40574685, 0.35985672, 0.3745116],  # apical
        13: [0.2, 0.2, 0.2],  # axon
        14: [0.2, 0.2, 0.2],  # axon
        15: [0.2, 0.2, 0.2],  # basal
        16: [0.2, 0.2, 0.2],  # basal
        17: [1.53278005, 1.13844323, 0.88405502],  # apical
        18: [0.66239756, 0.52660555, 0.52660555],  # apical
        19: [0.88405502, 0.70394671, 0.66747582],  # apical
        20: [0.66346329, 0.53578627, 0.53427696],  # apical
        21: [0.59090728, 0.49942747, 0.46960464],  # apical
        22: [0.46960464, 0.46642455, 0.47136903],  # apical
        23: [0.46960464, 0.42222959, 0.41856867],  # apical
        24: [0.38936657, 0.38004285, 0.38004285],  # apical
        25: [0.38936657, 0.45000905, 0.51071393],  # apical
    }
    for section_id, section in expected.sections.items():
        section.diameters = diameters[section_id]

    _compare_diameters(expected, small_morph)


def test_build(config, model_params, neuron, neuron_diametrized):
    """Test the main build function to diametrize a neuron"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_seed(config, model_params, neuron, neuron_diametrized):
    """Test main build without seed in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    seed = config[model].pop("seed")
    np.random.seed(seed)
    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    _compare_diameters(neuron_diametrized, neuron)


def test_build_multiple_models_warning(config, model_params, neuron, neuron_diametrized, caplog):
    """Test main build with multiple models in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["models"].append("SKIPPED MODEL")

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    module, level, entry = caplog.record_tuples[0]
    assert module == "diameter_synthesis.build_diameters"
    assert level == 30
    assert entry == "Several models provided, we will only use the first"

    _compare_diameters(neuron_diametrized, neuron)


def test_build_one_sample(test_data_path, config, model_params, neuron):
    """Test main build with 1 sample in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["n_samples"] = 1

    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_1_sample.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_no_sequential(test_data_path, config, model_params, neuron):
    """Test main build with no sequential in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "astrocyte"

    config[model] = config["generic"]
    config[model]["models"] = [model]
    config[model]["n_samples"] = 1
    model_params[model] = model_params["generic"]
    model_params[model][mtype]["sibling_ratios"]["apical"]["sequential"] = None

    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_astrocyte.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_small_n_tries_warning(config, model_params, neuron, caplog):
    """Test main build with trunk_max_tries = 1"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["trunk_max_tries"] = 1

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    assert len(caplog.record_tuples) == 8
    neurite_types = [
        "apical",
        "apical",
        "basal",
        "apical",
        "apical",
        "apical",
        "basal",
        "apical",
    ]
    for i in zip(neurite_types, caplog.record_tuples):
        neurite_type, (module, level, entry) = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == f"max tries attained with {neurite_type}"


def test_build_small_trunk_diam_warning(config, model_params, neuron, caplog):
    """Test trunk diam warning in main build"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    TRUNK_FRAC_DECREASE = build_diameters.TRUNK_FRAC_DECREASE
    try:
        build_diameters.TRUNK_FRAC_DECREASE = 999
        caplog.clear()
        caplog.set_level(logging.WARNING)
        build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])
    finally:
        build_diameters.TRUNK_FRAC_DECREASE = TRUNK_FRAC_DECREASE

    assert len(caplog.record_tuples) == 5
    for i in caplog.record_tuples:
        module, level, entry = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == "sampled trunk diameter < 0.01, so use 1 instead"


def test_select_model():
    """Test the _select_model function"""
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


@reset_random_seed
def test_sample_sibling_ratio(model_params):
    """Test the _sample_sibling_ratio function"""
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


@reset_random_seed
def test_sample_diameter_power_relation(model_params):
    """Test the _sample_diameter_power_relation function"""
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


@reset_random_seed
def test_sample_daughter_diameters(neuron, model_params):
    """Test the _sample_diameter_power_relation function"""
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


def test_diametrize_axon(config, model_params, small_morph):
    """Test the axon diametrizer"""

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
