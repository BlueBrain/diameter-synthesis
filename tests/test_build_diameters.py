import logging

import numpy as np
import morphio.mut
import pytest

from diameter_synthesis import build_diameters
from diameter_synthesis.exception import DiameterSynthesisError

from .testing_tools import _compare_diameters
from .testing_tools import reset_random_seed


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


def test_build_no_asymmetry_threshold(config, model_params, neuron):
    """Test main build without seed in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    config[model].pop("asymmetry_threshold")
    with pytest.raises(KeyError):
        build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])


def test_build_multiple_models_warning(config, model_params, neuron, neuron_diametrized, caplog):
    """Test main build without seed in config"""
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
    """Test main build without seed in config"""
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
    """Test main build without seed in config"""
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


def test_build_small_asymmetry_threshold(test_data_path, config, model_params, neuron):
    """Test the main build function to diametrize a neuron with a small asymmetry threshold"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"

    config[model]["asymmetry_threshold"]["basal"] = 0.1

    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    neuron_diametrized = morphio.mut.Morphology(
        test_data_path / "C030796A-P3_lite_diametrized_small_asymmetry_threshold.h5"
    )

    _compare_diameters(neuron_diametrized, neuron)


def test_build_small_n_tries_warning(config, model_params, neuron, caplog):
    """Test main build without seed in config"""
    neurite_types = ["basal", "apical"]
    mtype = "L5_TPC:A"
    model = "generic"
    config[model]["trunk_max_tries"] = 1

    caplog.clear()
    caplog.set_level(logging.WARNING)
    build_diameters.build(neuron, model_params[model][mtype], neurite_types, config[model])

    assert len(caplog.record_tuples) == 9
    neurite_types = [
        "basal",
        "apical",
        "basal",
        "apical",
        "basal",
        "apical",
        "basal",
        "basal",
        "apical",
    ]
    for i in zip(neurite_types, caplog.record_tuples):
        neurite_type, (module, level, entry) = i
        assert module == "diameter_synthesis.build_diameters"
        assert level == 30
        assert entry == f"max tries attained with {neurite_type}"


def test_build_small_trunk_diam_warning(config, model_params, neuron, caplog):
    """Test main build without seed in config"""
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

    assert len(caplog.record_tuples) == 2
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
        model_params[model][mtype], "basal", 0, mode="generic", asymmetry_threshold=0.3
    )
    assert d_generic == pytest.approx(0.8253851329090136)

    d_threshold = build_diameters._sample_sibling_ratio(
        model_params[model][mtype], "basal", 0, mode="threshold", asymmetry_threshold=0.3
    )
    assert d_threshold == pytest.approx(0.7244487906052952)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_sibling_ratio(
            model_params[model][mtype], "basal", 0, mode="UNKNOWN", asymmetry_threshold=0.3
        )


@reset_random_seed
def test_sample_diameter_power_relation(model_params):
    """Test the _sample_diameter_power_relation function"""
    mtype = "L5_TPC:A"
    model = "generic"
    d_generic = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", 0, mode="generic", asymmetry_threshold=0.3
    )
    assert d_generic == pytest.approx(4.942498757112429)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", 0, mode="threshold", asymmetry_threshold=0.3
    )
    assert d_threshold == pytest.approx(5.494060358805385)

    d_threshold = build_diameters._sample_diameter_power_relation(
        model_params[model][mtype], "basal", 0, mode="exact", asymmetry_threshold=0.3
    )
    assert d_threshold == pytest.approx(1)

    with pytest.raises(DiameterSynthesisError):
        build_diameters._sample_diameter_power_relation(
            model_params[model][mtype], "basal", 0, mode="UNKNOWN", asymmetry_threshold=0.3
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
    }
    d = build_diameters._sample_daughter_diameters(
        neuron.sections[0], model_params[model][mtype], param_tree
    )
    assert d == pytest.approx([0.7736807430277025, 0.9373572556376409])

    param_tree["with_asymmetry"] = False
    d_no_asymmetry = build_diameters._sample_daughter_diameters(
        neuron.sections[0], model_params[model][mtype], param_tree
    )
    assert d_no_asymmetry == pytest.approx([0.7732998023091868, 0.8796494557710439])
