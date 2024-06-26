"""Test the distribution_fitting module."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import logging

import pytest

from diameter_synthesis import distribution_fitting
from diameter_synthesis.exception import DiameterSynthesisError


def test_fit_distribution(config, caplog):
    """Test the fit_distribution function."""
    all_data = [
        [0.75, 0.25],
        [0.7, 0.5],
        [0.8, 0.001],
        [0.99, 0.033],
        [0.428, 0.254],
        [0.808, 0.006],
        [0.397, 0.537],
        [0.72, 0.223],
        [0.243, 0.02],
        [0.342, 0.244],
    ]
    data_1d = [i[0] for i in all_data]
    config["neurite_type"] = "basal_dendrite"

    # Test expon_rev distribution
    res_expon_rev = distribution_fitting.fit_distribution(data_1d, "expon_rev", extra_params=config)
    assert res_expon_rev == pytest.approx(
        {"loc": -0.99, "scale": 0.3722, "min": 0.28755, "max": 0.99, "num_value": 10.0}
    )

    # Test exponnorm distribution
    res_exponnorm = distribution_fitting.fit_distribution(data_1d, "exponnorm", extra_params=config)

    assert res_exponnorm == pytest.approx(
        {
            "a": 0.3,
            "loc": 0.550945,
            "scale": 0.224617,
            "min": 0.28755,
            "max": 0.9081,
            "num_value": 10,
        }
    )

    # Test gamma distribution
    res_gamma = distribution_fitting.fit_distribution(data_1d, "gamma", extra_params=config)
    assert res_gamma == pytest.approx(
        {"a": 4.0, "loc": -1e-09, "scale": 0.15445, "min": 0.28755, "max": 0.9081, "num_value": 10}
    )

    # Test unknown distribution
    with pytest.raises(DiameterSynthesisError):
        distribution_fitting.fit_distribution(all_data, "UNKNOWN", extra_params=config)

    # Test unknown attribute name
    with pytest.raises(DiameterSynthesisError):
        distribution_fitting.fit_distribution(
            data_1d, "expon_rev", attribute_name="UNKNOWN", extra_params=config
        )

    # Test with empty input data
    caplog.clear()
    caplog.set_level(logging.WARNING)
    config["name"] = "sibling_ratios"
    res_no_data = distribution_fitting.fit_distribution([], "expon_rev", extra_params=config)

    assert res_no_data == {
        "a": 0.0,
        "loc": 0.0,
        "scale": 0.0,
        "min": 0.0,
        "max": 0.1,
        "num_value": 0,
    }

    assert len(caplog.record_tuples) == 1
    for i in caplog.record_tuples:
        module, level, entry = i
        assert module == "diameter_synthesis.distribution_fitting"
        assert level == 30
        assert entry == "Not enough data to fit distribution sibling_ratios with 0 points"

    # Test A_MIN and A_MAX
    A_MIN = distribution_fitting.A_MIN
    A_MAX = distribution_fitting.A_MAX
    try:
        distribution_fitting.A_MIN = 999
        distribution_fitting.A_MAX = 999

        # Test exponnorm distribution
        res_exponnorm_MIN_MAX = distribution_fitting.fit_distribution(
            data_1d, "exponnorm", extra_params=config
        )
        assert res_exponnorm_MIN_MAX == pytest.approx(
            {
                "a": 999,
                "loc": 0.24197457,
                "scale": 0.00037496929,
                "min": 0.28755,
                "max": 0.9081,
                "num_value": 10,
            }
        )

        # Test gamma distribution
        res_gamma_MIN_MAX = distribution_fitting.fit_distribution(
            data_1d, "gamma", extra_params=config
        )
        assert res_gamma_MIN_MAX == pytest.approx(
            {
                "a": 999,
                "loc": -1e-09,
                "scale": 0.000618418,
                "min": 0.28755,
                "max": 0.9081,
                "num_value": 10,
            }
        )
    finally:
        distribution_fitting.A_MIN = A_MIN
        distribution_fitting.A_MAX = A_MAX

    # Test with attribute name
    config["features"] = {"TEST FEATURE": {"basal_dendrite": 0.5}}
    res_expon_rev_feature = distribution_fitting.fit_distribution(
        all_data, "expon_rev", attribute_name="TEST FEATURE", extra_params=config
    )
    assert res_expon_rev_feature == pytest.approx(
        {"loc": -0.99, "scale": 0.354875, "min": 0.27765, "max": 0.99, "num_value": 8.0}
    )


def test_sample_distribution():
    """Test the sample_distribution function."""
    model = {
        "params": {
            "a": 1.5,
            "loc": 1,
            "max": 4,
            "min": 0.5,
            "num_value": 500,
            "scale": 0.5,
        },
        "sequential": None,
    }

    # Test expon_rev distribution
    model["distribution"] = "expon_rev"
    model["params"]["loc"] = -1
    res_expon_rev = distribution_fitting.sample_distribution(model)
    assert res_expon_rev == pytest.approx(0.602062745918445)

    # Test exponnorm distribution
    model["distribution"] = "exponnorm"
    model["loc"] = 1
    res_exponnorm = distribution_fitting.sample_distribution(model)
    assert res_exponnorm == pytest.approx(0.5337329679223535)

    # Test gamma distribution
    model["distribution"] = "gamma"
    model["loc"] = 1
    res_gamma = distribution_fitting.sample_distribution(model)
    assert res_gamma == pytest.approx(0.5690896450847811)

    # Test unknown distribution
    model["distribution"] = "UNKNOWN"
    model["loc"] = 1
    with pytest.raises(DiameterSynthesisError):
        distribution_fitting.sample_distribution(model)

    # Test truncate fail
    model["distribution"] = "expon_rev"
    model["params"]["loc"] = 999
    with pytest.raises(DiameterSynthesisError):
        distribution_fitting.sample_distribution(model)


def test_evaluate_distribution():
    """Test the evaluate_distribution function."""
    # Test expon_rev distribution
    res_expon_rev = distribution_fitting.evaluate_distribution(
        -1, "expon_rev", {"loc": 1, "scale": 10}
    )
    assert res_expon_rev == pytest.approx(0.1)

    # Test exponnorm distribution
    res_expon_rev = distribution_fitting.evaluate_distribution(
        1, "exponnorm", {"a": 1, "loc": 1, "scale": 10}
    )
    assert res_expon_rev == pytest.approx(0.026157829186512344)

    # Test exponnorm distribution
    res_expon_rev = distribution_fitting.evaluate_distribution(
        1, "gamma", {"a": 1, "loc": 1, "scale": 10}
    )
    assert res_expon_rev == pytest.approx(0.1)

    # Test unknown distribution
    with pytest.raises(DiameterSynthesisError):
        distribution_fitting.evaluate_distribution(0, "UNKNOWN", {})
