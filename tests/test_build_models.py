"""Test the build_models module."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import dictdiffer
import pytest

from diameter_synthesis import build_models
from diameter_synthesis.exception import DiameterSynthesisError


def test_build(
    single_pop, model_params, model_data, empty_build_result, astro_model_params, astro_model_data
):
    """Test the build function."""
    # Test with generic model
    config = {
        "models": ["generic"],
        "neurite_types": ["basal_dendrite", "apical_dendrite"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    res = build_models.build(single_pop, config, with_data=True)
    res_models_params = build_models.build(single_pop, config, with_data=False)

    assert len(res) == 2
    assert dictdiffer.diff(res, [model_params, model_data], absolute_tolerance=1e-3)
    assert dictdiffer.diff(res_models_params, model_params, absolute_tolerance=1e-3)

    # Test with astrocyte model
    config_astrocyte = {
        "models": ["astrocyte"],
        "neurite_types": ["basal_dendrite", "apical_dendrite"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    res_astrocyte = build_models.build(single_pop, config_astrocyte, with_data=True)

    assert len(res_astrocyte) == 2

    assert dictdiffer.diff(
        res_astrocyte, [astro_model_params, astro_model_data], absolute_tolerance=1e-3
    )

    # Test with empty population
    res_empty = build_models.build([], config_astrocyte, with_data=True)
    assert res_empty == empty_build_result

    # Test with unknown model (should raise a DiameterSynthesisError exception)
    bad_config = {
        "models": ["UNKNOWN"],
        "neurite_types": ["basal_dendrite", "apical_dendrite"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    with pytest.raises(DiameterSynthesisError):
        build_models.build(single_pop, bad_config, with_data=True)
