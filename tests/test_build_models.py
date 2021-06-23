"""Test the build_models module."""

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

import pytest

from diameter_synthesis import build_models
from diameter_synthesis.exception import DiameterSynthesisError

from .testing_tools import compare_dicts


def test_build(single_pop, model_params, model_data, empty_build_result):
    """Test the build function."""
    # Test with generic model
    config = {
        "models": ["generic"],
        "neurite_types": ["basal", "apical"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    res = build_models.build(single_pop, config, with_data=True)
    res_models_params = build_models.build(single_pop, config, with_data=False)

    assert len(res) == 2
    compare_dicts(res, [model_params, model_data], precision=3)
    compare_dicts(res_models_params, model_params, precision=3)

    # Test with astrocyte model
    config_astrocyte = {
        "models": ["astrocyte"],
        "neurite_types": ["basal", "apical"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    res_astrocyte = build_models.build(single_pop, config_astrocyte, with_data=True)

    assert len(res_astrocyte) == 2
    compare_dicts(res_astrocyte, [model_params, model_data], precision=3)

    # Test with empty population
    res_empty = build_models.build([], config_astrocyte, with_data=True)
    assert res_empty == empty_build_result

    # Test with unknown model (should raise a DiameterSynthesisError exception)
    bad_config = {
        "models": ["UNKNOWN"],
        "neurite_types": ["basal", "apical"],
        "terminal_threshold": 2,
        "taper": {"max": 1e-6, "min": -0.01},
    }
    with pytest.raises(DiameterSynthesisError):
        build_models.build(single_pop, bad_config, with_data=True)
