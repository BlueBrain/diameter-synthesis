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

import morphio.mut

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


def test_build_diameters(simpler_config, simpler_model_params, small_morph, test_data_path):
    """Test the main build function to diametrize a neuron."""
    neurite_types = ["basal_dendrite", "apical_dendrite"]

    build_diameters.build(small_morph, neurite_types, simpler_model_params, simpler_config)
    expected = morphio.mut.Morphology(test_data_path / "simpler_morph_diametrized.asc")
    _compare_diameters(expected, small_morph)
