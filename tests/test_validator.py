"""Test the validators module."""

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
from jsonschema.exceptions import ValidationError

from diameter_synthesis import validators


def test_defaults(config, model_params):
    """Test the validation of JSON schemas."""
    mtype = "L5_TPC:A"
    model = "generic"

    validators.validate_model_params(model_params[model][mtype])
    validators.validate_config(config[model])

    config[model]["features"] = {"FEATURE 1": {"apical": 0.2}}
    validators.validate_config(config[model])

    config[model]["features"] = {"FEATURE 2": {"UNKNOWN": 0.2}}
    with pytest.raises(ValidationError):
        validators.validate_config(config[model])
