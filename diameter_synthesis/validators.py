"""Tools to validate the input parameters."""

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

import json

from jsonschema import validate
from pkg_resources import resource_stream

MODEL_SCHEMA = json.load(resource_stream("diameter_synthesis", "schemas/model_params.json"))
CONFIG_SCHEMA = json.load(resource_stream("diameter_synthesis", "schemas/config.json"))


def validate_model_params(data):
    """Validate model parameter dictionary."""
    validate(data, MODEL_SCHEMA)


def validate_config(data):
    """Validate config dictionary."""
    validate(data, CONFIG_SCHEMA)
