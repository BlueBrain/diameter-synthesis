"""Tools to validate the input parameters."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

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
