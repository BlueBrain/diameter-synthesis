"""Tools to validate the input parameters."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import json
from importlib import resources

from jsonschema import validate

MODEL_SCHEMA = json.load(
    (resources.files("diameter_synthesis") / "schemas" / "model_params.json").open(encoding="utf-8")
)
CONFIG_SCHEMA = json.load(
    (resources.files("diameter_synthesis") / "schemas" / "config.json").open(encoding="utf-8")
)


def validate_model_params(data):
    """Validate model parameter dictionary."""
    validate(data, MODEL_SCHEMA)


def validate_config(data):
    """Validate config dictionary."""
    validate(data, CONFIG_SCHEMA)
