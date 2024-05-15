"""Test the validators module."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import pytest
from jsonschema.exceptions import ValidationError

from diameter_synthesis import validators


def test_defaults(config, model_params):
    """Test the validation of JSON schemas."""
    validators.validate_model_params(model_params)
    validators.validate_config(config)

    config["features"] = {"FEATURE 1": {"apical_dendrite": 0.2}}
    validators.validate_config(config)

    config["features"] = {"FEATURE 2": {"UNKNOWN": 0.2}}
    with pytest.raises(ValidationError):
        validators.validate_config(config)
