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
