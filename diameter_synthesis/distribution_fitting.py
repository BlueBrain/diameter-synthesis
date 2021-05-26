"""Functions to fit distributions to parameters of diameter models. TO REVIEW!."""
import logging

import numpy as np

from diameter_synthesis.exception import DiameterSynthesisError

# pylint: disable=import-outside-toplevel

MIN_DATA_POINTS = 1  # minimum number of points to fit a distribution
A_MAX = 4
A_MIN = 0.3

N_BINS = 10
PERCENTILE = 5
MIN_SAMPLE_NUM = 10

MAX_TRUNCATE_TRIES = 100

L = logging.getLogger(__name__)
np.seterr(invalid="ignore", divide="ignore")


def _truncate(sample_func, min_value, max_value):
    """Ensure sample is within bounds."""
    sample = sample_func()
    n_tries = 0
    while sample > max_value or sample < min_value:
        sample = sample_func()
        n_tries += 1
        if n_tries >= MAX_TRUNCATE_TRIES:
            raise DiameterSynthesisError(
                f"Could not truncate the sample between {min_value} and {max_value} "
                f"(the last value was {sample})"
            )
    return sample


def fit_distribution(all_data, distribution, attribute_name=None, extra_params=None):
    """Fit a distribution from data.

    Args:
        all_data (list/array): list of data points to fit a distribution to
        distribution (str): Distribution name
        attribute_name (str): Name of additional attribute to fit
        extra_params (dict): Possible additional parameters for the fit

    Returns:
        dict: parameters of the fit
    """
    if attribute_name is not None:
        if extra_params is None:
            raise DiameterSynthesisError(
                "'extra_params' can not be None when 'attribute_name' is not None"
            )
        threshold = (
            extra_params.get("features", {})
            .get(attribute_name, {})
            .get(extra_params["neurite_type"])
        )
        if threshold is None:
            raise DiameterSynthesisError(
                "The threshold could not be retrived from the config in 'features'->"
                f"'{attribute_name}'->'{extra_params['neurite_type']}'"
            )
        attribute = np.asarray(all_data, dtype=np.float32)[:, 1]
        data = np.asarray(all_data, dtype=np.float32)[:, 0]
        data = data[attribute < threshold]
    else:
        data = all_data

    if len(data) < MIN_DATA_POINTS:
        L.warning(
            "Not enough data to fit distribution %s with %s points",
            extra_params["name"],
            len(data),
        )
        return {
            "a": 0.0,
            "loc": 0.0,
            "scale": 0.0,
            "min": 0.0,
            "max": 0.1,
            "num_value": len(data),
        }

    if distribution == "expon_rev":
        from scipy.stats import expon

        loc, scale = expon.fit(-np.array(data))
        return {
            "loc": float(loc),
            "scale": float(scale),
            "min": float(np.percentile(data, PERCENTILE)),
            "max": float(max(data)),
            "num_value": float(len(data)),
        }

    if distribution == "exponnorm":
        from scipy.stats import exponnorm

        var_a, loc, scale = exponnorm.fit(data)
        if var_a > A_MAX:
            var_a, loc, scale = exponnorm.fit(data, f0=A_MAX)
        if var_a < A_MIN:
            var_a, loc, scale = exponnorm.fit(data, f0=A_MIN)
        return {
            "a": float(var_a),
            "loc": float(loc),
            "scale": float(scale),
            "min": float(np.percentile(data, PERCENTILE)),
            "max": float(np.percentile(data, 100 - PERCENTILE)),
            "num_value": int(len(data)),
        }

    if distribution == "gamma":
        from scipy.stats import gamma

        var_a, loc, scale = gamma.fit(data, floc=-1e-9)
        if var_a > A_MAX:
            var_a, loc, scale = gamma.fit(data, f0=A_MAX, floc=-1e-9)
        if var_a < A_MIN:
            var_a, loc, scale = gamma.fit(data, f0=A_MIN, floc=-1e-9)
        return {
            "a": float(var_a),
            "loc": float(loc),
            "scale": float(scale),
            "min": float(np.percentile(data, PERCENTILE)),
            "max": float(np.percentile(data, 100 - PERCENTILE)),
            "num_value": int(len(data)),
        }

    raise DiameterSynthesisError("Distribution not understood")


def sample_distribution(model, rng=np.random):
    """Sample from a distribution."""
    if "a" in model["params"]:
        a_clip = np.clip(model["params"]["a"], A_MIN, A_MAX)

    if model["distribution"] == "constant":
        return model["params"]["value"]

    if model["distribution"] == "expon_rev":
        return _truncate(
            lambda: -(
                model["params"]["loc"] + model["params"]["scale"] * rng.standard_exponential()
            ),
            model["params"]["min"],
            model["params"]["max"],
        )

    if model["distribution"] == "exponnorm":
        return _truncate(
            lambda: (
                model["params"]["loc"]
                + model["params"]["scale"]
                * (rng.standard_exponential() * a_clip + rng.standard_normal())
            ),
            model["params"]["min"],
            model["params"]["max"],
        )

    if model["distribution"] == "gamma":
        return _truncate(
            lambda: model["params"]["loc"] + model["params"]["scale"] * rng.standard_gamma(a_clip),
            model["params"]["min"],
            model["params"]["max"],
        )

    raise DiameterSynthesisError("Distribution not understood")


def evaluate_distribution(value, distribution, params):
    """Evaluate the fit of a distribution."""
    if distribution == "expon_rev":
        from scipy.stats import expon

        return expon.pdf(-value, params["loc"], params["scale"])

    if distribution == "exponnorm":
        from scipy.stats import exponnorm

        return exponnorm.pdf(value, params["a"], params["loc"], params["scale"])

    if distribution == "gamma":
        from scipy.stats import gamma

        return gamma.pdf(value, params["a"], params["loc"], params["scale"])

    raise DiameterSynthesisError("Distribution not understood")
