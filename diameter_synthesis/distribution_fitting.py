"""Functions to fit distributions to parameters of diameter models"""
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

L = logging.getLogger(__name__)
np.seterr(invalid="ignore", divide="ignore")


def _truncate(sample_func, min_value, max_value):
    """ensure sample is within bounds"""
    sample = sample_func()
    while sample > max_value or sample < min_value:
        sample = sample_func()
    return sample


def fit_distribution(data, distribution, attribute_name=None, extra_params=None):
    """fit a distribution from data"""
    if attribute_name == "asymmetry_threshold":
        tpes = np.asarray(data, dtype=np.float32)[:, 1]
        data = np.asarray(data, dtype=np.float32)[:, 0]
        data = data[tpes < extra_params["threshold"][extra_params["neurite_type"]]]
    elif attribute_name is not None:
        raise DiameterSynthesisError(
            "attribute_name {} not implemented".format(attribute_name)
        )

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

        loc, scale = expon.fit(-data)
        return {
            "loc": loc,
            "scale": scale,
            "min": np.percentile(data, PERCENTILE),
            "max": max(data),
            "num_value": len(data),
        }

    if distribution == "exponnorm":
        from scipy.stats import exponnorm

        var_a, loc, scale = exponnorm.fit(data)
        if var_a > A_MAX:
            var_a, loc, scale = exponnorm.fit(data, f0=A_MAX)
        if var_a < A_MIN:
            var_a, loc, scale = exponnorm.fit(data, f0=A_MIN)
        return {
            "a": var_a,
            "loc": loc,
            "scale": scale,
            "min": np.percentile(data, PERCENTILE),
            "max": np.percentile(data, 100 - PERCENTILE),
            "num_value": len(data),
        }

    if distribution == "gamma":
        from scipy.stats import gamma

        var_a, loc, scale = gamma.fit(data, floc=-1e-9)
        if var_a > A_MAX:
            var_a, loc, scale = gamma.fit(data, f0=A_MAX, floc=-1e-9)
        if var_a < A_MIN:
            var_a, loc, scale = gamma.fit(data, f0=A_MIN, floc=-1e-9)
        return {
            "a": var_a,
            "loc": loc,
            "scale": scale,
            "min": np.percentile(data, PERCENTILE),
            "max": np.percentile(data, 100 - PERCENTILE),
            "num_value": len(data),
        }

    raise DiameterSynthesisError("Distribution not understood")


def sample_distribution(model):
    """sample from a distribution"""
    if "a" in model["params"]:
        a_clip = np.clip(model["params"]["a"], A_MIN, A_MAX)

    if model["distribution"] == "expon_rev":
        from scipy.stats import expon

        return _truncate(
            lambda: -expon.rvs(model["params"]["loc"], model["params"]["scale"]),
            model["params"]["min"],
            model["params"]["max"],
        )

    if model["distribution"] == "exponnorm":
        from scipy.stats import exponnorm

        return _truncate(
            lambda: exponnorm.rvs(
                a_clip, model["params"]["loc"], model["params"]["scale"]
            ),
            model["params"]["min"],
            model["params"]["max"],
        )

    if model["distribution"] == "gamma":
        from scipy.stats import gamma

        return _truncate(
            lambda: gamma.rvs(a_clip, model["params"]["loc"], model["params"]["scale"]),
            model["params"]["min"],
            model["params"]["max"],
        )

    raise DiameterSynthesisError("Distribution not understood")


def evaluate_distribution(value, distribution, params):
    """evaluate the fit of a distribution"""
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
