""" Functions to fit distributions to parameters of diameter models """
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.stats import expon, exponnorm, gamma, skewnorm

import diameter_synthesis.utils as utils
from diameter_synthesis.utils import A_MAX, A_MIN, MIN_DATA_POINTS, ROUND

##################################
# Distribution related functions #
##################################

N_BINS = 10
PERCENTILE = 5
MIN_SAMPLE_NUM = 10

np.seterr(invalid="ignore", divide="ignore")


def build_spline(val_x, val_y, weights):
    """ build a spline model and return parameters"""
    spl = UnivariateSpline(
        val_x,
        val_y,
        w=np.array(weights) / np.sum(weights),
        s=len(weights) * utils.SPLINE_SMOOTH,
        k=3,
    )
    return spl._eval_args  # pylint: disable=protected-access


def evaluate_spline(val_x, tck):
    """ evaluate a spline model from parameters"""
    spl = UnivariateSpline._from_tck(tck)  # pylint: disable=protected-access

    return spl(val_x)


def evaluate_distribution(val_x, distribution, params):
    """ evaluate the fit of a distribution"""

    if distribution == "expon_rev":

        return expon.pdf(-val_x, params["loc"], params["scale"])

    if distribution == "exponnorm":

        return exponnorm.pdf(val_x, params["a"], params["loc"], params["scale"])

    if distribution == "gamma":

        return gamma.pdf(val_x, params["a"], params["loc"], params["scale"])

    if distribution == "skewnorm":

        return skewnorm.pdf(val_x, params["a"], params["loc"], params["scale"])

    raise Exception("Distribution not understood")


def truncate(sample_func, min_value, max_value):
    """ truncatet a sampled distribution """
    sample = sample_func()

    while sample > max_value or sample < min_value:
        sample = sample_func()

    return sample


def sample_distribution_single(model):
    """ sample from a distribution (no slicing)"""
    params = model["params"]

    if model["distribution"] == "expon_rev":

        return truncate(
            lambda: -expon.rvs(params["loc"], params["scale"]),
            params["min"],
            params["max"],
        )

    if model["distribution"] == "exponnorm":

        return truncate(
            lambda: exponnorm.rvs(
                np.clip(params["a"], A_MIN, A_MAX), params["loc"], params["scale"]
            ),
            params["min"],
            params["max"],
        )

    if model["distribution"] == "gamma":

        return truncate(
            lambda: gamma.rvs(
                np.clip(params["a"], A_MIN, A_MAX), params["loc"], params["scale"]
            ),
            params["min"],
            params["max"],
        )

    if model["distribution"] == "skewnorm":

        return truncate(
            lambda: skewnorm.rvs(
                np.clip(params["a"], A_MIN, A_MAX), params["loc"], params["scale"]
            ),
            params["min"],
            params["max"],
        )

    raise Exception("Distribution not understood")


def sample_distribution(model, tpe=0):
    """ sample from a distribution"""
    if (
        not isinstance(model["sequential"], str)
        or model["sequential"] == "asymmetry_threshold"
    ):  # if no sequential fitting needed
        return sample_distribution_single(model)

    # if our sample is out of the bins used for fitting, restrict to these bins
    tpes = list(model["params"]["params_data"].keys())
    tpe_min = float(tpes[0])
    tpe_max = float(tpes[-1])
    tpe = np.clip(tpe, tpe_min, tpe_max)

    model_tpe = {}
    model_tpe["params"] = {}
    model_tpe["distribution"] = model["distribution"]
    try:
        model_tpe["params"]["a"] = evaluate_spline(tpe, model["params"]["a"])
    except BaseException:  # pylint: disable=broad-except
        pass
    model_tpe["params"]["loc"] = evaluate_spline(tpe, model["params"]["loc"])
    model_tpe["params"]["scale"] = evaluate_spline(tpe, model["params"]["scale"])
    model_tpe["params"]["min"] = evaluate_spline(tpe, model["params"]["min"])
    model_tpe["params"]["max"] = evaluate_spline(tpe, model["params"]["max"])

    try:
        return sample_distribution_single(model_tpe)
    except BaseException:  # pylint: disable=broad-except
        raise Exception("error in parameters for tpe ", tpe, " with model ", model_tpe)


def fit_distribution_single(data, distribution):
    """ generic function to fit a distribution with scipy (single slice)"""
    if len(data) > MIN_DATA_POINTS:

        if distribution == "expon_rev":
            loc, scale = expon.fit(-data)

            return {
                "loc": np.round(loc, ROUND),
                "scale": np.round(scale, ROUND),
                "min": np.round(np.percentile(data, PERCENTILE), ROUND),
                "max": np.round(max(data), ROUND),
                "num_value": len(data),
            }

        if distribution == "exponnorm":
            var_a, loc, scale = exponnorm.fit(data)
            # refit if we get crazy values for a
            if var_a > A_MAX:
                var_a, loc, scale = exponnorm.fit(data, f0=A_MAX)
            if var_a < A_MIN:
                var_a, loc, scale = exponnorm.fit(data, f0=A_MIN)

            return {
                "a": np.round(var_a, ROUND),
                "loc": np.round(loc, ROUND),
                "scale": np.round(scale, ROUND),
                "min": np.round(np.percentile(data, PERCENTILE), ROUND),
                "max": np.round(np.percentile(data, 100 - PERCENTILE), ROUND),
                "num_value": len(data),
            }

        if distribution == "skewnorm":
            var_a, loc, scale = skewnorm.fit(data)
            # refit if we get crazy values for a
            if var_a > A_MAX:
                var_a, loc, scale = skewnorm.fit(data, f0=A_MAX)
            if var_a < A_MIN:
                var_a, loc, scale = skewnorm.fit(data, f0=A_MIN)

            return {
                "a": np.round(var_a, ROUND),
                "loc": np.round(loc, ROUND),
                "scale": np.round(scale, ROUND),
                "min": np.round(np.percentile(data, PERCENTILE), ROUND),
                "max": np.round(np.percentile(data, 100 - PERCENTILE), ROUND),
                "num_value": len(data),
            }

        if distribution == "gamma":
            var_a, loc, scale = gamma.fit(data, floc=-1e-9)
            # refit if we get crazy values for a
            if var_a > A_MAX:
                var_a, loc, scale = gamma.fit(data, f0=A_MAX, floc=-1e-9)
            if var_a < A_MIN:
                var_a, loc, scale = gamma.fit(data, f0=A_MIN, floc=-1e-9)

            return {
                "a": np.round(var_a, ROUND),
                "loc": np.round(loc, ROUND),
                "scale": np.round(scale, ROUND),
                "min": np.round(np.percentile(data, PERCENTILE), ROUND),
                "max": np.round(np.percentile(data, 100 - PERCENTILE), ROUND),
                "num_value": len(data),
            }

        raise Exception("Distribution not understood")
    # if no data, return null parameters (for neurons without apical dentrites)
    return {
        "a": 0.0,
        "loc": 0.0,
        "scale": 0.0,
        "min": 0.0,
        "max": 0.1,
        "num_value": len(data),
    }


def fit_distribution(data, distribution, seq=None, extra_params=None):
    """ generic function to fit a distribution with scipy """

    threshold = extra_params["threshold"][extra_params["neurite_type"]]

    if not isinstance(seq, str):  # if no sequential fitting needed
        return fit_distribution_single(data, distribution)

    if seq == "asymmetry_threshold":

        if len(data) > 0:
            # here, enforcing a flat type is necessary, for some weird reason!
            tpes = np.asarray(data, dtype=np.float32)[:, 1]  # collect the type of point
            values = np.asarray(data, dtype=np.float32)[:, 0]  # collect the data itself

            values = values[tpes < threshold]

            return fit_distribution_single(values, distribution)
        return 0

    if len(data) > 0:
        tpes = np.asarray(data)[:, 1]  # collect the type of point
        values = np.asarray(data)[:, 0]  # collect the data itself

        # set the bins for estimating parameters if we can otherwise use two bins
        # to be able to fit later
        bins, _ = utils.set_bins(tpes, N_BINS, n_min=MIN_SAMPLE_NUM)

        params = {}
        for i in range(len(bins) - 1):
            data_tpe = values[
                (tpes >= bins[i]) & (tpes < bins[i + 1])
            ]  # select the values by its type
            params[
                np.round((bins[i + 1] + bins[i]) / 2.0, ROUND)
            ] = fit_distribution_single(data_tpe, distribution)
        return update_params_fit_distribution(params)

    return {
        "a": 0.0,
        "loc": 0.0,
        "scale": 0.0,
        "min": 0.0,
        "max": 0.1,
        "num_value": len(data),
        "params_data": {
            "a": 0.0,
            "loc": 0.0,
            "scale": 0.0,
            "min": 0.0,
            "max": 0.1,
            "num_value": len(data),
        },
    }


def update_params_fit_distribution(params_data):
    """ linear fit to model parameters as a function of a given quantity tpes_model
    and update the model dictionary with the fits of parameters """

    # update the parameters for evaluation with other values of tpes
    params = {}
    params["params_data"] = params_data  # save the data params for plotting

    if (
        len(params_data) > 1
    ):  # only try that if we have a sequence of more than two fits
        tpes_model = list(params_data)  # fancy python3 way to get dict.keys()
        params_values = list(params_data.values())

        try:
            var_as = np.array([v["a"] for v in params_values])
            # prevent large or small values of a from bad fits
            var_as = np.clip(var_as, A_MIN, A_MAX)
        except BaseException:  # pylint: disable=broad-except
            pass
        locs = np.array([v["loc"] for v in params_values])
        scales = np.array([v["scale"] for v in params_values])
        mins = np.array([v["min"] for v in params_values])
        maxs = np.array([v["max"] for v in params_values])
        weights = np.array([v["num_value"] for v in params_values])

        try:
            params["a"] = build_spline(tpes_model, var_as, weights)
        except BaseException:  # pylint: disable=broad-except
            pass
        params["loc"] = build_spline(tpes_model, locs, weights)
        params["scale"] = build_spline(tpes_model, scales, weights)
        params["min"] = build_spline(tpes_model, mins, weights)
        params["max"] = build_spline(tpes_model, maxs, weights)
        params["num_value"] = np.sum(weights)

    else:
        params["a"] = [0.0]
        params["loc"] = [0.0]
        params["scale"] = [0.0]
        params["min"] = 0.0
        params["max"] = 0.1
        params["num_value"] = 0

    return params
