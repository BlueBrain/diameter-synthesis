import os
import glob
import shutil
import json
from tqdm import tqdm

import numpy as np
from numpy.polynomial import polynomial as polynomial
from scipy.interpolate import UnivariateSpline

import neurom as nm
from neurom.core import iter_sections

import diameter_synthesis.utils as utils
from diameter_synthesis.utils import get_diameters, set_diameters, ROUND, MIN_DATA_POINTS, A_MAX, A_MIN

####################################
## Distribution related functions ##
####################################


def build_spline(x, y, w):
    """ build a spline model and return parameters"""
    spl = UnivariateSpline(x, y) #, w=np.array(w) / np.sum(w), s=len(w) * utils.SPLINE_SMOOTH, k=3)
    return spl._eval_args


def evaluate_spline(x, tck):
    """ evaluate a spline model from parameters"""
    spl = UnivariateSpline._from_tck(tck)
    return spl(x)


def evaluate_distribution(x, distribution, params):
    """ evaluate the fit of a distribution"""

    if distribution == 'expon_rev':
        from scipy.stats import expon

        return expon.pdf(1 - x, params['loc'], params['scale'])

    elif distribution == 'exponnorm':
        from scipy.stats import exponnorm

        return exponnorm.pdf(x, params['a'], params['loc'], params['scale'])

    elif distribution == 'gamma':
        from scipy.stats import gamma

        return gamma.pdf(x, params['a'], params['loc'], params['scale'])

    elif distribution == 'skewnorm':
        from scipy.stats import skewnorm

        return skewnorm.pdf(x, params['a'], params['loc'], params['scale'])

    else:
        raise Exception('Distribution not understood')


def truncate(sample_func, min_value, max_value):
    """ truncatet a sampled distribution """
    sample = sample_func()

    while sample > max_value or sample < min_value:
        sample = sample_func()

    return sample


def sample_distribution_single(model):
    """ sample from a distribution (no slicing)"""
    params = model['params']

    if model['distribution'] == 'expon_rev':
        from scipy.stats import expon

        return truncate(lambda: 1. - expon.rvs(params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'exponnorm':
        from scipy.stats import exponnorm

        return truncate(lambda: exponnorm.rvs(np.clip(params['a'], A_MIN, A_MAX), params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'gamma':
        from scipy.stats import gamma

        return truncate(lambda: gamma.rvs(np.clip(params['a'], A_MIN, A_MAX), params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'skewnorm':
        from scipy.stats import skewnorm

        return truncate(lambda: skewnorm.rvs(np.clip(params['a'], A_MIN, A_MAX), params['loc'], params['scale']), params['min'], params['max'])

    else:
        raise Exception('Distribution not understood')


def sample_distribution(model, tpe=0):
    """ sample from a distribution"""
    if not isinstance(model['sequential'], str) or model['sequential'] == 'asymmetry_threshold':  # if no sequential fitting needed
        return sample_distribution_single(model)

    else:

        # if our sample is out of the bins used for fitting, restrict to these bins
        tpes = [*model['params']['params_data'].keys()]
        tpe_min = float(tpes[0])
        tpe_max = float(tpes[-1])
        tpe = np.clip(tpe, tpe_min, tpe_max)

        model_tpe = {}
        model_tpe['params'] = {}
        model_tpe['distribution'] = model['distribution']
        try:
            model_tpe['params']['a'] = evaluate_spline(tpe, model['params']['a'])
        except:
            pass
        model_tpe['params']['loc'] = evaluate_spline(tpe, model['params']['loc'])
        model_tpe['params']['scale'] = evaluate_spline(tpe, model['params']['scale'])
        model_tpe['params']['min'] = evaluate_spline(tpe, model['params']['min'])
        model_tpe['params']['max'] = evaluate_spline(tpe, model['params']['max'])

        # TODO: update the hack to use the all data values if the fit failed (only basal)
        """
        if [*params.values()][0][0] == 0. or a < 0 or loc <0  or scale <0 or Min < 0 or Max<0:
            try:
                with open('./model_params_all.json', 'r') as f:
                    params_all = json.load(f)
                params_all = params_all['all_types']['M0']['trunk_diameter']['basal']['params']
            except Exception as e:
                print(e)
                raise Exception('Could not load params_all.json, please create it or have more cells in each class.')

            a = evaluate_spline(tpe, params_all['a'])
            loc = evaluate_spline(tpe, params_all['loc'])
            scale = evaluate_spline(tpe, params_all['scale'])
            Min = evaluate_spline(tpe, params_all['min'])
            Max = evaluate_spline(tpe, params_all['max'])
        """
        try:
            return sample_distribution_single(model_tpe)
        except:
            raise Exception('error in parameters for tpe ', tpe, ' with model ', model_tpe)


def fit_distribution_single(data, distribution, p=5):
    """ generic function to fit a distribution with scipy (single slice)"""
    if len(data) > MIN_DATA_POINTS:

        if distribution == 'expon_rev':
            from scipy.stats import expon
            loc, scale = expon.fit(1 - np.array(data))

            return {'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': 1, 'num_value': len(data)}

        elif distribution == 'exponnorm':
            from scipy.stats import exponnorm
            a, loc, scale = exponnorm.fit(data)
            # refit if we get crazy values for a
            if a > A_MAX:
                a, loc, scale = exponnorm.fit(data, f0=A_MAX)
            if a < A_MIN:
                a, loc, scale = exponnorm.fit(data, f0=A_MIN)

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100 - p), ROUND), 'num_value': len(data)}

        elif distribution == 'skewnorm':
            from scipy.stats import skewnorm
            a, loc, scale = skewnorm.fit(data)
            # refit if we get crazy values for a
            if a > A_MAX:
                a, loc, scale = skewnorm.fit(data, f0=A_MAX)
            if a < A_MIN:
                a, loc, scale = skewnorm.fit(data, f0=A_MIN)

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100 - p), ROUND), 'num_value': len(data)}

        elif distribution == 'gamma':
            from scipy.stats import gamma

            a, loc, scale = gamma.fit(data, floc=-1e-9)
            # refit if we get crazy values for a
            if a > A_MAX:
                a, loc, scale = gamma.fit(data, f0=A_MAX, floc=-1e-9)
            if a < A_MIN:
                a, loc, scale = gamma.fit(data, f0=A_MIN, floc=-1e-9)

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100 - p), ROUND), 'num_value': len(data)}

        else:
            raise Exception('Distribution not understood')
    else:
        # if no data, return null parameters (for neurons without apical dentrites)
        return {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1, 'num_value': len(data)}


def fit_distribution(data, distribution, min_sample_num=20, p=5, n_bins=10, seq=None, extra_params=None, name = 'test', threshold = 1.):
    """ generic function to fit a distribution with scipy """
    if not isinstance(seq, str):  # if no sequential fitting needed
        return fit_distribution_single(data, distribution, p=p)

    elif seq == 'asymmetry_threshold':

        if len(data)>0:
            tpes = np.asarray(data)[:, 1]  # collect the type of point
            values = np.asarray(data)[:, 0]  # collect the data itself

            values = values[tpes < threshold]

            return fit_distribution_single(values, distribution, p=p)
        else:
            return 0

    elif len(data) > 0:
        tpes = np.asarray(data)[:, 1]  # collect the type of point
        values = np.asarray(data)[:, 0]  # collect the data itself

        # set the bins for estimating parameters if we can otherwise use two bins to be able to fit later
        bins, num_values = utils.set_bins(tpes, n_bins, n_min=min_sample_num)

        params = {}
        for i in range(len(bins) - 1):
            data_tpe = values[(tpes >= bins[i]) & (tpes < bins[i + 1])]  # select the values by its type
            params[np.round((bins[i + 1] + bins[i]) / 2., ROUND)] = fit_distribution_single(data_tpe, distribution, p=p)
        return update_params_fit_distribution(params, orders=extra_params['orders'])
    else:
        return {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1, 'num_value': len(data), 'params_data': {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1, 'num_value': len(data)}}


def update_params_fit_distribution(params_data, orders={'a': 1, 'loc': 1, 'scale': 1, 'min': 1, 'max': 1}):
    """ linear fit to model parameters as a function of a given quantity tpes_model
    and update the model dictionary with the fits of parameters """

    # update the parameters for evaluation with other values of tpes
    params = {}
    params['params_data'] = params_data  # save the data params for plotting

    if len(params_data) > 1:  # only try that if we have a sequence of more than two fits
        tpes_model = [*params_data]  # fancy python3 way to get dict.keys()
        params_values = [*params_data.values()]

        try:
            As = np.array([v['a'] for v in params_values])
            # prevent large or small values of a from bad fits
            As = np.clip(As, A_MIN, A_MAX)
        except:
            pass
        locs = np.array([v['loc'] for v in params_values])
        scales = np.array([v['scale'] for v in params_values])
        mins = np.array([v['min'] for v in params_values])
        maxs = np.array([v['max'] for v in params_values])
        w = np.array([v['num_value'] for v in params_values])

        try:
            params['a'] = build_spline(tpes_model, As, w)
        except:
            pass
        params['loc'] = build_spline(tpes_model, locs, w)
        params['scale'] = build_spline(tpes_model, scales, w)
        params['min'] = build_spline(tpes_model, mins, w)
        params['max'] = build_spline(tpes_model, maxs, w)
        params['num_value'] = np.sum(w)

    else:
        params['a'] = [0.]
        params['loc'] = [0.]
        params['scale'] = [0.]
        params['min'] = 0.
        params['max'] = 0.1
        params['num_value'] = 0

    return params
