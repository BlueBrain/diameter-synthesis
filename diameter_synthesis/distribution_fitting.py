import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom as nm
from neurom.core import iter_sections

import diameter_synthesis.utils as utils 
from diameter_synthesis.utils import get_diameters, set_diameters 

####################################
## Distribution related functions ##
####################################

ROUND = 3 #number opf digits for the fitted parameters

def evaluate_distribution(x, model, tpes = []):
    """ evaluate the fit of a distribution"""

    params = model['params']

    if model['distribution'] == 'expon_rev':
        from scipy.stats import expon

        return expon.pdf(1 - x, params['loc'], params['scale'])

    elif model['distribution'] == 'exponnorm':
        from scipy.stats import exponnorm 

        return exponnorm.pdf(x, params['a'], params['loc'], params['scale'])


    elif model['distribution'] == 'gamma':
        from scipy.stats import gamma

        return gamma.pdf(x, params['a'], params['loc'], params['scale'])

    elif model['distribution'] == 'skewnorm':
        from scipy.stats import skewnorm 

        return skewnorm.pdf(x, params['a'], params['loc'], params['scale'])

    elif model['distribution'] == 'exponnorm_sequence':
        from scipy.stats import exponnorm

        fits = []
        for tpe in tpes:
            fits.append(exponnorm.pdf(x, np.poly1d(params['a'])(tpe), loc = np.poly1d(params['loc'])(tpe), scale = np.poly1d(params['scale'])(tpe)))

        return fits
    else:
        raise Exception('Distribution not understood')

def sample_distribution(model, tpe = 0):
    """ sample from a distribution"""

    def truncate(sample_func, min_value, max_value):
        sample = sample_func()

        while sample > max_value or sample < min_value:
            sample = sample_func()

        return sample 

    params = model['params']

    if model['distribution'] == 'expon_rev':
        from scipy.stats import expon

        return truncate(lambda: 1. - expon.rvs(params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'exponnorm':
        from scipy.stats import exponnorm 

        return truncate(lambda: exponnorm.rvs(params['a'], params['loc'], params['scale']), params['min'], params['max'])


    elif model['distribution'] == 'gamma':
        from scipy.stats import gamma

        return truncate(lambda: gamma.rvs(params['a'], params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'skewnorm':
        from scipy.stats import skewnorm 

        return truncate(lambda: skewnorm.rvs(params['a'], params['loc'], params['scale']), params['min'], params['max'])

    elif model['distribution'] == 'exponnorm_sequence':
        from scipy.stats import exponnorm

        params_all =  {                                                                          
                   "a": [
                            0.086, 
                            2.921
                        ], 
                        "loc": [
                            0.0, 
                            0.0
                        ], 
                        "max": [
                            0.45, 
                            2.068
                        ], 
                        "min": [
                            0.1, 
                            0.457
                        ], 
                        "scale": [
                            0.055, 
                            0.355
                        ]
                    }
 
        if tpe == 0:
            tpe = 1


        if params['a'][0]<0:
            params['a'][0] = 1

        a = np.poly1d(params['a'])(tpe)
        loc = np.poly1d(params['loc'])(tpe)

        if params['scale'][0]<0:
            params['scale'][0] = 1

        scale = np.poly1d(params['scale'])(tpe)
        Min = np.poly1d(params['min'])(tpe)
        Max = np.poly1d(params['max'])(tpe)

        #hack to use the all data values if the fit failed
        if [*params.values()][0] == [0.,0.] or np.array([a, loc, scale, Min, Max]).any()<0:

            a = np.poly1d(params_all['a'])(tpe)
            loc = np.poly1d(params_all['loc'])(tpe)
            scale = np.poly1d(params_all['scale'])(tpe)
            Min = np.poly1d(params_all['min'])(tpe)
            Max = np.poly1d(params_all['max'])(tpe)

        try: 
            return truncate(lambda: exponnorm.rvs(a, loc, scale), Min, Max)
        except:
            print('error in parameters for tpe',tpe, 'with params', [a, loc, scale, Min, Max], 'and', params)
    else:
        raise Exception('Distribution not understood')


def fit_distribution_params(params): 
    """ linear fit to model parameters as a function of a given quantity tppes_model """

    tpes_model = [*params] #fancy python3 way to get dict.keys()
    params_values = [*params.values()]
    As     = [v['a'] for v in params_values]
    locs   = [v['loc'] for v in params_values]
    scales = [v['scale'] for v in params_values]
    mins = [v['min'] for v in params_values]
    maxs = [v['max'] for v in params_values]
    z_As = np.polyfit(tpes_model, As, 1)
    z_locs = np.polyfit(tpes_model, locs, 1)
    z_scales = np.polyfit(tpes_model, scales, 1)
    z_mins = np.polyfit(tpes_model, mins, 1)
    z_maxs = np.polyfit(tpes_model, maxs, 1)
    
    return list(np.round(z_As, ROUND)), list(np.round(z_locs, ROUND)), list(np.round(z_scales, ROUND)), list(np.round(z_mins, ROUND)), list(np.round(z_maxs,ROUND))

def fit_distribution(data, distribution, floc = None, fa = None,  min_sample_num = 10, p = 5, n_bins = 10):
    """ generic function to fit a distribution with scipy """

    if len(data) > 0:

        if distribution == 'expon_rev':
            from scipy.stats import expon
            loc, scale = expon.fit(1 - np.array(data)) 

            return {'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': 1}

        elif distribution == 'exponnorm':
            from scipy.stats import exponnorm
            a, loc, scale = exponnorm.fit(data) 

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'skewnorm':
            from scipy.stats import skewnorm 
            a, loc, scale = skewnorm.fit(data) 

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'gamma':
            from scipy.stats import gamma
            if floc is not None:
                a, loc, scale = gamma.fit(data, floc = floc) 
            else:
                a, loc, scale = gamma.fit(data)

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'exponnorm_sequence':
            from scipy.stats import exponnorm 

            tpes = np.asarray(data)[:, 1] #collect the type of point (branching order for now)
            values = np.asarray(data)[:, 0] #collect the data itself

            #set the bins for estimating parameters
            bins = utils.set_bins(tpes, n_bins, n_min = min_sample_num)
            params = {}
            for i in range(n_bins-1):
                values_tpe = values[(tpes >= bins[i]) & (tpes < bins[i+1]) ] #select the values by its type 
                if len(values_tpe) > min_sample_num: #if enough points, try to fit (intermediate bins could be almost empty)
                    if floc is not None:
                        a, loc, scale = exponnorm.fit(values_tpe, floc = floc) 
                    elif fa is not None:
                        a, loc, scale = exponnorm.fit(values_tpe, f0 = fa) 
                    else:
                        a, loc, scale = exponnorm.fit(values_tpe)

                    params[(bins[i+1] + bins[i])/2.] = {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(values_tpe, p), ROUND), 'max': np.round(np.percentile(values_tpe, 100-p), ROUND)}

            return params

        else:
            raise Exception('Distribution not understood')
    else:
        # if no data, return null parameters (for neurons without apical dentrites)
        if distribution == 'exponnorm_sequence':
            return {0. : {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1 }}
        else:
            return {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1 }

def update_params_fit_distribution(data, model):
    """ update the model dictionary with the fits of parameters """

    if len(data) > 0:
        tpes = np.asarray(data)[:, 1] #collect the type of point (branching order for now)
        values = np.asarray(data)[:, 0] #collect the data itself

        z_As, z_locs, z_scales, z_mins, z_maxs = fit_distribution_params(model['params']) #fit them as function of tpes
        #update the parameters for evaluation with other values of tpes
        model['params_data'] = model['params'] #save the data params for plotting
        model['params'] = {'a': z_As, 'loc': z_locs, 'scale': z_scales, 'min': z_mins, 'max': z_maxs}

    else:
        # if no data, return null parameters (for neurons without apical dentrites)
        model['params_data'] = model['params'] #save the data params for plotting
        model['params'] = {'a': [0.], 'loc': [0.], 'scale': [0.], 'min': 0., 'max': 0.1 }

    return model



