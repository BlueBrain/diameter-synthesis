import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np
from numpy.polynomial import polynomial as polynomial

import neurom as nm
from neurom.core import iter_sections

import diameter_synthesis.utils as utils 
from diameter_synthesis.utils import get_diameters, set_diameters, ROUND, MIN_DATA_POINTS, A_MAX

####################################
## Distribution related functions ##
####################################

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
            fits.append(exponnorm.pdf(x, polynomial.polyval(tpe, params['a']), loc = polynomial.polyval(tpe, params['loc']), scale = polynomial.polyval(tpe, params['scale'])))

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
        
       
        if tpe == 0:
            tpe = 1

        if params['a'][0]<0:
            params['a'][0] = 1

        a = polynomial.polyval(tpe, params['a'])
        loc = polynomial.polyval(tpe, params['loc'])

        if params['scale'][0]<0:
            params['scale'][0] = 1

        scale = polynomial.polyval(tpe, params['scale'])
        Min = polynomial.polyval(tpe, params['min'])
        Max = polynomial.polyval(tpe, params['max'])

        #hack to use the all data values if the fit failed (only basal)
        if [*params.values()][0][0] == 0. or a < 0 or loc <0  or scale <0 or Min < 0 or Max<0:
            try:
                with open('./model_params_all.json', 'r') as f:
                    params_all = json.load(f)
                params_all = params_all['all_types']['M0']['trunk_diameter']['basal']['params']
            except Exception as e:
                print(e)
                raise Exception('Could not load params_all.json, please create it or have more cells in each class.')

            a = polynomial.polyval(tpe, params_all['a'])
            loc = polynomial.polyval(tpe, params_all['loc'])
            scale = polynomial.polyval(tpe, params_all['scale'])
            Min = polynomial.polyval(tpe, params_all['min'])
            Max = polynomial.polyval(tpe, params_all['max'])

        try: 
            return truncate(lambda: exponnorm.rvs(a, loc, scale), Min, Max)
        except:
            print('error in parameters for tpe',tpe, 'with params', [a, loc, scale, Min, Max], 'and', params)
    else:
        raise Exception('Distribution not understood')

def fit_distribution(data, distribution, floc = None, fa = None,  min_sample_num = 10, p = 5, n_bins = 10):
    """ generic function to fit a distribution with scipy """
    if len(data) > MIN_DATA_POINTS:

        if distribution == 'expon_rev':
            from scipy.stats import expon
            loc, scale = expon.fit(1 - np.array(data)) 

            return {'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': 1}

        elif distribution == 'exponnorm':
            from scipy.stats import exponnorm
            a, loc, scale = exponnorm.fit(data) 

            #refit if we get crazy values for a
            if a > utils.A_MAX:
                a, loc, scale = exponnorm.fit(data, f0 = A_MAX) 

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'skewnorm':
            from scipy.stats import skewnorm 
            a, loc, scale = skewnorm.fit(data) 

            #refit if we get crazy values for a
            if a > utils.A_MAX:
                a, loc, scale = skewnorm.fit(data, f0 = A_MAX) 

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'gamma':
            from scipy.stats import gamma
            if floc is not None:
                a, loc, scale = gamma.fit(data, floc = floc) 
                #refit if we get crazy values for a
                if a > utils.A_MAX:
                    a, loc, scale = gamma.fit(data, floc = floc, f0 = A_MAX) 
            else:
                a, loc, scale = gamma.fit(data)
                #refit if we get crazy values for a
                if a > utils.A_MAX:
                    a, loc, scale = gamma.fit(data, f0 = A_MAX) 

            return {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(data, p), ROUND), 'max': np.round(np.percentile(data, 100-p), ROUND)}

        elif distribution == 'exponnorm_sequence':
            from scipy.stats import exponnorm 

            tpes = np.asarray(data)[:, 1] #collect the type of point (branching order for now)
            values = np.asarray(data)[:, 0] #collect the data itself

            #set the bins for estimating parameters if we can otherwise use two bins to be able to fit later 
            bins, num_values = utils.set_bins(tpes, n_bins, n_min = min_sample_num)

            params = {}
            for i in range(len(bins)-1):
                values_tpe = values[(tpes >= bins[i]) & (tpes < bins[i+1]) ] #select the values by its type 
                if len(values_tpe)>0:
                    if floc is not None:
                        a, loc, scale = exponnorm.fit(values_tpe, floc = floc) 
                        #refit if we get crazy values for a
                        if a > utils.A_MAX:
                            a, loc, scale = exponnorm.fit(values_tpe, floc = floc, f0 = A_MAX) 
                    elif fa is not None:
                        a, loc, scale = exponnorm.fit(values_tpe, f0 = fa) 

                    else:
                        a, loc, scale = exponnorm.fit(values_tpe)
                        #refit if we get crazy values for a
                        if a > utils.A_MAX:
                            a, loc, scale = exponnorm.fit(values_tpe, f0 = A_MAX) 

                    params[np.round((bins[i+1] + bins[i])/2., ROUND)] = {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(values_tpe, p), ROUND), 'max': np.round(np.percentile(values_tpe, 100-p), ROUND), 'num_value': num_values[i]}

                else:
                    print("WARNING: could not fit anything, because of a lack of data points!")

            return params

        else:
            raise Exception('Distribution not understood')
    else:
        # if no data, return null parameters (for neurons without apical dentrites)
        if distribution == 'exponnorm_sequence':
            return {0. : {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1 , 'num_value': 0}}
        else:
            return {'a': 0., 'loc': 0., 'scale': 0., 'min': 0., 'max': 0.1 }

def update_params_fit_distribution(data, model, orders = {'a':1, 'loc':1, 'scale':1, 'min':1, 'max':1}): 
    """ linear fit to model parameters as a function of a given quantity tpes_model
    and update the model dictionary with the fits of parameters """

    #update the parameters for evaluation with other values of tpes
    model['params_data'] = model['params'] #save the data params for plotting

    if len(model['params'])>1: #only try that if we had enough samples for the first fits

        tpes_model = [*model['params']] #fancy python3 way to get dict.keys()
        params_values = [*model['params'].values()]
        As     = np.array([v['a'] for v in params_values])
        locs   = np.array([v['loc'] for v in params_values])
        scales = np.array([v['scale'] for v in params_values])
        mins = np.array([v['min'] for v in params_values])
        maxs = np.array([v['max'] for v in params_values])
        w = np.array([v['num_value'] for v in params_values])

        #prevent large values of a from bad fits
        As[As>A_MAX] = A_MAX
    
        z_As = polynomial.polyfit(tpes_model, As, orders['a'], w = w)
        z_locs = polynomial.polyfit(tpes_model, locs, orders['loc'], w = w)
        z_scales = polynomial.polyfit(tpes_model, scales, orders['scale'], w = w)
        z_mins = polynomial.polyfit(tpes_model, mins, orders['min'], w = w)
        z_maxs = polynomial.polyfit(tpes_model, maxs, orders['max'], w = w)
        
        model['params'] = {'a': list(np.round(z_As, ROUND)), 
                        'loc': list(np.round(z_locs, ROUND)),
                        'scale': list(np.round(z_scales, ROUND)),
                        'min': list(np.round(z_mins, ROUND)),
                        'max': list(np.round(z_maxs,ROUND))
                        }

    else:
        model['params'] = {'a': [0.], 'loc': [0.], 'scale': [0.], 'min': 0., 'max': 0.1 }

    return model



