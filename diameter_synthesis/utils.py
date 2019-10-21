import os, glob, shutil
import json
from tqdm import tqdm

import numpy as np

import neurom 

ROUND = 3 #number opf digits for the fitted parameters

##############################
## loading/saving functions ##
##############################

def load_morphologies_from_dict(morph_path, name_dict):
    """ Load the morphologies from a list of files """

    # just to get single progression bar 
    if len(name_dict)==1:
        tqdm_1 = True
        tqdm_2 = False
    else:
        tqdm_1 = False
        tqdm_2 = True 

    morphologies = {}
    for mtype in tqdm(name_dict, disable = tqdm_1):
        morphologies[mtype] = []
        for fname in tqdm(name_dict[mtype], disable = tqdm_2):
            name, ext = os.path.splitext(fname)
            if ext in {'.h5', '.asc', '.swc'}:
                neuron = neurom.load_neuron(morph_path + '/' + fname)
                morphologies[mtype].append([neuron, name])

    return morphologies 

def load_morphologies(morph_path, by_mtypes = True, n_morphs_max = None, n_mtypes_max = None, xml_file = './neuronDB.xml', ext = '.asc'):
    """ Load the morphologies from a directory, by mtypes or all at once """
    
    #load morphologies by mtypes, one listt of morphologies for each type
    if by_mtypes:

        #if not max number of mtypes, take all
        if not n_mtypes_max:
            n_mtypes_max =  1e10


        from xml.etree import ElementTree as ET
        FileDB = ET.parse(morph_path + xml_file)
        root = FileDB.findall('listing')[0]
        morphs = root.findall('morphology')
    
        name_dict = {}
        for m in morphs:
            try:
                # Define mtypes
                mtype = m.find('mtype').text
                # Define suntypes (if they exist)
                if m.find('msubtype').text:
                    mtype = mtype + ':' + m.find('msubtype').text

                #if it is a new mtype, add en entry to name_dict
                if mtype not in name_dict.keys() and len(name_dict) < n_mtypes_max:
                    name_dict[mtype] = [m.find('name').text + ext]
                elif mtype in name_dict.keys():
                    name_dict[mtype] += [m.find('name').text + ext]

            except:
                print('Failed to process', m)

    #load all the morphologies together
    else:
        name_dict = {}
        if n_morphs_max is not None:
            name_dict['all_types'] = os.listdir(morph_path)[:n_morphs_max]
        else:
            name_dict['all_types'] = os.listdir(morph_path)
    
    return load_morphologies_from_dict(morph_path, name_dict)
    

############################
## morphometric functions ##
############################

def sibling_ratios(neurite, method = 'mean'):
    """ compute the siblig ratios of a neurite"""

    s_ratios = []
    #loop over bifuraction points
    for bif_point in neurom.core.Tree.ibifurcation_point(neurite.iter_sections().next()):
        if len(bif_point.children) == 2:

            if method == 'mean':
                #take the average diameter of children to smooth noise out
                d1 = np.mean([2*p[3] for p in bif_point.children[0].points])
                d2 = np.mean([2*p[3] for p in bif_point.children[1].points])

            elif method == 'first':
                #take the first diameter, but subject to noise!
                d1 = 2*bif_point.children[0].points[0,3]
                d2 = 2*bif_point.children[1].points[0,3]

            else:
                raise Exception('Method for singling computation not understood!')

            s_ratios.append( np.min([d1,d2]) / np.max([d1,d2]) ) 

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return s_ratios

def Rall_deviations(neurite, method = 'mean'):
    """Returns the Rall deviation the diameters
       of the segments of a tree. """
    
    Rall_deviations = []
    for bif_point in neurom.core.Tree.ibifurcation_point(neurite.iter_sections().next()):
        if len(bif_point.children) == 2:

            if method == 'mean':

                d_0 = np.mean([2*p[3] for p in bif_point.points])
                d_1 = np.mean([2*p[3] for p in bif_point.children[0].points])
                d_2 = np.mean([2*p[3] for p in bif_point.children[1].points])

            elif method == 'first':

                d_0 = 2*bif_point.points[-1, 3]
                d_1 = 2*bif_point.children[0].points[0, 3]
                d_2 = 2*bif_point.children[1].points[0, 3]

            Rall_deviations.append( (d_1/d_0)**(3./2.) + (d_2/d_0)**(3./2.) )

        elif len(bif_point.children) > 2:
            raise Exception('Number of children is '+ str(len(bif_point.children)) + '!')

    return Rall_deviations

def terminal_diameters(neurite, method = 'mean', threshold = 0.8):
    """Returns the model for the terminations"""

    mean_radii = np.mean(neurite.points[:, 3])

    if method == 'mean':
        term_diam = [2. * np.mean(t.points[:, 3]) for t in neurom.core.Tree.ileaf(next(neurite.iter_sections())) if np.mean(t.points[:, 3]) < threshold * mean_radii]

    elif method == 'first':
        term_diam = [2. * t.points[-1, 3] for t in neurom.core.Tree.ileaf(next(neurite.iter_sections())) if t.points[-1, 3] < threshold * mean_radii]

    else:
        raise Exception('Method for singling computation not understood!')

    return term_diam

def trunk_diameter(neurite):

    trunk_diam =  2*neurite.root_node.points[0][3] 
    max_bo = np.max(neurom.get('section_term_branch_orders', neurite))

    return [[trunk_diam, max_bo], ]

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
            fits.append(exponnorm.pdf(x, np.poly1d(params['a'])(tpe), np.poly1d(params['loc'])(tpe), np.poly1d(params['scale'])(tpe)))

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

        return truncate(lambda: exponnorm.rvs(np.poly1d(params['a'])(tpe), np.poly1d(params['loc'])(tpe), np.poly1d(params['scale'])(tpe)), np.poly1d(params['min'])(tpe), np.poly1d(params['max'])(tpe))

        return fits
    else:
        raise Exception('Distribution not understood')


def fit_distribution_params(params): 
    tpes_model = params.keys()
    As     = [v['a'] for v in params.values()]
    locs   = [v['loc'] for v in params.values()]
    scales = [v['scale'] for v in params.values()]
    mins = [v['min'] for v in params.values()]
    maxs = [v['max'] for v in params.values()]

    print(As)
    z_As = np.polyfit(tpes_model, As, 1)
    z_locs = np.polyfit(tpes_model, locs, 1)
    z_scales = np.polyfit(tpes_model, scales, 1)
    z_mins = np.polyfit(tpes_model, mins, 1)
    z_maxs = np.polyfit(tpes_model, maxs, 1)
    
    return list(np.round(z_As, ROUND)), list(np.round(z_locs, ROUND)), list(np.round(z_scales, ROUND)), list(np.round(z_mins, ROUND)), list(np.round(z_maxs,ROUND))

def fit_distribution(data, distribution, floc = None, min_sample_num = 10, p = 5):
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

            params = {}
            for tpe in set(tpes):
                values_tpe = values[tpes==tpe] #select the values by its type 
                if len(values_tpe) > min_sample_num: #if enough points, try to fit
                    if floc is not None:
                        a, loc, scale = exponnorm.fit(values_tpe, floc = floc) 
                    else:
                        a, loc, scale = exponnorm.fit(values_tpe)

                    params[tpe] = {'a': np.round(a, ROUND), 'loc': np.round(loc, ROUND), 'scale': np.round(scale, ROUND), 'min': np.round(np.percentile(values_tpe, p), ROUND), 'max': np.round(np.percentile(values_tpe, 100-p), ROUND)}

            if len(params) < 2:
                #print('Not enough datapoints to fit the distribution with ', len(params), 'points.')
                params[0.] = {'a': 0., 'loc': 0., 'scale': 0., 'min': 0, 'max': 0}
                params[1.] = {'a': 0., 'loc': 0., 'scale': 0., 'min': 0, 'max': 0}

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


def save_neuron(neuron, model, folder):
        """ save the neuron morphology for later analysis """

        if not os.path.isdir(folder):
                os.mkdir(folder)

        neuron[0].write(folder + '/' + model + '_' + neuron[1] + '.asc')


