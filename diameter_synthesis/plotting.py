import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import diameter_synthesis.utils as utils 

colors = {'basal': 'b', 'apical': 'r'}

def plot_fit_distribution_params(model, neurite_types, fig_name = 'test', ext = '.png'):
    """ plot the fit of the parameter of the distributions """

    fig = plt.figure(figsize=(7,5))

    fig.subplots_adjust(hspace=0.0)
    ax1 = fig.add_subplot(311)

    for neurite_type in neurite_types:
        tpes_model = model[neurite_type]['params_data'].keys()
        ax1.plot(tpes_model, np.poly1d(model[neurite_type]['params']['a'])(tpes_model), c=colors[neurite_type])
        As     = [v['a'] for v in model[neurite_type]['params_data'].values()]
        ax1.plot(tpes_model, As, '+', c=colors[neurite_type])

    ax1.set_xlabel('max branching order')
    ax1.set_ylabel('a')

    ax2 = fig.add_subplot(312)

    for neurite_type in neurite_types:
        tpes_model = model[neurite_type]['params_data'].keys()
        ax2.plot(tpes_model, np.poly1d(model[neurite_type]['params']['loc'])(tpes_model), c=colors[neurite_type])
        locs   = [v['loc'] for v in model[neurite_type]['params_data'].values()]
        ax2.plot(tpes_model, locs, '+', c=colors[neurite_type])

    ax2.set_xlabel('max branching order')
    ax2.set_ylabel('loc')

    ax3 = fig.add_subplot(313)

    for neurite_type in neurite_types:
        tpes_model = model[neurite_type]['params_data'].keys()
        ax3.plot(tpes_model, np.poly1d(model[neurite_type]['params']['scale'])(tpes_model), c=colors[neurite_type])
        scales = [v['scale'] for v in model[neurite_type]['params_data'].values()]
        ax3.plot(tpes_model, scales, '+', c=colors[neurite_type])

    ax3.set_xlabel('max branching order')
    ax3.set_ylabel('scale')

    fig.savefig(fig_name + ext, bbox_inches='tight')

    plt.close(fig)


def plot_distribution_fit(data, model, neurite_types, fig_name = 'test', ext = '.png', figsize = (5,4)):
    """ Plot the data distribution and its fit """
    fig = plt.figure(figsize = figsize)
    for neurite_type in neurite_types:
        if len(data[neurite_type]) > 0:
            if isinstance(data[neurite_type][0], list): #if the fits are done as a function of a type

                tpes = np.asarray(data[neurite_type])[:, 1] #collect the type of point (branching order for now)
                values = np.asarray(data[neurite_type])[:, 0] #collect the data itself

                for tpe in set(tpes):

                    values_tpe = values[tpes==tpe]
                    bottom_shift = 0.2*tpe

                    plt.axhline(bottom_shift, lw=0.2, c='k')
                    plt.hist(values_tpe, bins = 50, log = False, density = True, histtype='bar', lw=0.5, bottom = bottom_shift, color=colors[neurite_type], alpha=0.2)

                    min_tpe =  np.poly1d(model[neurite_type]['params']['min'])(tpe)
                    max_tpe =  np.poly1d(model[neurite_type]['params']['max'])(tpe)
                    x = np.linspace(min_tpe, max_tpe)

                    #plt.axvline(min_tpe, ls='--', c='k')
                    #plt.axvline(max_tpe, ls='--', c='k')

                    plt.plot(x, bottom_shift + utils.evaluate_distribution(x, model[neurite_type], tpes = [tpe,])[0], c=colors[neurite_type], lw = 1, ls='--')#, label = neurite_type)
                plt.gca().set_ylim(0, bottom_shift + 1)

            else:

                min_val = model[neurite_type]['params']['min']
                max_val = model[neurite_type]['params']['max']

                if len(data[neurite_type]) >0:
                    plt.hist(data[neurite_type], bins = 50 , log = False, density = True, histtype='bar', lw=0.5, color=colors[neurite_type], alpha=0.5, range=[min_val*0.8, max_val*1.2])

                    x = np.linspace(min_val, max_val, 1000)
                    plt.plot(x, utils.evaluate_distribution(x, model[neurite_type]), c=colors[neurite_type], lw = 3, ls='--', label = neurite_type)

                plt.legend(loc='best')
                plt.gca().set_xlim(min_val*0.8, max_val*1.2)

    title_txt = 'Fit parameters:\n'
    for neurite_type in neurite_types:
        title_txt += neurite_type +': ' + str(model[neurite_type]['params']) 
        title_txt += '\n'

    plt.title(title_txt, fontsize = 8)

    plt.savefig(fig_name + ext, bbox_inches='tight')
    plt.close(fig)

    if isinstance(data[neurite_types[0]][0], list): #if the fits are done as a function of a type
        plot_fit_distribution_params(model, neurite_types, fig_name = fig_name+  '_param_fit', ext = ext)

def plot_fit_param_boxes(model_params, model = 'M0', neurite_type = 'basal', figname = 'test', ext = '.png', figsize = (6,3)):
    """ box plots for the fits of the model parameters """

    import collections 
    data = collections.OrderedDict() 

    mtype = model_params.keys()[0]
    for fit in  model_params[mtype][model]:
        if fit != 'trunk_diameter':
            for params in model_params[mtype][model][fit][neurite_type]['params']:
                data[fit + '_' + params] = []
        else:
            for params in model_params[mtype][model][fit][neurite_type]['params']:
                data[fit + '_' + params + '_0'] = []
                data[fit + '_' + params + '_1'] = []



    for mtype in model_params:
        for fit in  model_params[mtype][model]:
            if fit != 'trunk_diameter':
                for params in model_params[mtype][model][fit][neurite_type]['params']:
                    data[fit+'_' + params].append(model_params[mtype][model][fit][neurite_type]['params'][params])
            else:
                for params in model_params[mtype][model][fit][neurite_type]['params']:
                    data[fit+'_' + params + '_0'].append(model_params[mtype][model][fit][neurite_type]['params'][params][0])
                    data[fit+'_' + params + '_1'].append(model_params[mtype][model][fit][neurite_type]['params'][params][1])


    plt.figure(figsize = figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data)+1), data.keys(), rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + ext, bbox_inches = 'tight')

    data = collections.OrderedDict() 

    bos = ['0.0', '1.0', '2.0', '3.0', '4.0', '5.0','6.0','7.0','8.0',]
    mtype = model_params.keys()[0]
    for fit in  model_params[mtype][model]:
        if fit == 'trunk_diameter':
            for params in model_params[mtype][model][fit][neurite_type]['params_data']['0.0']:
                for bo in bos:
                    data[bo + '_' + params] = []


    print(data)
    for mtype in model_params:
        for fit in  model_params[mtype][model]:
            if fit == 'trunk_diameter':
                for bo in model_params[mtype][model][fit][neurite_type]['params_data']:
                    if bo in bos:
                        for params in model_params[mtype][model][fit][neurite_type]['params_data'][bo]:
                            data[bo + '_' + params].append(model_params[mtype][model][fit][neurite_type]['params_data'][bo][params])

    plt.figure(figsize = figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data)+1), data.keys(), rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + '_trunk' + ext, bbox_inches = 'tight')
