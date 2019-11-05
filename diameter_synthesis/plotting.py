import os, glob, shutil
import numpy as np

import matplotlib
matplotlib.use('Agg')
import pylab as plt

from neurom.core import iter_sections
from neurom import viewer

from diameter_synthesis.distribution_fitting import evaluate_distribution
import diameter_synthesis.utils as utils 

########################
## plotting functions ##
########################

colors = {'basal': 'b', 'apical': 'r'}

def plot_fit_distribution_params(model, neurite_types, fig_name = 'test', ext = '.png'):
    """ plot the fit of the parameter of the distributions """

    fig = plt.figure(figsize=(5,8))

    fig.subplots_adjust(hspace=0.0)
    ax1 = fig.add_subplot(511)

    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params_data']]
        ax1.plot(tpes_model, np.poly1d(model[neurite_type]['params']['a'])(tpes_model), c=colors[neurite_type])
        As     = [v['a'] for v in model[neurite_type]['params_data'].values()]
        ax1.plot(tpes_model, As, '+', c=colors[neurite_type])

    ax1.set_xlabel('max branching order')
    ax1.set_ylabel('a')

    ax2 = fig.add_subplot(512)

    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params_data']]
        ax2.plot(tpes_model, np.poly1d(model[neurite_type]['params']['loc'])(tpes_model), c=colors[neurite_type])
        locs   = [v['loc'] for v in model[neurite_type]['params_data'].values()]
        ax2.plot(tpes_model, locs, '+', c=colors[neurite_type])

    ax2.set_ylabel('loc')

    ax3 = fig.add_subplot(513)

    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params_data']]
        ax3.plot(tpes_model, np.poly1d(model[neurite_type]['params']['scale'])(tpes_model), c=colors[neurite_type])
        scales = [v['scale'] for v in model[neurite_type]['params_data'].values()]
        ax3.plot(tpes_model, scales, '+', c=colors[neurite_type])

    ax3.set_ylabel('scale')

    ax4 = fig.add_subplot(514)

    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params_data']]
        ax4.plot(tpes_model, np.poly1d(model[neurite_type]['params']['min'])(tpes_model), c=colors[neurite_type])
        mins = [v['min'] for v in model[neurite_type]['params_data'].values()]
        ax4.plot(tpes_model, mins, '+', c=colors[neurite_type])

    ax4.set_ylabel('min')

    ax5 = fig.add_subplot(515)

    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params_data']]
        ax5.plot(tpes_model, np.poly1d(model[neurite_type]['params']['max'])(tpes_model), c=colors[neurite_type])
        maxs = [v['max'] for v in model[neurite_type]['params_data'].values()]
        ax5.plot(tpes_model, maxs, '+', c=colors[neurite_type])

    ax5.set_xlabel('max branching order')
    ax5.set_ylabel('max')



    fig.savefig(fig_name + ext, bbox_inches='tight')
    plt.close(fig)


def plot_distribution_fit(data, model, neurite_types, fig_name = 'test', ext = '.png', figsize = (5,4), n_bins = 10):
    """ Plot the data distribution and its fit """
    fig = plt.figure(figsize = figsize)
    for neurite_type in neurite_types:
        if len(data[neurite_type]) > 0:
            if isinstance(data[neurite_type][0], list): #if the fits are done as a function of a type

                tpes = np.asarray(data[neurite_type])[:, 1] #collect the type of point (branching order for now)
                values = np.asarray(data[neurite_type])[:, 0] #collect the data itself
               
                bins = utils.set_bins(tpes, n_bins)
                min_val = 1e10
                max_val = -1e10

                for i in range(n_bins-1):
                    values_tpe = values[(tpes >= bins[i]) & (tpes < bins[i+1]) ] #select the values by its type 
                    tpe_mean = (bins[i]+bins[i+1])/2.
             
                    bottom_shift = tpe_mean

                    plt.axhline(bottom_shift, lw=0.2, c='k')

                    min_tpe = np.poly1d(model[neurite_type]['params']['min'])(tpe_mean)
                    max_tpe = np.poly1d(model[neurite_type]['params']['max'])(tpe_mean)

                    min_val = np.min([min_val, min_tpe])
                    max_val = np.max([max_val, max_tpe])
                    values_tpe = values_tpe[values_tpe<max_tpe*1.2] 
                    values_tpe = values_tpe[values_tpe>min_tpe*0.8] 

                    n, b = np.histogram(values_tpe, bins = 20)#, range=(min_tpe*0.8, max_tpe*1.2))
                    plt.bar(b[:-1], np.array(n)/np.max(n), width = b[1]-b[0], bottom = bottom_shift, color=colors[neurite_type], alpha = 0.5)
                    #plt.hist(values_tpe, bins = 50, log = False, density = True, histtype='bar', lw=0.5, bottom = bottom_shift, color=colors[neurite_type], alpha=0.2, range=(min_tpe*0.8, max_tpe*1.2))

                    x = np.linspace(min_tpe, max_tpe, 1000)
                    pdf = evaluate_distribution(x, model[neurite_type], tpes = [tpe_mean,])[0]
                    plt.plot(x, bottom_shift + pdf/np.max(pdf) , c=colors[neurite_type], lw = 1, ls='--')#, label = neurite_type)

                #plt.gca().set_ylim(0, bottom_shift + 1)
                #plt.gca().set_xlim(min_val*0.8, max_val*1.2)

            else:

                min_val = model[neurite_type]['params']['min']
                max_val = model[neurite_type]['params']['max']

                if len(data[neurite_type]) >0:
                    plt.hist(data[neurite_type], bins = 50 , log = False, density = True, histtype='bar', lw=0.5, color=colors[neurite_type], alpha=0.5, range=[min_val*0.8, max_val*1.2])

                    x = np.linspace(min_val, max_val, 1000)
                    plt.plot(x, evaluate_distribution(x, model[neurite_type]), c=colors[neurite_type], lw = 3, ls='--', label = neurite_type)

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

    mtype = [*model_params][0]
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
    plt.xticks(np.arange(1, len(data)+1), [*data], rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + ext, bbox_inches = 'tight')
    plt.close()

    data = collections.OrderedDict() 

    bos = ['0.0', '1.0', '2.0', '3.0', '4.0', '5.0','6.0','7.0','8.0',]
    mtype = [*model_params][0]
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
    plt.xticks(np.arange(1, len(data)+1), [*data], rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + '_trunk' + ext, bbox_inches = 'tight')
    plt.close()

def plot_neuron(neuron, folder, ext = '.png'):
    """ plot a neuron and save in a folder """

    if not os.path.isdir(folder):
            os.mkdir(folder)

    fig, ax = viewer.draw(neuron[0])
    plt.savefig(folder + '/' + neuron[1] + '_' + folder + ext, dpi = 500)
    plt.close()




from neurom.view import (plot_neuron, plot_neuron3d,
                   plot_tree, plot_tree3d,
                   plot_soma, plot_soma3d,
                   plot_dendrogram)

from neurom.view import common
from neurom.core import Tree, Neurite, Soma, Neuron

MODES = ('2d', '3d', 'dendrogram')

_VIEWERS = {
    'neuron_3d': plot_neuron3d,
    'neuron_2d': plot_neuron,
    'neuron_dendrogram': plot_dendrogram,
    'tree_3d': plot_tree3d,
    'tree_2d': plot_tree,
    'tree_dendrogram': plot_dendrogram,
    'soma_3d': plot_soma3d,
    'soma_2d': plot_soma
}

class ViewerError(Exception):
    '''Base class for viewer exceptions'''


class InvalidDrawModeError(ViewerError):
    '''Exception class to indicate invalid draw mode'''

class NotDrawableError(Exception):
    '''Exception class for things that aren't drawable'''

def draw_axis(obj, mode='2d', ax = None, **kwargs):
    '''Draw a morphology object

    Parameters:
        obj: morphology object to be drawn (neuron, tree, soma).
        mode (Optional[str]): drawing mode ('2d', '3d', 'dendrogram'). Defaults to '2d'.
        **kwargs: keyword arguments for underlying neurom.view.view functions.

    Raises:
        InvalidDrawModeError if mode is not valid
        NotDrawableError if obj is not drawable
        NotDrawableError if obj type and mode combination is not drawable

    Examples:

        >>> nrn = ... # load a neuron
        >>> fig, _ = viewer.draw(nrn)             # 2d plot
        >>> fig.show()
        >>> fig3d, _ = viewer.draw(nrn, mode='3d') # 3d plot
        >>> fig3d.show()
        >>> fig, _ = viewer.draw(nrn.neurites[0]) # 2d plot of neurite tree
        >>> dend, _ = viewer.draw(nrn, mode='dendrogram')

    '''

    if mode not in MODES:
        raise InvalidDrawModeError('Invalid drawing mode %s' % mode)

    if mode in ('2d', 'dendrogram'):
        if not ax:
            fig, ax = common.get_figure()
        else: 
            fig, ax_1 = common.get_figure()
        
    else:
        fig, ax = common.get_figure(params={'projection': '3d'})

    if isinstance(obj, Neuron):
        tag = 'neuron'
    elif isinstance(obj, (Tree, Neurite)):
        tag = 'tree'
    elif isinstance(obj, Soma):
        tag = 'soma'
    else:
        raise NotDrawableError('draw not implemented for %s' % obj.__class__)

    viewer = '%s_%s' % (tag, mode)
    try:
        plotter = _VIEWERS[viewer]
    except KeyError:
        raise NotDrawableError('No drawer for class %s, mode=%s' % (obj.__class__, mode))

    output_path = kwargs.pop('output_path', None)
    plotter(ax, obj, **kwargs)

    if mode != 'dendrogram':
        common.plot_style(fig=fig, ax=ax, **kwargs)

    if output_path:
        common.save_plot(fig=fig, output_path=output_path, **kwargs)

    return fig, ax


def plot_diameter_diff(neuron_name, morph_path, new_morph_path, model, neurite_types, folder):
    """ plot original morphology, new one and differences """    

    if not os.path.isdir(folder):
        os.mkdir(folder)

    neuron_orig     = utils.load_neuron(neuron_name, None, morph_path)
    neuron_diff_pos = utils.load_neuron(neuron_name, None, morph_path)
    neuron_diff_neg = utils.load_neuron(neuron_name, None, morph_path)
    neuron_new      = utils.load_neuron(neuron_name, model, new_morph_path)

    fig, axs = plt.subplots(2, 2, figsize=(10,10))

    draw_axis(neuron_orig, ax = axs[0,0])
    axs[0,0].set_title('Original neuron')

    draw_axis(neuron_new, ax= axs[0,1])
    axs[0,1].set_title('New neuron')

    neurite_types = ['basal', 'axon', 'apical']
    for neurite_type in neurite_types:
        neurites = [neurite for neurite in neuron_orig.neurites if neurite.type == utils.STR_TO_TYPES[neurite_type]]
        neurites_new = [neurite for neurite in neuron_new.neurites if neurite.type == utils.STR_TO_TYPES[neurite_type]]
        neurites_diff_neg = [neurite for neurite in neuron_diff_neg.neurites if neurite.type == utils.STR_TO_TYPES[neurite_type]]
        neurites_diff_pos = [neurite for neurite in neuron_diff_pos.neurites if neurite.type == utils.STR_TO_TYPES[neurite_type]]

        for i, neurite in enumerate(neurites):

            diam_orig = []
            for s in iter_sections(neurite):
                diam_orig.append(utils.get_diameters(s))

            diam_new = []
            for s in iter_sections(neurites_new[i]):
                diam_new.append(utils.get_diameters(s))

            for j, s in enumerate(iter_sections(neurites_diff_pos[i])):
                diff = diam_new[j] - diam_orig[j]
                diff_pos = diff.copy()
                diff_pos[diff_pos<0] = 0 
                utils.set_diameters(s, diff_pos)

            for j, s in enumerate(iter_sections(neurites_diff_neg[i])):
                diff = diam_new[j] - diam_orig[j]
                diff_neg = -diff.copy()
                diff_neg[diff_neg<0] = 0 
                utils.set_diameters(s, diff_neg)


    draw_axis(neuron_diff_pos, ax= axs[1,0])
    axs[1,0].set_title('Positive diameter differences')
    draw_axis(neuron_diff_neg, ax= axs[1,1])
    axs[1,1].set_title('Negative diameter differences')

    fig.savefig(folder + '/' + neuron_name + '_' + folder+'_'+ model + '.png', dpi = 500)
    plt.close('all')


