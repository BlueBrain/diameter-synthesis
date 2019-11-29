from neurom.core import Tree, Neurite, Soma, Neuron
from neurom.view import common
from neurom.view import (plot_neuron, plot_neuron3d,
                         plot_tree, plot_tree3d,
                         plot_soma, plot_soma3d,
                         plot_dendrogram)
import os
import glob
import shutil
import numpy as np
from scipy import stats
from numpy.polynomial import polynomial as polynomial
from scipy.interpolate import UnivariateSpline

import matplotlib
# matplotlib.use('Agg')
import pylab as plt

import neurom
from neurom.core import iter_sections
from neurom import viewer

from diameter_synthesis.types import STR_TO_TYPES
from diameter_synthesis.distribution_fitting import evaluate_distribution, evaluate_spline
import diameter_synthesis.utils as utils
from diameter_synthesis import io

########################
## plotting functions ##
########################

colors = {'basal': 'r', 'apical': 'm', 'axon': 'b'}


def plot_fit_distribution_params(model, neurite_types, fig_name='test', ext='.png'):
    """ plot the fit of the parameter of the distributions """

    fig = plt.figure(figsize=(5, 8))

    fig.subplots_adjust(hspace=0.0)
    try:
        ax1 = fig.add_subplot(511)
        for neurite_type in neurite_types:
            tpes_model = [*model[neurite_type]['params']['params_data']]

            x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
            ax1.plot(x, evaluate_spline(x, model[neurite_type]
                                        ['params']['a']), c=colors[neurite_type])

            As = np.array([v['a'] for v in model[neurite_type]['params']['params_data'].values()])
            w = np.array([v['num_value']
                          for v in model[neurite_type]['params']['params_data'].values()])

            # prevent large values of a from bad fitting
            As[As > utils.A_MAX] = utils.A_MAX
            ax1.scatter(tpes_model, As, s=w, edgecolors=colors[neurite_type], facecolors='none')

        ax1.set_xlabel('max path distance')
        ax1.set_ylabel('a')
        n_plot = 0
    except:
        n_plot = 101
        pass

    ax2 = fig.add_subplot(512 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params']['params_data']]

        x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax2.plot(x, evaluate_spline(x, model[neurite_type]
                                    ['params']['loc']), c=colors[neurite_type])

        locs = [v['loc'] for v in model[neurite_type]['params']['params_data'].values()]
        w = np.array([v['num_value']
                      for v in model[neurite_type]['params']['params_data'].values()])

        ax2.scatter(tpes_model, locs, s=w, edgecolors=colors[neurite_type], facecolors='none')
    ax2.set_ylabel('loc')

    ax3 = fig.add_subplot(513 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params']['params_data']]

        x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax3.plot(x, evaluate_spline(x, model[neurite_type]
                                    ['params']['scale']), c=colors[neurite_type])

        scales = [v['scale'] for v in model[neurite_type]['params']['params_data'].values()]
        ax3.scatter(tpes_model, scales, s=w, edgecolors=colors[neurite_type], facecolors='none')
    ax3.set_ylabel('scale')

    ax4 = fig.add_subplot(514 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params']['params_data']]

        x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax4.plot(x, evaluate_spline(x, model[neurite_type]
                                    ['params']['min']), c=colors[neurite_type])

        mins = [v['min'] for v in model[neurite_type]['params']['params_data'].values()]
        ax4.scatter(tpes_model, mins, s=w, edgecolors=colors[neurite_type], facecolors='none')
    ax4.set_ylabel('min')

    ax5 = fig.add_subplot(515 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = [*model[neurite_type]['params']['params_data']]

        x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax5.plot(x, evaluate_spline(x, model[neurite_type]
                                    ['params']['max']), c=colors[neurite_type])

        maxs = [v['max'] for v in model[neurite_type]['params']['params_data'].values()]
        ax5.scatter(tpes_model, maxs, s=w, edgecolors=colors[neurite_type], facecolors='none')

    ax5.set_xlabel('max branching order')
    ax5.set_ylabel('max')

    fig.savefig(fig_name + ext, bbox_inches='tight')
    plt.close(fig)


def plot_distribution_fit(data, model, neurite_types, fig_name='test', ext='.png', figsize=(5, 4), n_bins=10):
    """ Plot the data distribution and its fit """
    fig = plt.figure(figsize=figsize)

    for neurite_type in neurite_types:
        if isinstance(model[neurite_type]['sequential'], str):  # if the fits are done sequantially
            tpes = np.asarray(data[neurite_type])[:, 1]  # collect the type of point
            values = np.asarray(data[neurite_type])[:, 0]  # collect the data itself
            plt.figure()
            plt.hist(tpes, bins=500, log=True)
            plt.axvline(0.1, c='k')
            plt.axvline(0.3, c='r')
            plt.savefig('test' + ext)
            plt.figure()
            plt.scatter(values, tpes, s=5)
            plt.savefig('test2' + ext)
            bins, num_values = utils.set_bins(tpes, n_bins)

            min_val = 1e10
            max_val = -1e10

            for i in range(len(bins) - 1):
                values_tpe = values[(tpes >= bins[i]) & (tpes < bins[i + 1])
                                    ]  # select the values by its type
                tpe_mean = (bins[i] + bins[i + 1]) / 2.

                bottom_shift = tpe_mean

                height = (bins[i + 1] - bins[i]) / 2.

                plt.axhline(bottom_shift, lw=0.2, c='k')
                try:  # try to plot, may not be a fit to plot
                    min_tpe = evaluate_spline(tpe_mean, model[neurite_type]['params']['min'])
                    max_tpe = evaluate_spline(tpe_mean, model[neurite_type]['params']['max'])

                    min_val = np.min([min_val, min_tpe])
                    max_val = np.max([max_val, max_tpe])
                    values_tpe = values_tpe[values_tpe < max_tpe * 1.2]
                    values_tpe = values_tpe[values_tpe > min_tpe * 0.8]

                    n, b = np.histogram(values_tpe, bins=20)
                    plt.bar(b[1:], height * np.array(n) / np.max(n), width=b[1] - b[0],
                            bottom=bottom_shift, color=colors[neurite_type], alpha=0.5)

                    x = np.linspace(min_tpe, max_tpe, 1000)
                    params = {}
                    try:
                        params['a'] = evaluate_spline(tpe_mean, model[neurite_type]['params']['a'])
                    except:
                        pass
                    params['loc'] = evaluate_spline(tpe_mean, model[neurite_type]['params']['loc'])
                    params['scale'] = evaluate_spline(
                        tpe_mean, model[neurite_type]['params']['scale'])
                    pdf = evaluate_distribution(x, model[neurite_type]['distribution'], params)

                    plt.plot(x, bottom_shift + height * pdf / np.max(pdf),
                             c=colors[neurite_type], lw=1, ls='--')
                except Exception as e:  # if not fit, just plot the histrogams
                    print('plotting exception:', e)
                    n, b = np.histogram(values_tpe, bins=20)
                    plt.bar(b[:-1], height * np.array(n) / np.max(n), width=b[1] - b[0],
                            bottom=bottom_shift, color=colors[neurite_type], alpha=0.5)

        else:

            min_val = model[neurite_type]['params']['min']
            max_val = model[neurite_type]['params']['max']

            if len(data[neurite_type]) > 0:
                if len(np.shape(data[neurite_type])) > 1:
                    plt.hist(np.array(data[neurite_type])[:, 0], bins=50, log=False, density=True, histtype='bar',
                             lw=0.5, color=colors[neurite_type], alpha=0.5, range=[min_val * 0.8, max_val * 1.2])
                else:
                    plt.hist(data[neurite_type], bins=50, log=False, density=True, histtype='bar', lw=0.5,
                             color=colors[neurite_type], alpha=0.5, range=[min_val * 0.8, max_val * 1.2])

                x = np.linspace(min_val, max_val, 1000)
                plt.plot(x, evaluate_distribution(
                    x, model[neurite_type]['distribution'], model[neurite_type]['params']), c=colors[neurite_type], lw=3, ls='--', label=neurite_type)

            plt.legend(loc='best')
            plt.gca().set_xlim(min_val * 0.8, max_val * 1.2)

    title_txt = 'Fit parameters:\n'
    try:
        for neurite_type in neurite_types:
            params_title = {}
            for p in model[neurite_type]['params']:
                if p != 'params_data':
                    params_title[p] = model[neurite_type]['params'][p]
            title_txt += neurite_type + ': ' + str(params_title)
            title_txt += '\n'
    except:
        title_txt += ' no parameter could be fitted.'

    #plt.title(title_txt, fontsize = 8)

    plt.savefig(fig_name + ext, bbox_inches='tight')
    plt.close(fig)

    if isinstance(model[neurite_type]['sequential'], str):  # if the fits are done sequantially
        plot_fit_distribution_params(
            model, neurite_types, fig_name=fig_name + '_param_fit', ext=ext)


def plot_fit_param_boxes(model_params, model='M0', neurite_type='basal', figname='test', ext='.png', figsize=(6, 3)):
    """ box plots for the fits of the model parameters """

    import collections
    data = collections.OrderedDict()

    mtype = [*model_params][0]
    for fit in model_params[mtype][model]:
        if fit != 'trunk_diameter':
            for params in model_params[mtype][model][fit][neurite_type]['params']:
                data[fit + '_' + params] = []
        else:
            for params in model_params[mtype][model][fit][neurite_type]['params']:
                data[fit + '_' + params + '_0'] = []
                data[fit + '_' + params + '_1'] = []

    for mtype in model_params:
        for fit in model_params[mtype][model]:
            if fit != 'trunk_diameter':
                for params in model_params[mtype][model][fit][neurite_type]['params']:
                    data[fit + '_' + params].append(model_params[mtype]
                                                    [model][fit][neurite_type]['params'][params])
            else:
                for params in model_params[mtype][model][fit][neurite_type]['params']:
                    data[fit + '_' + params +
                         '_0'].append(model_params[mtype][model][fit][neurite_type]['params'][params][0])
                    data[fit + '_' + params +
                         '_1'].append(model_params[mtype][model][fit][neurite_type]['params'][params][1])

    plt.figure(figsize=figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data) + 1), [*data], rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + ext, bbox_inches='tight')
    plt.close()

    data = collections.OrderedDict()

    bos = ['0.0', '1.0', '2.0', '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', ]
    mtype = [*model_params][0]
    for fit in model_params[mtype][model]:
        if fit == 'trunk_diameter':
            for params in model_params[mtype][model][fit][neurite_type]['params_data']['0.0']:
                for bo in bos:
                    data[bo + '_' + params] = []

    for mtype in model_params:
        for fit in model_params[mtype][model]:
            if fit == 'trunk_diameter':
                for bo in model_params[mtype][model][fit][neurite_type]['params_data']:
                    if bo in bos:
                        for params in model_params[mtype][model][fit][neurite_type]['params_data'][bo]:
                            data[bo + '_' + params].append(model_params[mtype][model]
                                                           [fit][neurite_type]['params_data'][bo][params])

    plt.figure(figsize=figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data) + 1), [*data], rotation='vertical')
    #plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + '_trunk' + ext, bbox_inches='tight')
    plt.close()


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


def draw_axis(obj, mode='2d', ax=None, **kwargs):
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

    neuron_orig = io.load_neuron(neuron_name, None, morph_path)
    neuron_diff_pos = io.load_neuron(neuron_name, None, morph_path)
    neuron_diff_neg = io.load_neuron(neuron_name, None, morph_path)
    neuron_new = io.load_neuron(neuron_name, model, new_morph_path)

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    draw_axis(neuron_orig, ax=axs[0, 0])
    axs[0, 0].set_title('Original neuron')

    draw_axis(neuron_new, ax=axs[0, 1])
    axs[0, 1].set_title('New neuron')

    neurite_types = ['basal', 'axon', 'apical']
    for neurite_type in neurite_types:
        neurites = [neurite for neurite in neuron_orig.neurites if neurite.type ==
                    STR_TO_TYPES[neurite_type]]
        neurites_new = [neurite for neurite in neuron_new.neurites if neurite.type ==
                        STR_TO_TYPES[neurite_type]]
        neurites_diff_neg = [
            neurite for neurite in neuron_diff_neg.neurites if neurite.type == STR_TO_TYPES[neurite_type]]
        neurites_diff_pos = [
            neurite for neurite in neuron_diff_pos.neurites if neurite.type == STR_TO_TYPES[neurite_type]]

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
                diff_pos[diff_pos < 0] = 0
                utils.set_diameters(s, diff_pos)

            for j, s in enumerate(iter_sections(neurites_diff_neg[i])):
                diff = diam_new[j] - diam_orig[j]
                diff_neg = -diff.copy()
                diff_neg[diff_neg < 0] = 0
                utils.set_diameters(s, diff_neg)

    draw_axis(neuron_diff_pos, ax=axs[1, 0])
    axs[1, 0].set_title('Positive diameter differences')
    draw_axis(neuron_diff_neg, ax=axs[1, 1])
    axs[1, 1].set_title('Negative diameter differences')

    fig.savefig(folder + '/' + neuron_name + '_' + folder + '_' + model + '.png', dpi=500)
    plt.close('all')


def _create_data(feature1, feature2, original_cells, diametrized_cells, step_size, neurite_types):

    def feature_data(cell, neurite_type):
        nm_neurite_type = STR_TO_TYPES[neurite_type]
        return [neurom.get(feat, cell, neurite_type=nm_neurite_type) for feat in (feature1, feature2)]

    def create_paired_features(cell_list1, cell_list2, neurite_type):
        for cell1, cell2 in zip(cell_list1, cell_list2):
            yield feature_data(cell1, neurite_type), \
                feature_data(cell2, neurite_type)

    def create_bins(step_size, max_value):
        bins = np.arange(0., max_value + step_size, step_size)
        bin_centers = 0.5 * step_size + bins[:-1]
        return bin_centers, bins

    def find_upper_bound(pairs):
        return max(max(vals1.max(), vals2.max())for (vals1, _), (vals2, _) in pairs)

    def per_neurite_data(cell1, cells2, neurite_types):
        for n, neurite_type in enumerate(neurite_types):
            yield list(create_paired_features(original_cells, diametrized_cells, neurite_type))

    assert len(original_cells) == len(diametrized_cells)

    n_cells = len(original_cells)
    iter_neurite_data = per_neurite_data(original_cells, diametrized_cells, neurite_types)
    for n, data_pairs in enumerate(iter_neurite_data):

        upper_bound = find_upper_bound(data_pairs)
        bin_centers, bins = create_bins(step_size, upper_bound)

        stats1 = np.empty((n_cells, len(bin_centers)), dtype=np.float)
        stats2 = np.empty_like(stats1)

        for i, ((metric1, data1), (metric2, data2)) in enumerate(data_pairs):

            res1 = stats.binned_statistic(metric1, data1, statistic='sum',
                                          bins=bins, range=(0, upper_bound))
            res2 = stats.binned_statistic(metric2, data2, statistic='sum',
                                          bins=bins, range=(0, upper_bound))

            stats1[i] = np.cumsum(res1.statistic)
            stats2[i] = np.cumsum(res2.statistic)

        yield bin_centers, stats1, stats2


def plot_cumulative_distribution(original_cells, diametrized_cells, feature1, feature2, neurite_types, step_size=1.0):
    """
    Plot the cumulative distribution of feature2 with respect to the metric values determined via feature1

    Args:
        original_cells: list of NeuroM objects
        diametrized_cells: list of NeuroM objects
            The new cell with the changed diameters.
        neurite_types: string
            The list neurite types to be considered. e.g. ['basal', 'axon']
        step_size: float
            The step size of the cumulative histogram

    Examples of metric features (feature1):
        - segment_radial_distances
        - segment_path_distances (not implemented yet)

    Examples of cumulative distribution features (feature2):
        - segment_volumes
        - segment_surface_areas (not implemented yet)
    """

    assert len(original_cells) == len(diametrized_cells)

    data_generator = _create_data(feature1, feature2, original_cells,
                                  diametrized_cells, step_size, neurite_types)

    f, axes = plt.subplots(2, 1, figsize=(10, 10))

    for i, (bin_centers, stats1, stats2) in enumerate(data_generator):

        color = colors[neurite_types[i]]
        means = stats1.mean(axis=0)

        axes[0].plot(bin_centers, means, c=color, linestyle='-')

        if len(stats1) > 1:  # don't plot std if there is only one curve
            sdevs = stats1.mean(axis=0)
            axes[0].fill_between(bin_centers, means - sdevs, means + sdevs, color=color, alpha=0.2)

        means = stats2.mean(axis=0)

        axes[0].plot(bin_centers, means, c=color, linestyle='--')

        if len(stats2) > 1:
            sdevs = stats2.mean(axis=0)
            axes[0].fill_between(bin_centers, means - sdevs, means + sdevs, color=color, alpha=0.2)

        diffs = np.abs((stats2 - stats1) / stats1)

        diff_means = diffs.mean(axis=0)

        axes[1].plot(bin_centers, diff_means, c=color, linestyle='-')

        if len(diffs) > 1:
            diff_sdevs = diffs.mean(axis=0)
            axes[1].fill_between(bin_centers, diff_means - diff_sdevs,
                                 diff_means + diff_sdevs, color=color, alpha=0.2)
    return f, axes
