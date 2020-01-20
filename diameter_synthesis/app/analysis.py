""" Perform analysis """
import json
import os

import click
import pylab as plt

from diameter_synthesis.io import iter_morphology_filenames, load_morphology
from diameter_synthesis.plotting import plot_cumulative_distribution


def _ensure_dir(path):

    if not os.path.exists(path):
        os.makedirs(path)


def _split_prefix(neurom_feature_name):

    name_list = neurom_feature_name.split('_')

    prefix = name_list[0]
    basename = '_'.join(name_list[1:])

    return prefix, basename


@click.command(help=__doc__)
@click.option("--config", help="Configuration JSON file", required=True)
@click.option("-o", "--out-dir", help='Directory to output the analysis results', required=True)
def cmd(config, out_dir):
    """ cmd function """

    def make_figures(original_cells, diametrized_cells, feature1, feature2):

        prefix1, basename1 = _split_prefix(feature1)
        prefix2, basename2 = _split_prefix(feature2)

        assert prefix1 == prefix2

        fig, _ = plot_cumulative_distribution(
            original_cells, diametrized_cells, feature1, feature2, config['neurite_types'])

        figure_name = 'cumulative_{}_{}_{}'.format(prefix1, basename1, basename2)

        fig.savefig(os.path.join(out_dir, figure_name + '.svg'), bbox_inches='tight')
        plt.close(fig)

        if not os.path.isdir(os.path.join(out_dir, figure_name + '_individual')):
            os.mkdir(os.path.join(out_dir, figure_name + '_individual'))

        # for i, (original_cell, diametrized_cell) in
        #         enumerate(zip(original_cells, diametrized_cells)):
        #     f, axes = plot_cumulative_distribution([original_cell], [diametrized_cell],
        #         feature1, feature2, config['neurite_types'], auto_limit=False)
        #     fname = '{}_{}.svg'.format(figure_name, original_cell.name)
        #     f.savefig( os.path.join(out_dir, figure_name + '_individual/' , str(i) + '_' + fname),
        #         bbox_inches = 'tight')
        #     plt.close(f)

    with open(config, 'r') as filename:
        config = json.load(filename)

    out_dir += '_' + config['mtypes_sort']
    _ensure_dir(out_dir)

    if len(config['models']) > 1:
        print('multiple models provided, will only use the first in the list for analysis')
    # model_names = config['models'][0]

    filenames_original = list(iter_morphology_filenames(config['morph_path']))
    if len(list(filenames_original)) == 0:  # to deal with cells in mtypes folder
        filenames_original = []
        for folder in os.listdir(config['morph_path']):
            for fname in os.listdir(os.path.join(config['morph_path'], folder)):
                if fname.endswith(('.h5', '.asc', '.swc')):
                    filenames_original.append(os.path.join(folder, fname))

    filenames_diametrized = list(iter_morphology_filenames(config['new_morph_path']))

    # if os.path.exists(os.path.join(config['new_morph_path'], model_names + '_' + filename)))
    original_filepaths = (os.path.join(config['morph_path'], filename)
                          for filename in filenames_original)
    # if os.path.exists(os.path.join(config['new_morph_path'], filename)))
    diametrized_filepaths = (os.path.join(config['new_morph_path'], filename)
                             for filename in filenames_diametrized)

    print('Loading neurons...')
    original_cells = list(map(load_morphology, original_filepaths))
    diametrized_cells = list(map(load_morphology, diametrized_filepaths))
    print('...done.')
    # sort cells to ensure they are ordered the same
    original_cells = sorted(original_cells, key=lambda neuron: neuron.name)
    diametrized_cells = sorted(diametrized_cells, key=lambda neuron: neuron.name)

    feature_pairs = [
        ('segment_radial_distances', 'segment_volumes'),
        # ('section_radial_distances', 'section_areas'),
        ('section_path_distances', 'section_areas'),
        # ('section_branch_orders',    'section_areas'),
        # ('section_branch_orders',    'section_volumes')
    ]

    for feature1, feature2 in feature_pairs:
        make_figures(original_cells, diametrized_cells, feature1, feature2)
