""" Perform analysis
"""
import os
import click


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
#@click.option("--original-dir", help="Directory with input morphologies", required=True)
#@click.option("--diametrized-dir", help="Directory with diametrized morphologies", required=True)
@click.option("-o", "--out-dir", help='Directory to output the analysis results', required=True)
def cmd(config, out_dir):
    import json
    import pylab as plt
    from diameter_synthesis.io import load_morphology
    from diameter_synthesis.io import iter_morphology_filenames
    from diameter_synthesis.plotting import plot_cumulative_distribution

    def make_figures(original_cells, diametrized_cells, feature1, feature2):

        prefix1, basename1 = _split_prefix(feature1)
        prefix2, basename2 = _split_prefix(feature2)

        assert prefix1 == prefix2

        f, axes = plot_cumulative_distribution(original_cells, diametrized_cells, feature1, feature2, neurite_types)

        #ax.set_xlabel('{}'.format(basename1.replace('_', ' ').title()))
        #ax.set_ylabel('{}'.format(basename2.replace('_', ' ').title()))

        figure_name = 'cumulative_{}_{}_{}'.format(prefix1, basename1, basename2)

        f.savefig(os.path.join(out_dir, figure_name + '.svg'), bbox_inches = 'tight')
        plt.close(f)
        
        if not os.path.isdir(os.path.join(out_dir, figure_name + '_individual')):
            os.mkdir(os.path.join(out_dir, figure_name + '_individual'))  
        for i, (original_cell, diametrized_cell) in enumerate(zip(original_cells, diametrized_cells)):

            f, axes = plot_cumulative_distribution([original_cell], [diametrized_cell], feature1, feature2, neurite_types)

            #ax.set_xlabel('{}'.format(basename1.replace('_', ' ').title()))
            #ax.set_ylabel('{}'.format(basename2.replace('_', ' ').title()))

            fname = '{}_{}.svg'.format(figure_name, original_cell.name)

            f.savefig( os.path.join(out_dir, figure_name + '_individual/' , str(i) + '_' + fname), bbox_inches = 'tight')
            plt.close(f)
    _ensure_dir(out_dir)

    with open(config, 'r') as f:
        config = json.load(f)

    if len(config['models']) > 1:
        print('multiple models provided, will only use the first in the list for analysis')
    model_names = config['models'][0]
    
    original_dir = config['morph_path']
    diametrized_dir = config['new_morph_path']
    neurite_types = config['neurite_types']

    filenames = list(iter_morphology_filenames(original_dir))

    original_filepaths = (os.path.join(original_dir, filename) for filename in filenames if os.path.exists(os.path.join(diametrized_dir, model_names + '_' + filename)))
    diametrized_filepaths = (os.path.join(diametrized_dir, model_names + '_' + filename) for filename in filenames if os.path.exists(os.path.join(diametrized_dir, model_names + '_' + filename)))

    original_cells = list(map(load_morphology, original_filepaths))
    diametrized_cells = list(map(load_morphology, diametrized_filepaths))

    feature_pairs = [
            #('segment_radial_distances', 'segment_volumes'),
            #('section_radial_distances', 'section_areas'),
            ('section_path_distances',   'section_areas'),
            #('section_branch_orders',    'section_areas'),
            #('section_branch_orders',    'section_volumes')
            ]

    for feature1, feature2 in feature_pairs:
        make_figures(original_cells, diametrized_cells, feature1, feature2)

