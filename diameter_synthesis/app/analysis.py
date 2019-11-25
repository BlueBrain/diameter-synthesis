""" Perform analysis
"""
import os
import click


def _ensure_dir(path):

    if not os.path.exists(path):
        os.makedirs(path)


@click.command(help=__doc__)
@click.option("--config", help="Configuration JSON file", required=True)
@click.option("--original-dir", help="Directory with input morphologies", required=True)
@click.option("--diametrized-dir", help="Directory with diametrized morphologies", required=True)
@click.option("-o", "--out-dir", help='Directory to output the analysis results', required=True)
def cmd(config, original_dir, diametrized_dir, out_dir):
    import json
    from diameter_synthesis.io import load_morphology
    from diameter_synthesis.io import iter_morphology_filenames
    from diameter_synthesis.plotting import plot_cumulative_distribution

    _ensure_dir(out_dir)

    with open(config, 'r') as f:
        config = json.load(f)

    model_names = set(config['models'])
    neurite_types = config['neurite_types']

    filenames = list(iter_morphology_filenames(original_dir))

    original_filepaths = (os.path.join(original_dir, filename) for filename in filenames)
    diametrized_filepaths = (os.path.join(diametrized_dir, 'M0_' + filename) for filename in filenames)


    original_cells = list(map(load_morphology, original_filepaths))
    diametrized_cells = list(map(load_morphology, diametrized_filepaths))

    feature1 = 'segment_radial_distances'
    feature2 = 'segment_volumes'

    f, ax = plot_cumulative_distribution(original_cells, diametrized_cells, feature1, feature2, neurite_types)

    ax.set_xlabel('Radial Distance (um)')
    ax.set_ylabel('Volume um3')

    figure_name = 'cumulative_radial_distances.png'

    f.savefig(os.path.join(out_dir, figure_name))

    for original_cell, diametrized_cell in zip(original_cells, diametrized_cells):
        f, ax = plot_cumulative_distribution([original_cell], [diametrized_cell], feature1, feature2, neurite_types)

        figure_name = 'cumulative_radial_distances_' + original_cell.name + '.png'

        f.savefig(os.path.join(out_dir, figure_name))
