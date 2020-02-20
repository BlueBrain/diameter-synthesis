""" click module """
import click
import os
from tqdm import tqdm

from .main import run_models as run_models_main
from .main import run_diameters as run_diameters_main
from .utils import create_morphologies_dict
from .io import load_morphology


@click.group()
def cli():
    """ A tool to learn and generate diameters of neurons """


@cli.command("run_models")
@click.argument("config_file", type=click.Path(exists=True))
def run_models(config_file):
    """ Run the model extraction from config file"""
    run_models_main(config_file)


@cli.command("run_diameters")
@click.argument("config_file", type=click.Path(exists=True))
def run_diameters(config_file):
    """ Build new diameters from config file and diameter model"""
    run_diameters_main(config_file)


@cli.command("run_analysis")
@click.option("--config", help="Configuration JSON file", required=True)
@click.option(
    "-o", "--out-dir", help="Directory to output the analysis results", required=True
)
@click.option("--cumulative", help="Cumulative distribution plots", is_flag=True)
@click.option("--individual", help="Output a plot for each neuron", is_flag=True)
@click.option("--violin", help="Violin distribution plots", is_flag=True)
def run_analysis(config, out_dir, cumulative, individual, violin):
    """produce figures for validation/analysis"""
    if cumulative:
        from .plotting import (
            cumulative_analysis,
        )  # pylint: disable=import-outside-toplevel

        cumulative_analysis(config, out_dir, individual)
    if violin:
        from .plotting import violin_analysis  # pylint: disable=import-outside-toplevel

        violin_analysis(config, out_dir)


@cli.command("plot_diff")
@click.argument("original_folder", type=click.Path(exists=True))
@click.argument("diametrized_folder", type=click.Path(exists=True))
@click.argument("plot_folder", type=click.Path())
@click.option("-n", "--ncells", help="Violin distribution plots")
def plot_diff(original_folder, diametrized_folder, plot_folder, ncells=10):
    """produce figures for validation/analysis"""
    from .plotting import plot_diameter_diff  # pylint: disable=import-outside-toplevel

    ncells = int(ncells)

    morphologies_dict_original = create_morphologies_dict(
        original_folder, mtypes_sort="mtype"
    )

    morphologies_dict_diametrized = create_morphologies_dict(
        diametrized_folder, mtypes_sort="mtype"
    )

    for mtype in tqdm(morphologies_dict_original):
        count = 0
        for neuron_name in morphologies_dict_original[mtype]:
            neuron_new = load_morphology(os.path.join(diametrized_folder, neuron_name))
            neurite_types = ["basal", "apical"]
            plot_diameter_diff(
                os.path.splitext(neuron_name)[0],
                original_folder,
                neuron_new,
                "",
                neurite_types,
                plot_folder,
                ext=".png",
            )

            count += 1
            if count > ncells:
                break
