"""Click app."""
import logging
import os
from pathlib import Path

import click
import morphio
from tqdm import tqdm

from .utils import create_morphologies_dict

morphio.set_maximum_warnings(0)
L = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


# pylint: disable=import-outside-toplevel,redefined-outer-name


@click.group()
def cli():
    """Cli to learn and generate diameters of neurons."""


@cli.command("run_models")
@click.argument("config_file", type=click.Path(exists=True))
@click.option("--plot", is_flag=True)
@click.option("--ext", default=".png")
def run_models(config_file, plot=False, ext=".png"):
    """Run the model extraction from config file."""
    from .main import run_models

    run_models(config_file, plot=plot, ext=ext)


@cli.command("run_diameters")
@click.argument("config_file", type=click.Path(exists=True))
@click.argument("models_params_file", type=click.Path(exists=True))
def run_diameters(config_file, models_params_file):
    """Build new diameters from config file and diameter model."""
    from .main import run_diameters

    run_diameters(config_file, models_params_file)


@cli.command("plot_diff")
@click.argument("original_folder", type=click.Path(exists=True))
@click.argument("diametrized_folder", type=click.Path(exists=True))
@click.argument("plot_folder", type=click.Path())
@click.option("-n", "--ncells", help="max number of cells to plot")
def plot_diff(original_folder, diametrized_folder, plot_folder, ncells=None):
    """Plot original and new neurons as well as their differences."""
    from .plotting import plot_diameter_diff

    if ncells is not None:
        ncells = int(ncells)
    else:
        ncells = -1

    neurite_types = ["basal", "apical"]
    import neurom as nm

    morphologies_dict = create_morphologies_dict(original_folder)

    if not Path(plot_folder).exists():
        os.mkdir(plot_folder)

    for mtype in morphologies_dict:
        L.info("Plot mtype %s", mtype)

        plot_folder_mtype = Path(plot_folder) / mtype
        if not Path(plot_folder_mtype).exists:
            os.mkdir(plot_folder_mtype)

        for neuron in tqdm(morphologies_dict[mtype][:ncells]):
            neuron_new = nm.load_neuron(Path(diametrized_folder) / neuron.name)

            plot_diameter_diff(
                neuron, neuron_new, neurite_types, plot_folder_mtype, ext=".png",
            )


@cli.command("run_analysis")
@click.option("--config", help="Configuration JSON file", required=True)
@click.option(
    "-o", "--out-dir", help="Directory to output the analysis results", required=True
)
@click.option("--cumulative", help="Cumulative distribution plots", is_flag=True)
@click.option("--individual", help="Output a plot for each neuron", is_flag=True)
@click.option("--violin", help="Violin distribution plots", is_flag=True)
def run_analysis(config, out_dir, cumulative, individual, violin):
    """Produce figures for validation/analysis."""
    if cumulative:
        from .plotting import cumulative_analysis

        cumulative_analysis(config, out_dir, individual)

    if violin:
        from .plotting import violin_analysis

        violin_analysis(config, out_dir)
