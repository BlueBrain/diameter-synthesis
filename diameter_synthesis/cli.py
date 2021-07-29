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
@click.option("-e", "--ext", help="figures extention")
def plot_diff(original_folder, diametrized_folder, plot_folder, ncells=None, ext=".png"):
    """Plot original and new neurons as well as their differences."""
    from .plotting import plot_diameter_diff

    neurite_types = ["basal", "apical", "axon"]
    import neurom as nm

    morphologies_dict = create_morphologies_dict(original_folder)

    if not Path(plot_folder).exists():
        os.mkdir(plot_folder)

    for mtype in morphologies_dict:
        L.info("Plot mtype %s", mtype)

        plot_folder_mtype = Path(plot_folder) / mtype
        if not Path(plot_folder_mtype).exists:
            os.mkdir(plot_folder_mtype)

        if ncells is not None:
            neurons = morphologies_dict[mtype][: int(ncells)]
        else:
            neurons = morphologies_dict[mtype]

        for neuron in tqdm(neurons):
            neuron_new = nm.load_morphology(Path(diametrized_folder) / neuron.name)

            plot_diameter_diff(neuron, neuron_new, neurite_types, plot_folder_mtype, ext=ext)


@cli.command("run_analysis")
@click.option("--orig-path", help="Path to original cells", required=True)
@click.option("--diam-path", help="Path to diametrized cells", required=True)
@click.option("--mtypes-file", help="Path to mtypes file", required=False)
@click.option("-o", "--out-dir", help="Directory to output the analysis results", required=True)
@click.option("--cumulative", help="Cumulative distribution plots", is_flag=True)
@click.option("--individual", help="Output a plot for each neuron", is_flag=True)
@click.option("--violin", help="Violin distribution plots", is_flag=True)
@click.option("-e", "--ext", help="Figures extention")
def run_analysis(
    orig_path, diam_path, out_dir, cumulative, individual, violin, mtypes_file=None, ext=".png"
):
    """Produce figures for validation/analysis."""
    if not cumulative and not violin:
        raise ValueError("Should at least set one of '--cumulative' or '--violin' to True.")
    if cumulative:
        from .plotting import cumulative_analysis

        cumulative_analysis(
            orig_path,
            diam_path,
            Path(out_dir) / "basal",
            individual,
            mtypes_file=mtypes_file,
            neurite_types=["basal"],
            ext=ext,
        )
        cumulative_analysis(
            orig_path,
            diam_path,
            Path(out_dir) / "axon",
            individual,
            mtypes_file=mtypes_file,
            neurite_types=["axon"],
            ext=ext,
        )
        cumulative_analysis(
            orig_path,
            diam_path,
            Path(out_dir) / "apical",
            individual,
            mtypes_file=mtypes_file,
            neurite_types=["apical"],
            ext=ext,
        )
    if violin:
        from .plotting import violin_analysis

        violin_analysis(orig_path, diam_path, out_dir, mtypes_file=mtypes_file)
