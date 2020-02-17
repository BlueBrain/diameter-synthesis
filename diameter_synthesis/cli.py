""" click module """
import click

from .main import run_models as run_models_main
from .main import run_diameters as run_diameters_main


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
        from .plotting import cumulative_analysis  # pylint: disable=import-outside-toplevel
        cumulative_analysis(config, out_dir, individual)
    if violin:
        from .plotting import violin_analysis  # pylint: disable=import-outside-toplevel
        violin_analysis(config, out_dir)
