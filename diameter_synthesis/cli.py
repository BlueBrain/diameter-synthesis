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
