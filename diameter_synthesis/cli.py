import click

@click.group()
def cli():
    pass

@cli.command('run_models')
@click.argument('config_file', type=click.Path(exists=True))
def run_models(config_file):
    """ Run the model extraction from config file"""
    from .main import run_models
    run_models(config_file)

def run_diameters(config_file):
    """ Build new diameters from config file and diameter model"""
    from .main import run_diameters
    run_diameters(config_file)

