""" Long description diameter synthesis
"""

import logging
import click

from diameter_synthesis.app import analysis

from diameter_synthesis.app.logger import setup_logging
#from diameter_synthesis.version import VERSION

#diable warnings
import morphio
morphio.set_maximum_warnings(0)



@click.group('diameter-synthesis', help=__doc__.format(esc='\b'))
@click.option("-v", "--verbose", count=True, help="-v for INFO, -vv for DEBUG")
#@click.version_option(VERSION)
def app(verbose=0):
    # pylint: disable=missing-docstring
    level = {
        0: logging.WARNING,
        1: logging.INFO,
        2: logging.DEBUG,
    }[verbose]
    setup_logging(level)


app.add_command(name='analysis', cmd=analysis.cmd)
