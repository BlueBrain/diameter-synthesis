"""
Logging utilities.
"""
import logging


L = logging.getLogger('diameter-synthesis')


def setup_logging(level):
    """ Setup application logger. """
    logging.basicConfig(
        format="%(asctime)s;%(levelname)s;%(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
        level=logging.WARNING
    )
    L.setLevel(level)
