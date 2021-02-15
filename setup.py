#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

if sys.version_info < (2, 7):
    sys.exit("Sorry, Python < 2.7 is not supported")

VERSION = imp.load_source("", "diameter_synthesis/version.py").__version__

setup(
    name="diameter-synthesis",
    author="BlueBrain Cells",
    author_email="bbp-ou-cells@groupes.epfl.ch",
    version=VERSION,
    description="",
    license="BBP-internal-confidential",
    install_requires=[
        "click>=7.0",
        "matplotlib>=2.2.0",
        "morphio>=2.3.4",
        "neurom>=1.4.2",
        "numpy>=1.15.0",
        "pandas>=0.24.0",
        "scipy>=0.13.3",
        "seaborn",
    ],
    entry_points={
        "console_scripts": ["diameter-synthesis=diameter_synthesis.cli:cli"],
    },
    packages=find_packages(exclude="tests"),
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
    ],
)
