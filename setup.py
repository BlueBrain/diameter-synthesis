#!/usr/bin/env python

import imp
import sys

from setuptools import setup, find_packages

if sys.version_info < (2, 7):
    sys.exit("Sorry, Python < 2.7 is not supported")

VERSION = imp.load_source("", "diameter_synthesis/version.py").__version__

setup(
    name="diameter-synthesis",
    author="BlueBrain NSE",
    author_email="bbp-ou-nse@groupes.epfl.ch",
    version=VERSION,
    description="",
    url="https://bbpteam.epfl.ch/project/issues/projects/NSETM/issues",
    download_url="ssh://bbpcode.epfl.ch/nse/model-management-synthesis",
    license="BBP-internal-confidential",
    install_requires=[
        'click>=7.0',
        'scipy>=1.2.0',
        'h5py>=2.9.0',
        'neurom',
        'voxcell>=2.5.5',
    ],
    #dependency_links=[,],
    entry_points={
        'console_scripts': ['diameter-synthesis=diameter_synthesis.main:cli'],
    },
    packages=find_packages(),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
    ],
)
