"""Setup for the diameter-synthesis package."""
import imp
import sys

from setuptools import find_packages
from setuptools import setup

if sys.version_info < (2, 7):
    sys.exit("Sorry, Python < 2.7 is not supported")

VERSION = imp.load_source("", "diameter_synthesis/version.py").__version__

test_reqs = [
    "diff_pdf_visually>=1.5.1",
    "matplotlib<=3.3.4",
    "decorator",
    "matplotlib<3.4",  # PDF rendering is slightly modified with matplotlib >= 3.4
    "mock",
    "pytest",
    "pytest-cov",
    "pytest-html",
    "pytest-xdist",
]

setup(
    name="diameter-synthesis",
    author="BlueBrain Cells",
    author_email="bbp-ou-cells@groupes.epfl.ch",
    version=VERSION,
    description="Diametrize cells",
    license="BBP-internal-confidential",
    url="https://bbpteam.epfl.ch/documentation/projects/diameter-synthesis",
    project_urls={
        "Tracker": "https://bbpteam.epfl.ch/project/issues/projects/CELLS/issues",
        "Source": "ssh://bbpcode.epfl.ch/cells/diameter-synthesis",
    },
    install_requires=[
        "click>=7.0",
        "jsonschema>=3",
        "matplotlib>=2.2.0",
        "morphio>=2.3.4",
        "neurom>=3.0,<4.0",
        "numpy>=1.15.0",
        "pandas>=0.24.0",
        "scipy>=0.13.3",
    ],
    extras_require={
        "plot": "seaborn>=0.11.1",
        "test": test_reqs,
    },
    tests_require=test_reqs,
    entry_points={
        "console_scripts": ["diameter-synthesis=diameter_synthesis.cli:cli"],
    },
    python_requires=">=3.6",
    packages=find_packages(exclude="tests"),
    include_package_data=True,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
