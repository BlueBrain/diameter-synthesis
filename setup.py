"""Setup for the diameter-synthesis package."""

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from pathlib import Path

from setuptools import find_namespace_packages
from setuptools import setup

reqs = [
    "click>=7.0",
    "jsonschema>=3",
    "matplotlib>=3.4",
    "morphio>=3.3.4",
    "neurom>=3.0,<4.0",
    "numpy>=1.22.0",
    "pandas>=1.0.5",
    "scipy>=1.6",
]

doc_reqs = [
    "m2r2",
    "sphinx",
    "sphinx-bluebrain-theme",
    "sphinx-jsonschema",
    "sphinx_click",
]

test_reqs = [
    "diff_pdf_visually>=1.7.0",
    "decorator>=4",
    "mock>=3",
    "pytest>=6",
    "pytest-click>=1",
    "pytest-console-scripts>=1.3",
    "pytest-cov>=3",
    "pytest-html>=2",
    "pytest-xdist>=2",
]

setup(
    name="diameter-synthesis",
    author="Blue Brain Project, EPFL",
    description="Diametrize cells.",
    long_description=Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url="https://diameter-synthesis.readthedocs.io",
    project_urls={
        "Tracker": "https://github.com/BlueBrain/diameter-synthesis/issues",
        "Source": "git@github.com:BlueBrain/diameter-synthesis.git",
    },
    license="GNU General Public License v3.0",
    packages=find_namespace_packages(include=["diameter_synthesis*"]),
    python_requires=">=3.8",
    use_scm_version=True,
    setup_requires=[
        "setuptools_scm",
    ],
    install_requires=reqs,
    extras_require={
        "docs": doc_reqs,
        "plot": "seaborn>=0.11.1",
        "test": test_reqs,
    },
    entry_points={
        "console_scripts": [
            "diameter-synthesis=diameter_synthesis.cli:main",
        ],
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
