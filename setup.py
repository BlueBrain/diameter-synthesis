"""Setup for the diameter-synthesis package."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

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
    "pandas>=1.1",
    "scipy>=1.6",
    "tqdm>=4.50",
]

doc_reqs = [
    "docutils<0.21",
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
        "Source": "https://github.com/BlueBrain/diameter-synthesis",
    },
    license="Apache License 2.0",
    packages=find_namespace_packages(include=["diameter_synthesis*"]),
    python_requires=">=3.8",
    use_scm_version=True,
    setup_requires=[
        "setuptools_scm",
    ],
    install_requires=reqs,
    extras_require={
        "docs": doc_reqs,
        "plot": "seaborn>=0.13",
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
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
