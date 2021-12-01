"""Setup for the diameter-synthesis package.

Copyright (C) 2021  Blue Brain Project, EPFL

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
from setuptools import find_packages
from setuptools import setup

# Read the contents of the README file
with open("README.md", encoding="utf-8") as f:
    README = f.read()

reqs = [
    "click>=7.0",
    "jsonschema>=3",
    "matplotlib>=2.2.0",
    "morphio>=2.3.4",
    "neurom>=3.0,<4.0",
    "numpy>=1.15.0",
    "pandas>=0.24.0",
    "scipy>=0.13.3",
]

doc_reqs = [
    "m2r2",
    "sphinx",
    "sphinx-bluebrain-theme",
    "sphinx-jsonschema",
    "sphinx_click",
]

test_reqs = [
    "diff_pdf_visually>=1.5.1",
    "decorator",
    "matplotlib",
    "mock",
    "pytest",
    "pytest-cov",
    "pytest-html",
    "pytest-xdist",
]

setup(
    name="diameter-synthesis",
    author="Blue Brain Project, EPFL",
    description="Diametrize cells",
    long_description=README,
    long_description_content_type="text/markdown",
    license="GPLv3",
    url="https://github.com/BlueBrain/diameter-synthesis",
    project_urls={
        "Tracker": "https://github.com/BlueBrain/diameter-synthesis/issues",
        "Source": "https://github.com/BlueBrain/diameter-synthesis",
    },
    install_requires=reqs,
    extras_require={
        "docs": doc_reqs,
        "plot": "seaborn>=0.11.1",
        "test": test_reqs,
    },
    tests_require=test_reqs,
    entry_points={
        "console_scripts": ["diameter-synthesis=diameter_synthesis.cli:cli"],
    },
    python_requires=">=3.8",
    use_scm_version=True,
    setup_requires=[
        "setuptools_scm",
    ],
    packages=find_packages(exclude="tests"),
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
