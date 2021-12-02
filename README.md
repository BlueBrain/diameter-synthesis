[![Version](https://img.shields.io/pypi/v/diameter-synthesis)](https://github.com/BlueBrain/diameter-synthesis/releases)
[![Build status](https://github.com/BlueBrain/diameter-synthesis/actions/workflows/run-tox.yml/badge.svg?branch=main)](https://github.com/BlueBrain/diameter-synthesis/actions)
[![Codecov.io](https://codecov.io/github/BlueBrain/diameter-synthesis/coverage.svg?branch=main)](https://codecov.io/github/BlueBrain/diameter-synthesis?branch=main)
[![License](https://img.shields.io/badge/License-GPLv3-blue)](https://github.com/BlueBrain/diameter-synthesis/blob/main/LICENSE.txt)
[![Documentation status](https://readthedocs.org/projects/diameter-synthesis/badge/?version=latest)](https://diameter-synthesis.readthedocs.io/)
[![DOI](https://img.shields.io/badge/DOI-10.1101/2020.04.15.040410-blue)](https://doi.org/10.1101/2020.04.15.040410)


# Diameter synthesis

This code aims at generating synthetic diameters for neurons, with parameters learned from a set of biological neurons.


## Installation

Use pip:

```bash
pip install diameter-synthesis
```

## Main usage

### Step 1: Building models

In folder `example`, you first have to modify `create_jsons.py` to suit your needs.

You have the following important parameters for the dict `extract_models_params`:

- `morph_path`: path to morphology files
- `mtypes_sort`: how to learn distributions: `all` to use all together, `mtypes` to use by mtypes , `super_mtypes` to use home made cells types (see `diameter_types` below)
- `models`: to create several models (for now they are all the same, just different realisation of random numbers)
- `neurite_types`: types of neurite to learn parameters for
- `extra_params`: dict of additional model parameters

### Step 2: Building diameters

Then simply run `./run_models.sh` to create the models (saved in a json file).

In `create_jsons.py`, the dict `generate_diameters_params` needs to be updated, too, with entries matching the previous dict.
The path in `new_morph_path` will be where the new morphologies will be saved.

Then run `./run_diamters.sh` to generate diameters.


## Additional scripts

Several additional scripts in folder `scripts`:

- `diameter-checks`: run the diameter-check code (bluepymm) on the biological and sampled cells
- `diameter_types`: cluster mtypes using distributions of surface areas (uses two privates repositories a the moment)
- `extract_morphometrics`: from bio and sample cells, extracts and plot distribution of surface area and diameter as a function of branch order and path lengths
- `extract_morphologies`: from a cell release, find the ones that can be run through diameter-check
- `plot_morphologies`: plot all morphologies in mtype folders


## Examples

The `examples` folder contains a simple example that will fetch morphologies from [neuromorpho.org](http://neuromorpho.org), learn a diameter model, rediametrize these morphologies, and perform some analysis of the results to compare original and diametrized morphologies.
This example can simply be run using the following command:
```bash
./run.sh
```


## Funding & Acknowledgment

The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government’s ETH Board of the Swiss Federal Institutes of Technology.

For license and authors, see `LICENSE.txt` and `AUTHORS.md` respectively.

Copyright © 2021 Blue Brain Project/EPFL
