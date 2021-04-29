Diameter synthesis
==================

This code is a work in progress and will aimed at generating synthetic diameters for neurons, with parameters learned from a set of biological neurons.


Installation
------------

The usual

.. code:: bash

	pip install -e .

Main usage
-----------

Step 1: Building models
~~~~~~~~~~~~~~~~~~~~~~~

In folder `example`, you first have to modify `create_jsons.py` to suit your needs.

You have the following important parameters for the dict `extract_models_params`:

- `morph_path`: path to morphology files
- `mtypes_sort`: how to learn distributions: `all` to use all together, `mtypes` to use by mtypes , `super_mtypes` to use home made cells types (see `diameter_types` below)
- `models`: to create several models (for now they are all the same, just differen realisation of random numbers)
- `neurite_types`: types of neurite to learn parameters for
- `extra_params`: dict of additional model parameters

Step 2: Building diameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

Then simply run `./run_models.sh` to create the models (saved in a json file).

In `create_jsons.py`, the dict `generate_diameters_params` needs to be updated, too, with entries matching the previous dict.
The path in `new_morph_path` will be where the new morphologies will be saved.

Then run `./run_diamters.sh` to generate diameters.


Additional scripts
------------------

Several additional scripts in folder `scripts`:

- `diameter-checks`: run the diameter-check code (bluepymm) on the biological and sampled cells
- `diameter_types`: cluster mtypes using distributions of surface areas (uses two privates repositories a the moment)
- `extract_morphometrics`: from bio and sample cells, extracts and plot distribution of surface aread and diameter as a function of branch order and path lengths
- `extract_morphologies`: from a cell release, find the ones that can be run through diameter-check
- `plot_morphologies`: plot all morphologies in mtype folders

