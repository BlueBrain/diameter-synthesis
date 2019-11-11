#!/bin/bash

set -e

rm -rf figures 

#module purge all
#. ~/diam/bin/activate

#module load neuron/7.6.8/python2/serial

export OMP_NUM_THREADS=1

python create_jsons.py

diameter-synthesis run_models reextract_models_params.json

