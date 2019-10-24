#!/bin/bash

set -e

rm -rf new_neurons
rm -rf original_neurons


module purge all
. ~/diam/bin/activate

module load neuron/7.6.8/python2/serial

export OMP_NUM_THREADS=1

python create_jsons.py

diameter-synthesis run_diameters generate_diameters_params.json

