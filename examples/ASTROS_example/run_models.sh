#!/bin/bash

set -e

rm -rf figures 

#module purge all
#. ~/base/bin/activate

#module load neuron/7.6.8/python2/serial

export OMP_NUM_THREADS=1

python create_jsons.py

#python -m cProfile -o profile_model model_run.py

diameter-synthesis run_models extract_models_params.json

