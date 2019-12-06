#!/bin/bash

set -e

rm -rf tmp
rm -rf output

module purge all
. ~/diam/bin/activate

module load neuron/7.6.8/python2/serial

export OMP_NUM_THREADS=1

diameter-check run mm-config.json final.json emodel_etype_map.json
diameter-check plot
