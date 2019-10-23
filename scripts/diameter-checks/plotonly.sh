#!/bin/bash

set -e

module purge all
. ~/diam/bin/activate

module load neuron/7.6.8/python2/serial

export OMP_NUM_THREADS=1

diameter-check plot
