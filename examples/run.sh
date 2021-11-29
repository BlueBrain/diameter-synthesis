#!/bin/bash

# fetch some morphologies from neuromorpho.org
python get_morphologies.py

# learn diameter model
diameter-synthesis run_models model_config.json

# rediametrize morphologies
diameter-synthesis run_diameters diametrizer_config.json diameter_model.json

# plot differences between diameterds
diameter-synthesis  plot_diff --orig-path=morphologies --diam-path=diametrized_morphologies --out-dir=./analysis

# plot morphometrics
diameter-synthesis run_analysis --orig-path=morphologies --diam-path=diametrized_morphologies --out-dir=./analysis --violin

# plot cumulatice diameter distributions
diameter-synthesis run_analysis --orig-path=morphologies --diam-path=diametrized_morphologies --out-dir=./analysis --cumulative
