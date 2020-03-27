#!/bin/bash

python create_jsons.py
diameter-synthesis run_diameters diametrizer_params.json model_params_mtypes.json

