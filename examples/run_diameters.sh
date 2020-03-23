#!/bin/bash

python create_jsons.py
diameter-synthesis run_diameters generate_diameters_params.json model_params_mtypes.json

