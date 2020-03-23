#!/bin/bash
python create_jsons.py
diameter-synthesis run_models extract_models_params.json --plot

