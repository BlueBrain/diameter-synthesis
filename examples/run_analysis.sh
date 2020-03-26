#!/bin/bash

diameter-synthesis run_analysis --config=generate_diameters_params.json --out-dir=./analysis_both --violin
diameter-synthesis run_analysis --config=generate_diameters_params.json --out-dir=./analysis_both --cumulative
