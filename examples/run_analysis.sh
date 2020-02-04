#!/bin/bash

set -e

diameter-synthesis run_analysis --config=generate_diameters_params.json --out-dir=./analysis_both --violin
