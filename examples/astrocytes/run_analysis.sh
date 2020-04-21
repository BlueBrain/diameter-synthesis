#!/bin/bash

#diameter-synthesis run_analysis --orig-path=./morphologies --diam-path=./diametrized_morphologies --out-dir=./analysis_both --violin
diameter-synthesis run_analysis --orig-path=./morphologies --diam-path=./diametrized_morphologies --out-dir=./analysis_both --cumulative --individual
