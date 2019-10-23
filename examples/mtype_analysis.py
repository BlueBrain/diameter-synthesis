import os, glob, shutil 
import click
import json
from tqdm import tqdm

import numpy as np
import matplotlib
matplotlib.use('Agg')
import pylab as plt


from diameter_synthesis.plotting import plot_fit_param_boxes

model_file = 'model_params.json'
with open(model_file, 'r') as f: 
    model_params = json.load(f)    

#which model and neurite type to plot
model = 'M0'
neurite_type = 'basal'

plot_fit_param_boxes(model_params, model, neurite_type, 'test', ext = '.png', figsize = (10,4))
