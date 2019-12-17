import numpy as np
import matplotlib.pyplot as plt

import bglibpy
import neurom as nrm
from bglibpy import bluepy
import os, re, json

from utils import *

bglibpy.set_verbose(100)

morph_dir_bio = '../../examples/05_RepairUnravel-asc/'
morph_dir_diam = 'new_neurons/'


with open('gids.json', 'r') as f:
    gids = json.load(f)

neurite_type = 'apical'

errors_impedance = []
errors_transfer = []
k = 0 

for cell_name in gids:
    if k==100:
        break
    k+=1
    gid = gids[cell_name][0]

    is_transfer = True

    folder_traces = 'figure_traces'
    if not os.path.isdir(folder_traces):
        os.mkdir(folder_traces)

    data_bio = run(morph_dir_bio, gid, cell_name, neurite_type=neurite_type, is_transfer=is_transfer, tpe='bio', folder=folder_traces)
    data_diam = run(morph_dir_diam, gid, cell_name, neurite_type=neurite_type, is_transfer=is_transfer, tpe='diam', folder=folder_traces)

    folder = 'figure_transfer_'+neurite_type
    if not os.path.isdir(folder):
        os.mkdir(folder)

    error_transfer = plot_impedance(data_bio, data_diam, cell_name, folder, log=True)
    errors_transfer.append(error_transfer)

    is_transfer = False

    data_bio = run(morph_dir_bio, gid, cell_name, neurite_type=neurite_type, is_transfer=is_transfer, tpe='bio', folder=folder_traces)
    data_diam = run(morph_dir_diam, gid, cell_name, neurite_type=neurite_type, is_transfer=is_transfer, tpe='diam', folder=folder_traces)

    folder = 'figure_impedance_'+neurite_type
    if not os.path.isdir(folder):
        os.mkdir(folder)

    error_impedance = plot_impedance(data_bio, data_diam, cell_name, folder, log=True)
    errors_impedance.append(error_impedance)


    plt.figure()
    plt.hist(errors_impedance, bins=50)
    plt.xlabel('L2 error on impedances')
    plt.savefig('errors_impedances.png', bbox_inches= 'tight')

    plt.figure()
    plt.hist(errors_transfer, bins=50)
    plt.xlabel('L2 error on transfer impedances')
    plt.savefig('errors_transfer.png', bbox_inches= 'tight')

    plt.figure()
    plt.scatter(errors_impedance, errors_transfer)
    plt.xlabel('L2 error on impedances')
    plt.ylabel('L2 error on transfer impedances')
    plt.savefig('errors_both_scatter.png', bbox_inches= 'tight')
