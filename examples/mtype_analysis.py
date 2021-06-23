""" """

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, glob, shutil
import click
import json
from tqdm import tqdm

import numpy as np
import matplotlib

matplotlib.use("Agg")
import pylab as plt


from diameter_synthesis.plotting import plot_fit_param_boxes

model_file = "model_params.json"
with open(model_file, "r") as f:
    model_params = json.load(f)

# which model and neurite type to plot
model = "M0"
neurite_type = "basal"

plot_fit_param_boxes(
    model_params, model, neurite_type, "test", ext=".png", figsize=(10, 4)
)
