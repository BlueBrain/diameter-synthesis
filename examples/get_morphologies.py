"""Fetch some barrel morphologies from neuromorpho.org

Copyright (C) 2021  Blue Brain Project, EPFL

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import pandas as pd
import os
import re
import json
from urllib.request import urlopen, Request
from pathlib import Path

NEUROMORPHO_URL = "http://neuromorpho.org"


def get_swc_by_neuron_index(neuronIndex, folder="morphologies"):
    """Download a neuron by index and store it into a SWC file
    Keyword arguments:
    neronIndex -- the neuron index in the database

    Adapted from https://github.com/NeuroBox3D/neuromorpho/blob/master/rest_wrapper/rest_wrapper.py
    """

    url = "%s/api/neuron/id/%i" % (NEUROMORPHO_URL, neuronIndex)
    req = Request(url)
    response = urlopen(req)
    neuron_name = json.loads(response.read().decode("utf-8"))["neuron_name"]
    url = "%s/neuron_info.jsp?neuron_name=%s" % (NEUROMORPHO_URL, neuron_name)
    html = urlopen(url).read().decode("utf-8")
    p = re.compile(r"<a href=dableFiles/(.*)>Morphology File \(Standardized\)</a>", re.MULTILINE)
    m = re.findall(p, html)
    for match in m:
        file_name = match.replace("%20", " ").split("/")[-1]
        response = urlopen("%s/dableFiles/%s" % (NEUROMORPHO_URL, match))
        with open(folder + "/" + file_name, "w") as f:
            f.write(response.read().decode("utf-8"))


if __name__ == "__main__":
    brainRegion = "barrel"
    cell_type = "interneuron"
    numNeurons = 500

    url = "%s/api/neuron/select?q=brain_region:%s&size=%i" % (
        NEUROMORPHO_URL,
        brainRegion,
        numNeurons,
    )
    req = Request(url)
    response = urlopen(req)
    neurons = json.loads(response.read().decode("utf-8"))
    df = pd.DataFrame(neurons["_embedded"]["neuronResources"])
    df["type"] = df["cell_type"].apply(lambda t: t[-1])
    df = df[df["type"] == "principal cell"]
    folder = "morphologies"
    if not Path(folder).exists():
        os.mkdir(folder)
    for gid in df.index:
        get_swc_by_neuron_index(df.loc[gid, "neuron_id"], folder=folder)
