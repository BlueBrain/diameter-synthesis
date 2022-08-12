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
import json
import os
import re
import ssl
from pathlib import Path
from urllib.request import Request
from urllib.request import urlopen

import pandas as pd

NEUROMORPHO_URL = "http://neuromorpho.org"


def get_url(url):
    """Open a URL without SSL verification."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    ctx.set_ciphers("DEFAULT@SECLEVEL=1")
    req = Request(url)
    response = urlopen(req, context=ctx)  # pylint: disable=consider-using-with
    return response


def get_swc_by_neuron_index(neuronIndex, folder="morphologies"):
    """Download a neuron by index and store it into a SWC file
    Keyword arguments:
    neronIndex -- the neuron index in the database

    Adapted from https://github.com/NeuroBox3D/neuromorpho/blob/master/rest_wrapper/rest_wrapper.py
    """

    url = "%s/api/neuron/id/%i" % (NEUROMORPHO_URL, neuronIndex)
    response = get_url(url)
    neuron_name = json.loads(response.read().decode("utf-8"))["neuron_name"]
    url = "%s/neuron_info.jsp?neuron_name=%s" % (NEUROMORPHO_URL, neuron_name)
    html = get_url(url).read().decode("utf-8")
    p = re.compile(r"<a href=dableFiles/(.*)>Morphology File \(Standardized\)</a>", re.MULTILINE)
    m = re.findall(p, html)
    for match in m:
        file_name = match.replace("%20", " ").split("/")[-1]
        response = get_url("%s/dableFiles/%s" % (NEUROMORPHO_URL, match))
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
    response = get_url(url)
    neurons = json.loads(response.read().decode("utf-8"))
    df = pd.DataFrame(neurons["_embedded"]["neuronResources"])
    df["type"] = df["cell_type"].apply(lambda t: t[-1])
    df = df[df["type"] == "principal cell"]
    folder = "morphologies"
    if not Path(folder).exists():
        os.mkdir(folder)
    for gid in df.index:
        get_swc_by_neuron_index(df.loc[gid, "neuron_id"], folder=folder)
