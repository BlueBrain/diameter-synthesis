""" Type mappings between NeuroM and MorphIO """
import json

import numpy as np
from morphio import SectionType
from neurom import NeuriteType

STR_TO_TYPES = {
    "apical": SectionType.apical_dendrite,
    "basal": SectionType.basal_dendrite,
    "axon": SectionType.axon,
}


TYPE_TO_STR = {
    SectionType.apical_dendrite: "apical",
    SectionType.basal_dendrite: "basal",
    SectionType.axon: "axon",
    SectionType.soma: "soma",
}


NEUROM_TYPE_TO_STR = {
    NeuriteType.apical_dendrite: "apical",
    NeuriteType.basal_dendrite: "basal",
    NeuriteType.soma: "soma",
    NeuriteType.axon: "axon",
}


STR_TO_NEUROM_TYPES = {
    "apical": NeuriteType.apical_dendrite,
    "basal": NeuriteType.basal_dendrite,
    "soma": NeuriteType.soma,
    "axon": NeuriteType.axon,
}


class NumpyEncoder(json.JSONEncoder):
    """To encode numpy arrays"""

    def default(self, o):  # pylint: disable=method-hidden
        """encoder"""
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.integer):
            return int(o)
        return json.JSONEncoder.default(self, o)
