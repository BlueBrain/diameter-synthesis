""" Type mappings between NeuroM and MorphIO """
from neurom import NeuriteType
from morphio import SectionType


STR_TO_TYPES = {'apical': SectionType.apical_dendrite,
                'basal': SectionType.basal_dendrite,
                'axon': SectionType.axon}


TYPE_TO_STR = {SectionType.apical_dendrite: 'apical',
               SectionType.basal_dendrite: 'basal',
               SectionType.axon: 'axon',
               SectionType.soma: 'soma'}


NEUROM_TYPE_TO_STR = {NeuriteType.apical_dendrite: 'apical',
                      NeuriteType.basal_dendrite: 'basal',
                      NeuriteType.soma: 'soma',
                      NeuriteType.axon: 'axon'}

STR_TO_NEUROM_TYPES = {'apical': NeuriteType.apical_dendrite,
                       'basal': NeuriteType.basal_dendrite,
                       'soma': NeuriteType.soma,
                       'axon': NeuriteType.axon}
