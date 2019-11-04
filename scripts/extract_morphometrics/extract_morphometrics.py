import diameter_synthesis.utils as utils 
import json, os, shutil

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

import neurom as nm

from neurom.core import iter_sections
from neurom import NeuriteType

NEUROM_TYPE_TO_STR = {NeuriteType.apical_dendrite: 'apical',
                      NeuriteType.basal_dendrite: 'basal',
                      NeuriteType.soma: 'soma',
                      NeuriteType.axon: 'axon'}


STR_TO_NEUROM_TYPES = {'apical': NeuriteType.apical_dendrite,
                       'basal': NeuriteType.basal_dendrite,
                       'soma': NeuriteType.soma,
                       'axon': NeuriteType.axon}


def diameter_distance(neuron):
    areas = []
    dists = []
    diams = []
    for neurite in neuron.neurites:
            for neurite_type in neurite_types:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:
                    for sec in iter_sections(neurite):
                        diams += [np.mean(utils.get_diameters(sec))]

                    areas += list(nm.get('section_areas', neurite))
                    dists += list(nm.get('section_path_distances', neurite))

    
    return areas, dists, diams

def plot_morphometrics(morphologies, morphologies_new, folder='figures'):
    for morph in morphologies:

        dists_tot = []
        diams_tot = []
        dists_tot_new = []
        diams_tot_new = []


        if len(morphologies_new[morph])>0:
            plt.figure(figsize=(5,4))
            for i, neuron in enumerate(morphologies[morph]):
                areas, dists, diams = diameter_distance(neuron[0])

                diams_tot += diams 
                dists_tot += dists 

            for i, neuron in enumerate(morphologies_new[morph]):
                areas_new, dists_new, diams_new = diameter_distance(neuron[0])

                diams_tot_new += diams_new 
                dists_tot_new += dists_new 

            #counts, xbins, ybins, images = plt.hist2d(dists_tot, np.log10(diams_tot), bins = 50)#, density=True)#,s=2,label='original')
            counts, xbins, ybins = np.histogram2d(dists_tot, np.log10(diams_tot), bins = 30)#, density=True)#,s=2,label='original')
            cmax = np.max(counts)
            plt.contour(counts.transpose(),extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=1, cmap = plt.cm.rainbow, levels = [0.1*cmax, 0.3*cmax])

            plt.scatter(dists_tot, np.log10(diams_tot), s=0.1, color='k', label='original')

            plt.scatter(dists_tot_new, np.log10(diams_tot_new), s=8,color = 'r', label='new')
            #plt.xscale('log')
            #plt.yscale('log')
            plt.xlabel('path distance')
            plt.ylabel('log10(section mean diameters)')
            plt.legend()
            plt.savefig(folder+'/'+morph+'.png', bbox_inches='tight')
            plt.close()

if __name__ == '__main__':
    config_file = 'config.json'

    with open(config_file, 'r') as f:
        config = json.load(f)


    neurite_types = ['basal',]
    folder = 'figures'

    if os.path.isdir(folder): 
        shutil.rmtree(folder)
        os.mkdir(folder)
    else:
        os.mkdir(folder)

    #Load morphologies
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], by_mtypes = config['by_mtypes'], n_mtypes_max = config['n_mtypes_max'])
    morphologies_new = utils.load_morphologies(config['new_morph_path'], n_morphs_max = config['n_morphs_max'], by_mtypes = config['by_mtypes'], n_mtypes_max = config['n_mtypes_max'], prefix='M1_')

    #plot the diameter vs distance
    plot_morphometrics(morphologies, morphologies_new, folder)
