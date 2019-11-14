import diameter_synthesis.utils as utils 
import json, os, shutil

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np

#diable warnings
import morphio
morphio.set_maximum_warnings(0)

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

def extract_morphometrics(neuron, neurite_types):
    """ extract some morphometrics of a neuron """

    areas = []
    dists = []
    diams = []
    bos = []
    for neurite in neuron.neurites:
            for neurite_type in neurite_types:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:
                    for sec in iter_sections(neurite):
                        diams += [np.mean(utils.get_diameters(sec))]

                    areas += list(nm.get('section_areas', neurite))
                    bos += list(nm.get('section_branch_orders', neurite))
                    dists += list(nm.get('section_path_distances', neurite))

    
    return areas, dists, diams, bos

def plot_morphometrics(morphologies, morphologies_new, neurite_types, folder='figures'):
    """ plot morphometrics per mtypes""" 

    for morph in morphologies:

        dists_tot = []
        bos_tot = []
        diams_tot = []
        areas_tot = []

        dists_tot_new = []
        bos_tot_new = []
        diams_tot_new = []
        areas_tot_new = []


        if len(morphologies_new[morph])>0:

            for i, neuron in enumerate(morphologies[morph]):
                areas, dists, diams, bos = extract_morphometrics(neuron[0], neurite_types)

                diams_tot += diams 
                areas_tot += areas 
                dists_tot += dists 
                bos_tot += bos 

            for i, neuron in enumerate(morphologies_new[morph]):
                areas_new, dists_new, diams_new, bos_new = extract_morphometrics(neuron[0], neurite_types)

                diams_tot_new += diams_new 
                areas_tot_new += areas_new 
                dists_tot_new += dists_new 
                bos_tot_new += bos_new 


            plt.figure(figsize=(5,4))
            areas_bo = []
            for bo in set(bos_tot):
                areas_bo.append(np.array(areas_tot)[np.array(bos_tot)==bo])

            plt.boxplot(areas_bo)
            color = dict(color='b')
            plt.boxplot(areas_bo, boxprops=color, capprops=color, whiskerprops=color,flierprops=dict(markeredgecolor='b'), medianprops=color, meanprops=color)

            areas_bo_new = []
            for bo in set(bos_tot_new):
                areas_bo_new.append(np.array(areas_tot_new)[np.array(bos_tot_new)==bo])
            color = dict(color='r')
            plt.boxplot(areas_bo_new, boxprops=color, capprops=color, whiskerprops=color,flierprops=dict(markeredgecolor='r'), medianprops=color, meanprops=color)
            plt.gca().set_xlim([0,np.max(bos_tot)+2])
            plt.axvline(-10,c='b',  label='original diameters')
            plt.axvline(-10,c='r',  label='new diameters')
            plt.xlabel('branching order')
            plt.ylabel('section surface area')
            plt.legend()
            plt.savefig(folder+'/branching_areas/branching_areas_'+morph+'.png', bbox_inches='tight')
            plt.close()


            plt.figure(figsize=(5,4))
            diams_bo = []
            for bo in set(bos_tot):
                diams_bo.append(np.array(diams_tot)[np.array(bos_tot)==bo])

            plt.boxplot(diams_bo)
            color = dict(color='b')
            plt.boxplot(diams_bo, boxprops=color, capprops=color, whiskerprops=color,flierprops=dict(markeredgecolor='b'), medianprops=color, meanprops=color)

            diams_bo_new = []
            for bo in set(bos_tot_new):
                diams_bo_new.append(np.array(diams_tot_new)[np.array(bos_tot_new)==bo])
            color = dict(color='r')
            plt.boxplot(diams_bo_new, boxprops=color, capprops=color, whiskerprops=color,flierprops=dict(markeredgecolor='r'), medianprops=color, meanprops=color)
            plt.gca().set_xlim([0,np.max(bos_tot)+2])
            plt.axvline(-10,c='b',  label='original diameters')
            plt.axvline(-10,c='r',  label='new diameters')
            plt.xlabel('branching order')
            plt.ylabel('section mean diameters')
            plt.legend()
            plt.savefig(folder+'/branching_diams/branching_diams_'+morph+'.png', bbox_inches='tight')
            plt.close()


            plt.figure(figsize=(5,4))

            counts, xbins, ybins = np.histogram2d(dists_tot_new, np.log10(areas_tot_new), bins = 30)
            counts /= np.max(counts)
            plt.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=1, cmap = plt.cm.rainbow, levels = [0.01, 0.1, 0.3])

            plt.colorbar(label='sampled areas')

            plt.scatter(dists_tot, np.log10(areas_tot), s=0.1, color='k', label='original surface areas')

            plt.axvline(100)
            plt.xlabel('path distance')
            plt.ylabel('log10(section mean surface area)')
            plt.legend()
            plt.savefig(folder+'/path_areas/path_areas_'+morph+'.png', bbox_inches='tight')
            plt.close()


            plt.figure(figsize=(5,4))

            counts, xbins, ybins = np.histogram2d(dists_tot_new, np.log10(diams_tot_new), bins = 30)
            counts /= np.max(counts)
            plt.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]],
                linewidths=1, cmap = plt.cm.rainbow, levels = [0.01, 0.1, 0.3])

            plt.colorbar(label='sampled diameters')

            plt.scatter(dists_tot, np.log10(diams_tot), s=0.1, color='k', label='original diameters')

            plt.axvline(100)
            plt.xlabel('path distance')
            plt.ylabel('log10(section mean diameters)')
            plt.legend()
            plt.savefig(folder+'/path_diams/path_diams_'+morph+'.png', bbox_inches='tight')
            plt.close()

        else:
            print("len=0")

if __name__ == '__main__':
    config_file = 'config.json'

    with open(config_file, 'r') as f:
        config = json.load(f)


    shutil.copy(config['morph_path']+'/neuronDB.xml', config['new_morph_path']+'/neuronDB.xml') 
    neurite_types = config['neurite_types']
    folder = 'figures_'+config['new_morph_path']

    if os.path.isdir(folder): 
        shutil.rmtree(folder)
        os.mkdir(folder)
        os.mkdir(folder+'branching_diams')
        os.mkdir(folder+'branching_areas')
        os.mkdir(folder+'path_diams')
        os.mkdir(folder+'path_areas')
    else:
        os.mkdir(folder)
        os.mkdir(folder+'branching_diams')
        os.mkdir(folder+'branching_areas')
        os.mkdir(folder+'path_diams')
        os.mkdir(folder+'path_areas')

    #Load morphologies
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'])
    morphologies_new = {}
    for prefix in config['models']:
        print(prefix)
        morphologies_new_tmp = utils.load_morphologies(config['new_morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'], prefix=prefix+'_')
        for morph in morphologies_new_tmp:
            if morph in morphologies_new:
                morphologies_new[morph] += morphologies_new_tmp[morph]
            else:
                morphologies_new[morph] = morphologies_new_tmp[morph]

    #plot the diameter vs distance
    plot_morphometrics(morphologies, morphologies_new, neurite_types, folder)
