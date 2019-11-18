import diameter_synthesis.utils as utils 
import json, os, shutil

import matplotlib
matplotlib.use('Agg')
import pylab as plt
import numpy as np
import scipy.stats as st

from tqdm import tqdm

#diable warnings
import morphio
morphio.set_maximum_warnings(0)

import neurom as nm

from neurom.core import iter_sections
from neurom import NeuriteType

import diameter_synthesis.morph_functions as morph_funcs

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
    for neurite in neuron.neurites:
            for neurite_type in neurite_types:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:

                    areas += list(nm.get('section_areas', neurite))
                    dists += list(nm.get('section_path_distances', neurite))

    return dists, areas

def morphometric_distances(morphologies, neurite_types):
    """ plot morphometrics per mtypes""" 

    dists_tot = {} 
    areas_tot = {}

    for morph in morphologies:
        areas_tot[morph] = []
        dists_tot[morph] = []

    print('extract parameters')
    for morph in tqdm(morphologies):

        for i, neuron in enumerate(morphologies[morph]):
            dists, areas = extract_morphometrics(neuron[0], neurite_types)

            dists_tot[morph] += dists 
            areas_tot[morph] += areas 

    return  dists_tot, areas_tot


if __name__ == '__main__':
    config_file = 'config.json'

    with open(config_file, 'r') as f:
        config = json.load(f)

    shutil.copy(config['morph_path']+'/neuronDB.xml', config['new_morph_path']+'/neuronDB.xml') 

    #Load morphologies
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'])

    #plot the areas vs path distances
    dists, areas = morphometric_distances(morphologies, config['neurite_types'])

    area_distances = np.zeros([len(areas), len(areas)])

    if not os.path.isdir('distributions'):
        os.mkdir('distributions')
    #for i, val1 in enumerate(areas.values()):
    for i, mtype in enumerate(areas):
        val1 = areas[mtype]
        val1 = np.array(val1)[np.array(dists[mtype])<100]

        plt.figure()
        plt.hist(val1, range = (0,10**2.5), density=True, bins = 50) 
        plt.savefig('distributions/dist_'+mtype+'.png')
        plt.close()

        for j, mtype2 in enumerate(areas):
            #diam_distances[i,j] = st.wasserstein_distance(val1, val2) 
            val2 = areas[mtype2]
            val2 = np.array(val2)[np.array(dists[mtype2])<100]
            area_distances[i,j] = st.energy_distance(val1, val2) 

    import scipy.cluster.hierarchy as hierarchy
    import scipy.spatial.distance as distance

    linked = hierarchy.linkage(distance.squareform(area_distances), 'ward')

    plt.figure(figsize=(15,15))
    R = hierarchy.dendrogram(linked,
                    orientation='left',
                    labels=list(areas.keys()),
                    no_plot = False,
                    distance_sort='descending',
                    no_labels=False,
                    show_leaf_counts=False)

    plt.savefig('wass_dendro_both.png',bbox_inches='tight')

    ordering = R['leaves']
 
    plt.figure(figsize=(15,15))
    plt.imshow(area_distances[np.ix_(ordering,ordering)])

    ax = plt.gca()
    ax.set_xticks(np.arange(len(areas)))
    ax.set_xticklabels(np.array(list(areas.keys()))[ordering])
    ax.set_yticks(np.arange(len(areas)))
    ax.set_yticklabels(np.array(list(areas.keys()))[ordering])
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")


    plt.colorbar()
    plt.savefig('wass_dist_both.png',bbox_inches='tight')

    import pygenstability.pygenstability as pgs
    import networkx as nx
    
    #A = 1/np.sqrt(1.+area_distances**2)
    A = 1/(1.+area_distances)
    np.fill_diagonal(A, 0) #remove diagonal 
    G = nx.Graph(A)

    from RMST import RMST
    G = RMST(G, gamma = 0.5, n_cpu = 1)

    A = nx.to_numpy_matrix(G)
    plt.figure()
    plt.imshow(A) 
    plt.colorbar()
    plt.savefig('A.png')

    louvain_runs = 100
    precision = 1e-8

    #continuous_combinatorial
    stability = pgs.PyGenStability(G, 'continuous_combinatorial', louvain_runs , precision)
    stability.use_spectral_gap = True 

    #number of cpu for parallel compuations
    stability.n_processes_louv = 5
    stability.n_processes_mi = 5

    #scan over a time interval
    times = np.logspace(-1.0, 0.5, 100)

    stability.post_process = False #apply the postprocessing
    stability.n_neigh = len(times) #but here, this is supposed to only look for a few neighbour times for the postprocess

    stability.scan_stability(times, disp=False)

    stability.plot_scan(time_axis=True)

    plt.savefig('louv_clust.png',bbox_inches='tight')

    stability.run_single_stability(time = 10**(-0.2) )
    stability.print_single_result(1, 1)
   
    comm_id = np.array(stability.single_stability_result['community_id'])
    
    for i in set(comm_id):
        ids = np.argwhere(comm_id==i).flatten()

        if len(ids)==1: 
            if i < np.max(comm_id):
                comm_id[ids] = i+1 
            else:
                comm_id[ids] = i-1

    ordering = [] 
    for i in set(comm_id):
            ordering += list(np.argwhere(comm_id==i).flatten())

    super_mtypes = {}
    for i in set(comm_id):
        super_mtypes[str(i)] = []
    
    #convert the community id to string to be used as super mtypes
    comm_id_str = []
    for c in comm_id:
        comm_id_str.append(str(c))
    super_mtypes = dict(zip(list(areas.keys()), comm_id_str ))

    with open('super_mtypes.json', 'w') as json_file:
        json.dump(super_mtypes, json_file, sort_keys=True, indent=4)

    plt.figure(figsize=(15,15))
    
    plt.imshow(A[np.ix_(ordering,ordering)])

    colors = []     
    cs = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5','C6', 'C7', 'C8',  ]
    for i in comm_id[ordering]:
        colors.append(cs[i % len(cs)])

    ax = plt.gca()
    ax.set_xticks(np.arange(len(areas)))
    ax.set_xticklabels(np.array(list(areas.keys()))[ordering])

    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), colors):
        ticklabel.set_color(tickcolor)

    ax.set_yticks(np.arange(len(areas)))
    ax.set_yticklabels(np.array(list(areas.keys()))[ordering])
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
        ticklabel.set_color(tickcolor)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")


    plt.colorbar()
    plt.savefig('dist_louv.png',bbox_inches='tight')



