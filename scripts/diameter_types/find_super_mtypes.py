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

def extract_morphometrics(neuron):
    """ extract some morphometrics of a neuron """

    areas = []
    dists = []
    diams = []
    bos = []
    for neurite in neuron.neurites:
            for neurite_type in neurite_types:
                if neurite.type == STR_TO_NEUROM_TYPES[neurite_type]:
                    #for sec in iter_sections(neurite):
                    #    diams += [utils.get_mean_diameter(sec),]

                    areas += list(nm.get('section_areas', neurite))
                    bos += list(nm.get('section_branch_orders', neurite))
                    dists += list(nm.get('section_path_distances', neurite))

                    #diams += morph_funcs.Rall_deviations(neurite)
                    #diams += morph_funcs.sibling_ratios(neurite)
                    diams += morph_funcs.trunk_diameter(neurite)[0]
                    #diams += morph_funcs.terminal_diameters(neurite)

    
    #return areas, dists, diams, bos
    return diams, dists, areas, bos

def morphometric_distances(morphologies):
    """ plot morphometrics per mtypes""" 

    dists_tot = {} 
    bos_tot = {} 
    diams_tot = {} 
    areas_tot = {}

    for morph in morphologies:
        diams_tot[morph] = []
        areas_tot[morph] = []
        dists_tot[morph] = []
        bos_tot[morph] = []

    print('extract parameters')
    for morph in tqdm(morphologies):

        for i, neuron in enumerate(morphologies[morph]):
            areas, dists, diams, bos = extract_morphometrics(neuron[0])

            diams_tot[morph] += diams 
            areas_tot[morph] += areas 
            dists_tot[morph] += dists 
            bos_tot[morph] += bos 

    return diams_tot, areas_tot, dists_tot, bos_tot


if __name__ == '__main__':
    config_file = 'config.json'

    with open(config_file, 'r') as f:
        config = json.load(f)


    shutil.copy(config['morph_path']+'/neuronDB.xml', config['new_morph_path']+'/neuronDB.xml') 
    neurite_types = ['basal', 'apical']

    #Load morphologies
    morphologies = utils.load_morphologies(config['morph_path'], n_morphs_max = config['n_morphs_max'], mtypes_sort = config['mtypes_sort'], n_mtypes_max = config['n_mtypes_max'])

    #plot the diameter vs distance
    diams, areas, dists, bos = morphometric_distances(morphologies)

    diam_distances = np.zeros([len(diams), len(diams)])
    for i, d1 in enumerate(diams.values()):
        for j, d2 in enumerate(diams.values()):
            #diam_distances[i,j] = st.wasserstein_distance(d1, d2) 

            #diam_distances[i,j] += st.wasserstein_distance(list(areas.values())[i], list(areas.values())[j]) 
            diam_distances[i,j] = st.energy_distance(d1, d2) #+ st.energy_distance(list(areas.values())[i], list(areas.values())[j]) 

    import scipy.cluster.hierarchy as hierarchy
    import scipy.spatial.distance as distance

    linked = hierarchy.linkage(distance.squareform(diam_distances), 'ward')

    plt.figure(figsize=(15,15))
    R = hierarchy.dendrogram(linked,
                    orientation='left',
                    labels=list(diams.keys()),
                    no_plot = False,
                    distance_sort='descending',
                    no_labels=False,
                    show_leaf_counts=False)

    plt.savefig('wass_dendro_both.png',bbox_inches='tight')

    ordering = R['leaves']
 
    plt.figure(figsize=(15,15))
    plt.imshow(diam_distances[np.ix_(ordering,ordering)])

    ax = plt.gca()
    ax.set_xticks(np.arange(len(diams)))
    ax.set_xticklabels(np.array(list(diams.keys()))[ordering])
    ax.set_yticks(np.arange(len(diams)))
    ax.set_yticklabels(np.array(list(diams.keys()))[ordering])
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")


    plt.colorbar()
    plt.savefig('wass_dist_both.png',bbox_inches='tight')

    import pygenstability.pygenstability as pgs
    import networkx as nx
    
   # diam_distances[diam_distances<1e-8] = 1e-8  
    print(np.max(diam_distances))
    #th = 0.5*np.max(diam_distances)
    #diam_distances[diam_distances>th] = th
    A = 1/np.sqrt(1.+diam_distances**2)
    np.fill_diagonal(A, 0) #remove diagonal 
    G = nx.Graph(A)

    from RMST import RMST
    len(G.edges())
    G = RMST(G, gamma = 2.0, n_cpu = 1)
    len(G.edges())

    A = nx.to_numpy_matrix(G)
    plt.figure()
    plt.imshow(A) 
    plt.colorbar()
    plt.savefig('A.png')

    louvain_runs = 50
    precision = 1e-8

    #continuous_combinatorial
    stability = pgs.PyGenStability(G, 'continuous_combinatorial', louvain_runs , precision)
    stability.use_spectral_gap = True 

    #number of cpu for parallel compuations
    stability.n_processes_louv = 10
    stability.n_processes_mi = 10

    stability.post_process = False #apply the postprocessing

    #scan over a time interval
    times = np.logspace(-1.5, 0.0, 100)
    stability.n_neigh = len(times) #but here, this is supposed to only look for a few neighbour times for the postprocess

    stability.scan_stability(times, disp=False)

    stability.plot_scan(time_axis=True)

    plt.savefig('louv_clust.png',bbox_inches='tight')

    stability.run_single_stability(time = 10**(-0.7) )
    stability.print_single_result(1, 1)
    
   
    comm_id = np.array(stability.single_stability_result['community_id'])
     
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
    super_mtypes = dict(zip(list(diams.keys()), comm_id_str ))

    with open(config['morph_path']+'/super_mtypes.json', 'w') as json_file:
        json.dump(super_mtypes, json_file, sort_keys=True, indent=4)

    diam_distances = A
    plt.figure(figsize=(15,15))
    
    plt.imshow(diam_distances[np.ix_(ordering,ordering)])

    colors = []     
    cs = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5','C6', 'C7', 'C8',  ]
    for i in comm_id[ordering]:
        colors.append(cs[i % len(cs)])

    ax = plt.gca()
    ax.set_xticks(np.arange(len(diams)))
    ax.set_xticklabels(np.array(list(diams.keys()))[ordering])

    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), colors):
        ticklabel.set_color(tickcolor)

    ax.set_yticks(np.arange(len(diams)))
    ax.set_yticklabels(np.array(list(diams.keys()))[ordering])
    for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
        ticklabel.set_color(tickcolor)

    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", rotation_mode="anchor")


    plt.colorbar()
    plt.savefig('dist_louv.png',bbox_inches='tight')



