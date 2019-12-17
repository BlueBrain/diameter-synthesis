import numpy as np
import matplotlib.pyplot as plt

import bglibpy
import neurom as nrm
import os, re

from diameter_synthesis.utils import set_bins 

bglibpy.set_verbose(100)

def compute_impedence(sim, icell, freq=100, neurite_type='basal', is_transfer=True):
    area = 0
    imp = sim.neuron.h.Impedance()

    if neurite_type == 'basal':
        sections = icell.basal
    elif neurite_type == 'apical':
        sections = icell.apical
    else:
        raise Exception('Unknown neurite type!')

    areas = []
    distances = []
    transfers = []
    transfer_phases = []
    for section in sections:
        for seg in section:

            distance = sim.neuron.h.distance(seg.x, sec=section)

            imp.loc(seg.x, sec=section)
            imp.compute(freq, 1)
            if is_transfer:
                transfer = imp.transfer(0.0, sec=icell.soma[0])
            else:
                transfer = imp.transfer(seg.x, sec=section)

            transfer_phase = imp.transfer_phase(seg.x, sec=section)

            distances.append(distance)
            transfers.append(transfer)
            transfer_phases.append(transfer_phase)
    
    return np.asarray([distances, transfers, transfer_phases])

def plot_trace(ssim, gid):
    plt.plot(ssim.get_time_trace(), ssim.get_voltage_trace(gid))

def run(morph_dir, gid, morph_name, neurite_type, is_transfer=True, tpe='bio', ext='.svg', folder='figure_traces'):

    morph_names = {}
    morph_names[gid] = morph_name

    ssim = initialize_sim(morph_dir)
    ssim.instantiate_gids([gid,], morph_names=morph_names)#, add_synapses=True, add_stimuli=True)#, add_minis=True, add_replay=True)

    t_stop = 20
    ssim.cells[gid].add_voltage_clamp(t_stop, -50, rs=0.01)
    ssim.run(t_stop = t_stop)
    
    plt.figure()
    plot_trace(ssim, gid)
    plt.axis([0,t_stop, -65, 0])
    plt.savefig(os.path.join(folder, 'impedance_' + morph_name + '_' + tpe + ext), bbox_inches= 'tight')
    plt.close()

    icell = ssim.cells[gid].cell

    return compute_impedence(bglibpy, icell, neurite_type=neurite_type, is_transfer=is_transfer)

def initialize_sim(morph_dir):
    bc = '/gpfs/bbp.cscs.ch/project/proj64/circuits/' \
        'S1HL-200um/20171002/simulations/003/BlueConfig'

    ssim = bglibpy.SSim(bc)
    ssim.morph_dir = morph_dir
    
    return ssim

def get_cells_gid(ssim, list_cells, mtypes=None):
    gids = {}
    gids_all = ssim.get_gids_of_mtypes(mtypes)

    for c in list_cells:
        name_cell = os.path.splitext(c)[0]
        for gid in gids_all:
            circuit_name_cell = ssim.fetch_morph_name(gid)
            name_extract = re.search('dend-(.*)_axon',circuit_name_cell)
            if name_extract:
                circuit_name_cell = name_extract.group(1)

            if circuit_name_cell == name_cell:
                if name_cell in gids:
                    gids[name_cell] += [int(gid)]
                else:
                    gids[name_cell] = [int(gid)]

    return gids

def bin_data(bins, data):
    
    mean_data = []
    mean_bins = []
    for i in range(len(bins)-1):
        mean_data.append(np.mean(data[1][( data[0]>bins[i]) & (data[0]< bins[i+1])]))
        mean_bins.append(np.mean([bins[i], bins[i+1]]))

    return np.array(mean_bins), np.array(mean_data)


def plot_impedance(data_bio, data_diam, cell_name, folder, log=False, ext='.svg'):

    f, axes = plt.subplots(2, 1, figsize=(8, 7))

    axes[0].scatter(data_bio[0], data_bio[1], 10, marker='o', label='biological cell', c='C0') #data_bio[2])
    axes[0].scatter(data_diam[0], data_diam[1], 10, marker='o', label='diametrized cell', c='C1') #data_bio[2])
     
    bins, values = set_bins(data_bio[0], 20, 5)

    mean_bins, mean_data_bio = bin_data(bins, data_bio)
    axes[0].plot(mean_bins, mean_data_bio, c='C0', lw=2)

    mean_bins, mean_data_diam = bin_data(bins, data_diam)
    axes[0].plot(mean_bins, mean_data_diam, c='C1', lw=2)
    scale = np.mean(mean_data_diam)
    mean_data_bio /= scale 
    mean_data_diam /= scale

    axes[0].legend()
    axes[0].set_xlabel('distances from soma')
    axes[0].set_ylabel('transfer impedance')

    if log:
        axes[0].set_yscale('log')

    diff_imped = (data_bio[1,:len(data_diam[1])] - data_diam[1,:len(data_bio[1])])/scale

    axes[1].scatter(data_bio[0][:min(len(data_diam[1]),len(data_bio[1]))], diff_imped, 5, c='C2')
    axes[1].plot(mean_bins, mean_data_bio - mean_data_diam, c='C2', lw=2)
    axes[1].axhline(0, c='k')

    axes[1].set_xlabel('distances from soma')
    axes[1].set_ylabel('normalized difference in transfer impedance')

    #error = np.linalg.norm(diff_imped) / len(data_bio[1])
    last = int(len(mean_data_diam)/2)
    error = np.mean( (mean_data_bio-mean_data_diam)[last:])

    #axes[1].axhline(error, c='r')
    #axes[1].axvline(mean_bins[last], c='k')
    axes[0].set_title('normalized discrepancy = ' + str(np.round(error,4)))

    f.savefig(os.path.join(folder,'impedance_' + cell_name + ext), bbox_inches= 'tight')
    plt.close()

    return error 

