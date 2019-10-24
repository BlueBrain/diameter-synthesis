import os, glob, shutil
import json

########################################
## This scripts selects morphologies  ##
## that can be used in diameter-check ##
########################################

def find_etypes_emodels(morph, emodel_morph_map, emodel_etype_map):

    etypes_emodels = []
    for emodel, emodel_morph in emodel_morph_map.items():
        if morph == emodel_morph:
            etypes_emodels.append((emodel_etype_map[emodel], emodel))

    if len(etypes_emodels) == 0:
        return False
    else:
        return True


def read_emodel_etype_map(emodel_etype_map_path):
    with open(emodel_etype_map_path) as emodel_etype_map_file:
        emodel_etype = json.load(emodel_etype_map_file)

    emodel_etype_map = {}
    for emodel, emodel_info in emodel_etype.items():
        emodel_etype_map[emodel] = emodel_info['etype']

    return emodel_etype_map

def read_final_json(final_json_path):
    with open(final_json_path) as final_json_file:
        final_json = json.load(final_json_file)

    emodel_morph_map = {}
    for emodel, emodel_info in final_json.items():
        morph = os.path.splitext(
            os.path.basename(
                emodel_info['morph_path']))[0]
        emodel_morph_map[emodel] = morph

    return emodel_morph_map


def extract_test_morphologies(config_file, final_json, emodel_etype_map):

    with open(config_file) as f:
        config = json.load(f)
    

    emodel_morph_map = read_final_json(final_json)
    emodel_etype_map = read_emodel_etype_map(emodel_etype_map)

    if not os.path.isdir(config['morph_path']):
        os.mkdir(config['morph_path'])

    #copy the neuronDB.xml file to the new folder for later analysis
    shutil.copy(config['rep_morph_path']+'neuronDB.xml', config['morph_path']+'neuronDB.xml')

    fnames = []
    n = 0 
    for fname in os.listdir(config['rep_morph_path']):
        name, ext = os.path.splitext(fname)
        if ext in {'.h5', '.asc', '.swc'} and find_etypes_emodels(os.path.splitext(fname)[0], emodel_morph_map, emodel_etype_map): #check if the etype exists for this neuron
            shutil.copyfile(config['rep_morph_path']+fname, config['morph_path']+fname)
            n+=1
        else:
            print('no corresponding emodel for '+ name)
    print('Found ' + str(n) + ' matching morphologies.')


if __name__ == "__main__":
    config_file = 'config.json'
    final_json = 'final.json'
    emodel_etype_map = 'emodel_etype_map.json'

    extract_test_morphologies(config_file, final_json, emodel_etype_map)
