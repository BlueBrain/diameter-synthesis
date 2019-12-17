import os, shutil, re

#orig_dir = '../../examples/PC_example/new_morphologies_basal/'
orig_dir ='/gpfs/bbp.cscs.ch/project/proj82/home/arnaudon/2017_model/2017_PC_cells/new_morphologies_mtypes/'
model = 'M1'
for name in os.listdir(orig_dir):
    morph = re.search(model + '_(.*)', name)
    if morph:
        shutil.copy(os.path.join(orig_dir, name), os.path.join('new_neurons', morph.group(1)))



