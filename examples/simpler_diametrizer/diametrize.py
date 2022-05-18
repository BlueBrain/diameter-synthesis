from morph_tool.morphdb import MorphDB
import yaml
from tqdm import tqdm
import numpy as np
from neurom import load_morphology, NeuriteType
from neurom import view

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from diameter_synthesis.simpler_diametrizer import make_model, plot_fit, diametrize

if __name__ == "__main__":
    path = "/gpfs/bbp.cscs.ch/project/proj83/home/gevaert/morph-release/morph_release_old_code-2020-07-27/output/06_RepairUnravel-asc/"
    df = MorphDB.from_neurondb(path + "neuronDB.xml", morphology_folder=path).df

    coeffs, all_lengths, all_diams, residues = make_model(df)
    yaml.dump(coeffs, open("diameter_fits.yaml", "w"))

    plot_fit(coeffs, all_lengths, all_diams, residues)

    with PdfPages(f"diametrized.pdf") as pdf:
        for mtype in tqdm(sorted(df.mtype.unique())):
            _df = df[df.mtype == mtype]
            for gid in _df.index[:1]:
                m = load_morphology(df.loc[gid, "path"])

                plt.figure()
                view.plot_morph(m, neurite_type=NeuriteType.basal_dendrite, ax=plt.gca())
                view.plot_morph(m, neurite_type=NeuriteType.apical_dendrite, ax=plt.gca())
                diametrize(m, coeffs[mtype])
                _m = m.transform(lambda x: x + np.array([500, 0, 0]))
                view.plot_morph(_m, neurite_type=NeuriteType.basal_dendrite, ax=plt.gca())
                view.plot_morph(_m, neurite_type=NeuriteType.apical_dendrite, ax=plt.gca())
                plt.axis("equal")
                plt.tight_layout()
                plt.gca().set_title(mtype)
                pdf.savefig()
                plt.close()
                m.write("diametrized_morphs/" + df.loc[gid, "name"] + ".asc")
