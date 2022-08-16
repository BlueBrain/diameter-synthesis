"""Example to diametrize single morphologie."""
import matplotlib.pyplot as plt
import numpy as np
from morphio.mut import Morphology
from neurom import NeuriteType
from neurom import load_morphology
from neurom import view

from diameter_synthesis.main import diametrize_single_neuron

if __name__ == "__main__":
    morph = Morphology("PATH_TO_MORPHOLOGY")
    smooth_morph = diametrize_single_neuron(morph)
    morph = load_morphology(morph)
    new_morph = load_morphology(smooth_morph)

    plt.figure()
    ax = plt.gca()
    view.plot_morph(
        morph.transform(lambda p: p - np.array([500, 0, 0])),
        ax=ax,
        plane="xy",
        soma_outline=False,
        realistic_diameters=False,
        neurite_type=NeuriteType.apical_dendrite,
    )
    view.plot_morph(
        new_morph.transform(lambda p: p + np.array([500, 0, 0])),
        ax=ax,
        plane="xy",
        soma_outline=False,
        realistic_diameters=False,
        neurite_type=NeuriteType.apical_dendrite,
    )
    view.plot_morph(
        morph.transform(lambda p: p - np.array([500, 0, 0])),
        ax=ax,
        plane="xy",
        soma_outline=False,
        realistic_diameters=False,
        neurite_type=NeuriteType.basal_dendrite,
    )
    view.plot_morph(
        new_morph.transform(lambda p: p + np.array([500, 0, 0])),
        ax=ax,
        plane="xy",
        soma_outline=False,
        realistic_diameters=False,
        neurite_type=NeuriteType.basal_dendrite,
    )
    plt.autoscale()
    plt.axis("equal")
    plt.savefig("compare_rediametetrized.pdf", bbox_inches="tight")
