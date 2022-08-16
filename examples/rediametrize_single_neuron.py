"""Example to diametrize single morphologie."""
from diameter_synthesis.main import diametrize_single_neuron
from morphio.mut import Morphology
import matplotlib.pyplot as plt
from neurom import view, NeuriteType
from neurom import load_morphology
import numpy as np

if __name__ == "__main__":
    morph = Morphology("C060110A2.asc")
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
