"""Simpler diametrizer."""
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from neurom import NeuriteType
from neurom import iter_sections
from neurom.core.morphology import Morphology
from neurom.core.morphology import Section
from numpy.polynomial import Polynomial
from scipy.stats import pearsonr


def section_path_length(section, cache):
    """Path length from section to root."""
    if section.id not in cache:
        cache[section.id] = sum(s.length for s in section.iupstream())
    return cache[section.id]


def _map_sections(fun, neurite, iterator_type=Section.ipreorder):
    """Map `fun` to all the sections."""
    return list(map(fun, iterator_type(neurite.root_node)))


def terminal_path_lengths(neurite, cache):
    """Get the path lengths to each terminal point."""
    return _map_sections(partial(section_path_length, cache=cache), neurite, Section.ileaf)


def build_simpler_model(morphologies, config, fit_orders=None):
    """Build diameter model."""
    neurite_types = config.get("neurite_types")
    if neurite_types is None:
        neurite_types = ["basal_dendrite", "apical_dendrite"]
    if fit_orders is None:
        fit_orders = {"basal_dendrite": 1, "apical_dendrite": 2, "axon": 1}

    coeffs = {}
    residues = {}
    all_diams = {n_type: [] for n_type in neurite_types}
    all_lengths = {n_type: [] for n_type in neurite_types}

    for neurite_type in neurite_types:
        for m in morphologies:
            diams = []
            lengths = []
            for neurite in m.neurites:
                cache = {}
                if neurite.type == getattr(NeuriteType, neurite_type):
                    for section in iter_sections(neurite):
                        diams.append(2 * np.mean(section.points[:, 3]))
                        tip_length = max(
                            section_path_length(_section, cache) for _section in section.ipreorder()
                        )
                        lengths.append(
                            tip_length - section_path_length(section, cache) + section.length
                        )
            lengths = np.array(lengths)
            diams = np.array(diams)
            lengths /= lengths.max()
            all_lengths[neurite_type] += lengths.tolist()
            all_diams[neurite_type] += diams.tolist()
        p, extra = Polynomial.fit(
            all_lengths[neurite_type],
            all_diams[neurite_type],
            fit_orders[neurite_type],
            full=True,
        )
        residues[neurite_type] = float(extra[0][0])
        coeffs[neurite_type] = p.convert().coef.tolist()
    return coeffs, [all_lengths, all_diams, residues]


def plot_model(coeffs, pdf, title_str, all_lengths, all_diams, residues):
    """Plot diameter model."""
    color = {"basal_dendrite": "r", "apical_dendrite": "m", "axon": "b"}
    plt.figure(figsize=(5, 4))
    for neurite_type, coeff in coeffs.items():
        xd = all_lengths[neurite_type]
        yd = all_diams[neurite_type]
        plt.scatter(xd, yd, s=0.5, marker=".", c=color[neurite_type])
        x = np.linspace(0, 1, 100)
        p = Polynomial(coeff)
        plt.plot(x, p(x), c=color[neurite_type], lw=2, label=neurite_type)
        res = np.round(residues[neurite_type], 1)
        coef = [np.round(c, 3) for c in coeff]
        pears = np.around(pearsonr(xd, yd)[0], 2)
        title_str += f"""
        {neurite_type}: fit: {coef}, residual: {res}, pearson: {pears}"""
    plt.suptitle(title_str, fontsize=5)
    plt.xlabel("normalise path distance to most distal tip")
    plt.ylabel("mean section radius")
    plt.legend()
    pdf.savefig()
    plt.close()


# pylint: disable=unused-argument
def simpler_diametrizer(morphology, neurite_types, model_params, diam_params=None, rng=np.random):
    """Diametrize a morphology."""
    if not isinstance(neurite_types, list):
        neurite_types = [neurite_types]
    neurite_types = [
        getattr(NeuriteType, neurite_type) if isinstance(neurite_type, str) else neurite_type
        for neurite_type in neurite_types
    ]
    _morphology = Morphology(morphology)
    for neurite in _morphology.neurites:
        if neurite.type in neurite_types and neurite.type.name in model_params:
            p = Polynomial(model_params[neurite.type.name])
            cache = {}
            max_len = max(terminal_path_lengths(neurite, cache))
            for section in iter_sections(neurite):
                tip_length = max(
                    section_path_length(_section, cache) for _section in section.ipreorder()
                )
                section_length = tip_length - section_path_length(section, cache) + section.length
                diam = p(section_length / max_len)
                morphology.sections[section.id].diameters = len(section.points) * [diam]
            for section in iter_sections(neurite):
                diam_start = morphology.sections[section.id].diameters[0]
                if section.children:
                    diam_end = max(
                        morphology.sections[sec.id].diameters[0] for sec in section.children
                    )
                    diam_end = 0.5 * (diam_start + diam_end)
                else:
                    diam_end = p(0)
                morphology.sections[section.id].diameters = diam_start + (
                    diam_end - diam_start
                ) * np.linspace(0, 1, len(section.points))
