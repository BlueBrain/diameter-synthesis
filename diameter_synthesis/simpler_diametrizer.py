"""Simpler diametrizer."""
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from neurom import NeuriteType
from neurom import iter_sections
from neurom import morphmath
from neurom.core.morphology import Morphology
from neurom.core.morphology import Section
from numpy.polynomial import Polynomial
from scipy.stats import pearsonr


def _s_length(s_id, points, cache):
    """Section length with cache."""
    if s_id not in cache:
        cache[s_id] = morphmath.section_length(points)
    return cache[s_id]


def section_path_length(section, cache):
    """Path length from section to root."""
    return sum(_s_length(s.id, s.points, cache) for s in section.iupstream())


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
                        diams.append(np.mean(section.points[:, 3]))
                        tip_length = max(
                            section_path_length(_section, cache) for _section in section.ipreorder()
                        )
                        lengths.append(
                            tip_length
                            - section_path_length(section, cache)
                            + _s_length(section.id, section.points, cache)
                        )
            if lengths:
                lengths = np.array(lengths)
                diams = np.array(diams)
                lengths /= lengths.max()
                all_lengths[neurite_type] += lengths.tolist()
                all_diams[neurite_type] += diams.tolist()
        if all_lengths[neurite_type]:
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


def _update_diameters(section, diameters):
    """Shortcut to update diameters with neurom v3."""
    points = section.points
    points[:, 3] = diameters
    section.points = points


# pylint: disable=unused-argument
def simpler_diametrizer(morphology, coeffs, neurite_types, config=None, rng=np.random):
    """Diametrize a morphology."""
    morphology = Morphology(morphology)

    for neurite in morphology.neurites:
        if neurite.type.name in coeffs:
            coeff = coeffs[neurite.type.name]
            p = Polynomial(coeff)
            cache = {}
            max_len = max(terminal_path_lengths(neurite, cache))
            for section in iter_sections(neurite):
                tip_length = max(
                    section_path_length(_section, cache) for _section in section.ipreorder()
                )
                section_len = (
                    tip_length
                    - section_path_length(section, cache)
                    + _s_length(section.id, section.points, cache)
                )
                diam = p(section_len / max_len)
                _update_diameters(section, len(section.points) * [diam])
            for section in iter_sections(neurite):
                diam_start = section.points[0, 3]
                if section.children:
                    diam_end = max(sec.points[0, 3] for sec in section.children)
                    diam_end = 0.5 * (diam_start + diam_end)
                else:
                    diam_end = p(0)
                _update_diameters(
                    section,
                    diam_start + (diam_end - diam_start) * np.linspace(0, 1, len(section.points)),
                )
    return morphology
