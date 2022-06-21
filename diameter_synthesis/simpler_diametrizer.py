"""Simpler diametrizer."""
from scipy.stats import pearsonr
from tqdm import tqdm
import matplotlib.pyplot as plt

import numpy as np
from neurom import load_morphology, iter_sections, NeuriteType

from numpy.polynomial import Polynomial
from matplotlib.backends.backend_pdf import PdfPages
from neurom import morphmath
from neurom.core.morphology import Section
from functools import partial


def s_length(s_id, points, cache):
    if s_id not in cache:
        cache[s_id] = morphmath.section_length(points)
    return cache[s_id]


def section_path_length(section, cache):
    """Path length from section to root."""
    return sum(s_length(s.id, s.points, cache) for s in section.iupstream())


def _map_sections(fun, neurite, iterator_type=Section.ipreorder):
    """Map `fun` to all the sections."""
    return list(map(fun, iterator_type(neurite.root_node)))


def terminal_path_lengths(neurite, cache):
    """Get the path lengths to each terminal point."""
    return _map_sections(partial(section_path_length, cache=cache), neurite, Section.ileaf)


def build_model(morphologies, config, neurite_types=None, fit_orders=None):
    neurite_types = config["neurite_types"]
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
                            + s_length(section.id, section.points, cache)
                        )
            if len(lengths):
                lengths = np.array(lengths)
                diams = np.array(diams)
                lengths /= lengths.max()
                all_lengths[neurite_type] += list(lengths)
                all_diams[neurite_type] += list(diams)
        if len(all_lengths[neurite_type]):
            p, extra = Polynomial.fit(
                all_lengths[neurite_type],
                all_diams[neurite_type],
                fit_orders[neurite_type],
                full=True,
            )
            residues[neurite_type] = extra[0][0]
            coeffs[neurite_type] = p.convert().coef.tolist()
    return coeffs, [all_lengths, all_diams, residues]


def plot_model(coeffs, pdf, title_str, all_lengths, all_diams, residues):
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


def make_model(df, neurite_types=None, fit_orders=None):
    mtypes = sorted(df.mtype.unique())
    if neurite_types is None:
        neurite_types = ["basal_dendrite", "apical_dendrite"]
    if fit_orders is None:
        fit_orders = {"basal_dendrite": 1, "apical_dendrite": 2, "axon": 1}
    coeffs = {mtype: {} for mtype in mtypes}
    residues = {mtype: {} for mtype in mtypes}
    all_diams = {mtype: {n_type: [] for n_type in neurite_types} for mtype in mtypes}
    all_lengths = {mtype: {n_type: [] for n_type in neurite_types} for mtype in mtypes}

    for neurite_type in neurite_types:
        for mtype in tqdm(mtypes):
            _df = df[df.mtype == mtype]
            for gid in _df.index:
                diams = []
                lengths = []
                m = load_morphology(df.loc[gid, "path"])
                for neurite in m.neurites:
                    cache = {}
                    if neurite.type == getattr(NeuriteType, neurite_type):
                        for section in iter_sections(neurite):
                            diams.append(np.mean(section.points[:, 3]))
                            tip_length = max(
                                section_path_length(_section, cache)
                                for _section in section.ipreorder()
                            )
                            lengths.append(
                                tip_length
                                - section_path_length(section, cache)
                                + s_length(section.id, section.points, cache)
                            )
                if len(lengths):
                    lengths = np.array(lengths)
                    diams = np.array(diams)
                    lengths /= lengths.max()
                    all_lengths[mtype][neurite_type] += list(lengths)
                    all_diams[mtype][neurite_type] += list(diams)
            if len(all_lengths[mtype][neurite_type]):
                p, extra = Polynomial.fit(
                    all_lengths[mtype][neurite_type],
                    all_diams[mtype][neurite_type],
                    fit_orders[neurite_type],
                    full=True,
                )
                residues[mtype][neurite_type] = extra[0][0]
                coeffs[mtype][neurite_type] = p.convert().coef.tolist()
    return coeffs, all_lengths, all_diams, residues


def plot_fit(coeffs, all_lengths, all_diams, residues, figure_name="diameter_prof.pdf"):
    color = {"basal_dendrite": "r", "apical_dendrite": "m", "axon": "b"}
    with PdfPages(figure_name) as pdf:
        for mtype, _coeffs in coeffs.items():
            plt.figure(figsize=(5, 4))
            title_str = f"""mtype: {mtype}"""
            for neurite_type, coeff in _coeffs.items():
                xd = all_lengths[mtype][neurite_type]
                yd = all_diams[mtype][neurite_type]
                plt.scatter(xd, yd, s=0.5, marker=".", c=color[neurite_type])
                x = np.linspace(0, 1, 100)
                p = Polynomial(coeff)
                plt.plot(x, p(x), c=color[neurite_type], lw=2, label=neurite_type)
                res = np.round(residues[mtype][neurite_type], 1)
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


def simpler_diametrizer(morphology, coeffs, neurite_types, config, rng=np.random):
    """Diametrize a morphology."""
    from neurom.core.morphology import Morphology
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
                    + s_length(section.id, section.points, cache)
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

def diametrize(morphology, coeffs):
    """Diametrize a morphology."""
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
                    + s_length(section.id, section.points, cache)
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
