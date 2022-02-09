"""Plotting functions."""

# Copyright (C) 2021  Blue Brain Project, EPFL
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import logging
import multiprocessing
import os
import re
from copy import copy
from functools import partial
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import neurom as nm
import numpy as np
import pandas
from matplotlib.backends.backend_pdf import PdfPages
from neurom import APICAL_DENDRITE
from neurom import AXON
from neurom import BASAL_DENDRITE
from neurom import COLS
from neurom import get
from neurom import iter_sections
from neurom.geom import bounding_box
from neurom.view import matplotlib_impl
from scipy import stats
from tqdm import tqdm

from diameter_synthesis import utils
from diameter_synthesis.distribution_fitting import evaluate_distribution

# pylint: disable=too-many-statements,too-many-locals,too-many-arguments

matplotlib.use("Agg")
L = logging.getLogger(__name__)
COLORS = {"basal": "r", "apical": "m", "axon": "b"}

CUMULATIVE_FEATURE_PAIRS = [
    ("section_path_distances", "section_volumes"),
    ("section_path_distances", "section_areas"),
]

VIOLIN_FEATURES_LIST = ["segment_radii", "section_areas", "section_volumes"]

VIOLIN_FEATURES_NAME = ["Segment radii", "Section areas", "Section volumes"]

NEURITE_STR_TO_TYPES = {
    "basal": nm.NeuriteType.basal_dendrite,
    "apical": nm.NeuriteType.apical_dendrite,
    "axon": nm.NeuriteType.axon,
}


def _compute_neurite_diff(neuron_orig, neuron_new, neuron_diff_pos, neuron_diff_neg, neurite_types):
    """Compute the differences between neurite diameters."""
    for neurite_type in neurite_types:
        neurites_orig = [
            neurite
            for neurite in neuron_orig.neurites
            if neurite.type == NEURITE_STR_TO_TYPES[neurite_type]
        ]
        neurites_new = [
            neurite
            for neurite in neuron_new.neurites
            if neurite.type == NEURITE_STR_TO_TYPES[neurite_type]
        ]
        neurites_diff_neg = [
            neurite
            for neurite in neuron_diff_neg.neurites
            if neurite.type == NEURITE_STR_TO_TYPES[neurite_type]
        ]
        neurites_diff_pos = [
            neurite
            for neurite in neuron_diff_pos.neurites
            if neurite.type == NEURITE_STR_TO_TYPES[neurite_type]
        ]

        for neurite_orig, neurite_new, neurite_diff_pos, neurite_diff_neg in zip(
            neurites_orig, neurites_new, neurites_diff_pos, neurites_diff_neg
        ):
            orig = []
            new = []
            for section_orig, section_new in zip(
                iter_sections(neurite_orig), iter_sections(neurite_new)
            ):
                orig.append(section_orig.points)
                new.append(section_new.points)

            for j, section in enumerate(iter_sections(neurite_diff_pos)):
                _diff = new[j][:, COLS.R] - orig[j][:, COLS.R]
                _diff[_diff < 0] = 0
                diff = copy(orig[j])
                diff[:, COLS.R] = _diff
                section.points = diff

            for j, section in enumerate(iter_sections(neurite_diff_neg)):
                _diff = -(new[j][:, COLS.R] - orig[j][:, COLS.R])
                _diff[_diff < 0] = 0
                diff = copy(orig[j])
                diff[:, COLS.R] = _diff
                section.points = diff


def plot_diameter_diff(neuron_name, neuron_new, neurite_types, folder, ext=".png"):
    """Plot original morphology, new one and differences."""
    if not Path(folder).exists():
        os.mkdir(folder)

    neuron_orig = nm.load_morphology(neuron_name)
    neuron_diff_pos = nm.load_morphology(neuron_name)
    neuron_diff_neg = nm.load_morphology(neuron_name)
    _compute_neurite_diff(neuron_orig, neuron_new, neuron_diff_pos, neuron_diff_neg, neurite_types)

    bbox = bounding_box(neuron_orig)

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    matplotlib_impl.plot_morph(neuron_orig, axs[0, 0])
    axs[0, 0].set_title("Original neuron")
    axs[0, 0].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[0, 0].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[0, 0].set_aspect("equal")

    matplotlib_impl.plot_morph(neuron_new, axs[0, 1])
    axs[0, 1].set_title("New neuron")
    axs[0, 1].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[0, 1].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[0, 1].set_aspect("equal")

    matplotlib_impl.plot_morph(neuron_diff_pos, axs[1, 0])
    axs[1, 0].set_title("Positive diameter differences")
    axs[1, 0].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[1, 0].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[1, 0].set_aspect("equal")

    matplotlib_impl.plot_morph(neuron_diff_neg, axs[1, 1])
    axs[1, 1].set_title("Negative diameter differences")
    axs[1, 1].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[1, 1].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[1, 1].set_aspect("equal")

    fig.savefig((Path(folder) / neuron_name.stem).with_suffix(ext), dpi=500)
    plt.close()


def _plot_attribute_scatter(data, model, neurite_types, fig_name, figsize, ext):
    """Plot scatter of parameters with attribute if any attributes."""
    plt.figure(figsize=figsize)
    save_plot = False
    max_val = -1e10
    min_val = 1e10
    for neurite_type in neurite_types:
        if len(model[neurite_type]) > 0:
            if model[neurite_type]["sequential"] == "asymmetry_threshold":
                save_plot = True
                tpes = np.asarray(data[neurite_type])[:, 1]
                values = np.asarray(data[neurite_type])[:, 0]
                plt.scatter(values, tpes, s=5, c=COLORS[neurite_type], alpha=0.5)
                plt.axhline(0.2, c="k")

                min_val = min(min_val, model[neurite_type]["params"]["min"])
                max_val = max(max_val, model[neurite_type]["params"]["max"])

    if Path(fig_name).name == "sibling_ratios":
        plt.gca().set_xlim(0.0, 1.0)
    else:
        plt.gca().set_xlim(min_val * 0.5, max_val * 2.0)

    if save_plot:
        plt.savefig(str(fig_name) + "_scatter" + ext)

    plt.close()


def plot_distribution_fit(data, model, neurite_types, fig_name="test", ext=".png", figsize=(5, 4)):
    """Plot the data distribution and its fit."""
    _plot_attribute_scatter(data, model, neurite_types, fig_name, figsize, ext)

    plt.figure()
    for neurite_type in neurite_types:
        if len(model[neurite_type]) > 0:
            min_val = model[neurite_type]["params"]["min"]
            max_val = model[neurite_type]["params"]["max"]
            if Path(fig_name).name in ["sibling_ratios", "tapers"]:
                hist_range = [min_val, max_val]
            else:
                hist_range = [min_val * 0.5, max_val * 2.0]

            if len(data[neurite_type]) > 0:
                if len(np.shape(data[neurite_type])) > 1:
                    data_hist = np.array(data[neurite_type])[:, 0]
                else:
                    data_hist = data[neurite_type]

                plt.hist(
                    data_hist,
                    bins=30,
                    density=True,
                    color=COLORS[neurite_type],
                    alpha=0.5,
                    range=hist_range,
                )

                values = np.linspace(min_val, max_val, 1000)
                plt.plot(
                    values,
                    evaluate_distribution(
                        values,
                        model[neurite_type]["distribution"],
                        model[neurite_type]["params"],
                    ),
                    c=COLORS[neurite_type],
                    lw=3,
                    ls="--",
                    label=neurite_type,
                )
                plt.legend(loc="best")
            else:
                L.warning("No data to plot")

    plt.savefig(str(fig_name) + ext, bbox_inches="tight")
    plt.close()


def _create_data(
    feature1, feature2, original_cells, diametrized_cells, step_size, neurite_types
):  # noqa, pylint: disable=too-many-locals,too-many-arguments
    def feature_data(cell, neurite_type):
        nm_neurite_type = NEURITE_STR_TO_TYPES[neurite_type]
        return [get(feat, cell, neurite_type=nm_neurite_type) for feat in (feature1, feature2)]

    def create_paired_features(cell_list1, cell_list2, neurite_type):
        for cell1, cell2 in zip(cell_list1, cell_list2):
            yield feature_data(cell1, neurite_type), feature_data(cell2, neurite_type)

    def create_bins(step_size, max_value):
        bins = np.arange(0.0, max_value + step_size, step_size)
        bin_centers = 0.5 * step_size + bins[:-1]
        return bin_centers, bins

    def find_upper_bound(pairs):
        return max(max(max(vals1), max(vals2)) for (vals1, _), (vals2, _) in pairs)

    def per_neurite_data(original_cells, diametrized_cells, neurite_types):
        for _, neurite_type in enumerate(neurite_types):
            data = list(create_paired_features(original_cells, diametrized_cells, neurite_type))
            if len(data[0][0][0]) > 0:
                yield data

    assert len(original_cells) == len(diametrized_cells)

    n_cells = len(original_cells)
    iter_neurite_data = per_neurite_data(original_cells, diametrized_cells, neurite_types)
    for _, data_pairs in enumerate(iter_neurite_data):

        try:
            upper_bound = find_upper_bound(data_pairs)
        except BaseException:  # pylint: disable=broad-except
            L.exception("failed to find upper bound, most likely due to no data points")
            upper_bound = 200

        bin_centers, bins = create_bins(step_size, upper_bound)

        stats1 = np.empty((n_cells, len(bin_centers)), dtype=float)
        stats2 = np.empty_like(stats1)

        for i, ((metric1, data1), (metric2, data2)) in enumerate(data_pairs):

            res1 = stats.binned_statistic(
                metric1, data1, statistic="sum", bins=bins, range=(0, upper_bound)
            )
            res2 = stats.binned_statistic(
                metric2, data2, statistic="sum", bins=bins, range=(0, upper_bound)
            )

            stats1[i] = np.cumsum(res1.statistic)
            stats2[i] = np.cumsum(res2.statistic)

        yield bin_centers, stats1, stats2


def plot_cumulative_distribution(
    original_cells,
    diametrized_cells,
    feature1,
    feature2,
    neurite_types,
    step_size=1.0,
    auto_limit=True,
):
    """Plot the cumulative distribution of features.

    It plots feature2 with respect to the metric values determined via feature1.

    Args:
        original_cells: list of NeuroM objects.
        diametrized_cells (list): The new cells with the changed diameters.
        feature1: the metric feature.
        feature2: the cumulative distribution feature.
        neurite_types (string): The list neurite types to be considered. e.g. ['basal', 'axon'].
        step_size (float): The step size of the cumulative histogram.
        auto_limit (bool): automatically compute limits.

    Examples of metric features (feature1):
        - segment_radial_distances
        - segment_path_distances (not implemented yet)

    Examples of cumulative distribution features (feature2):
        - segment_volumes
        - segment_surface_areas (not implemented yet)
    """
    assert len(original_cells) == len(diametrized_cells)

    data_generator = _create_data(
        feature1, feature2, original_cells, diametrized_cells, step_size, neurite_types
    )

    fig, axes = plt.subplots(3, 1, figsize=(5, 12))

    for (bin_centers, stats1, stats2) in data_generator:

        means = stats1.mean(axis=0)
        color = "C0"
        axes[0].plot(bin_centers, means, c=color, linestyle="-", lw=3, label="original cells")

        if len(stats1) > 1:
            for st1 in stats1:
                axes[0].plot(bin_centers, st1, c=color, linestyle="-", lw=0.5, alpha=0.2)

            sdevs = stats1.std(axis=0)
            axes[0].plot(bin_centers, means - sdevs, c=color, linestyle="--", lw=3)
            axes[0].plot(bin_centers, means + sdevs, c=color, linestyle="--", lw=3)

        means = stats2.mean(axis=0)

        color = "C1"
        axes[0].plot(
            bin_centers,
            means,
            c=color,
            linestyle="-",
            lw=3,
            label="rediametrized cells",
        )
        axes[0].legend(loc="best")
        if len(stats2) > 1:
            for st2 in stats2:
                axes[0].plot(bin_centers, st2, c=color, linestyle="-", lw=0.3, alpha=0.5)

            sdevs = stats2.std(axis=0)
            axes[0].plot(bin_centers, means - sdevs, c=color, linestyle="--", lw=3)
            axes[0].plot(bin_centers, means + sdevs, c=color, linestyle="--", lw=3)

        axes[0].set_xlabel("path distance")
        axes[0].set_ylabel("cumulative section areas")

        axes[0].set_xlabel("path distance")
        axes[0].set_ylabel("cumulative section areas")

        stats1[stats1 == 0] = 1
        diffs = stats1 - stats2  # / stats1

        diff_means = diffs.mean(axis=0)

        color = "C2"
        if len(diffs) > 1:
            for dfs in diffs:
                axes[1].plot(bin_centers, dfs, c=color, linestyle="-", lw=0.3, alpha=0.5)

            diff_sdevs = diffs.std(axis=0)
            axes[1].plot(bin_centers, diff_means - diff_sdevs, c=color, linestyle="--", lw=3)
            axes[1].plot(bin_centers, diff_means + diff_sdevs, c=color, linestyle="--", lw=3)

        axes[1].plot(bin_centers, diff_means, c=color, linestyle="-", lw=3)
        axes[1].axhline(0, ls="--", c="k")

        axes[1].set_xlabel("path distance")
        axes[1].set_ylabel("difference in cumulative section areas")

        if auto_limit:
            lim_min = 0.5 * np.min(stats1[:, -1])
            lim_max = np.max(stats1[:, -1])
        else:
            lim_min = 5000
            lim_max = 30000

        axes[2].scatter(stats1[:, -1], stats2[:, -1], c=color, marker="o", s=5)

        var_x = np.arange(lim_min, lim_max)

        axes[2].plot(var_x, var_x, ls="-", c="k")
        axes[2].set_xlim(lim_min, lim_max)
        axes[2].set_ylim(lim_min, lim_max)
        axes[2].set_title(
            "L2 error = " + str(np.around(np.linalg.norm(stats2[:, -1] - stats1[:, -1]), 1)),
            loc="left",
        )

        axes[2].set_xlabel("total surface area of original cells")
        axes[2].set_ylabel("total surface area of diametrized cells")

    return fig, axes


def _split_prefix(neurom_feature_name):
    """Split prefix."""
    name_list = neurom_feature_name.split("_")

    prefix = name_list[0]
    basename = "_".join(name_list[1:])

    return prefix, basename


def make_cumulative_figures(
    original_cells,
    diametrized_cells,
    feature1,
    feature2,
    neurite_types,
    out_dir,
    individual=False,
    figname_prefix="",
    ext=".png",
):
    """Make plots for cumulative distributions for a pair of features."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix1, basename1 = _split_prefix(feature1)
    prefix2, basename2 = _split_prefix(feature2)

    assert prefix1 == prefix2

    fig, _ = plot_cumulative_distribution(
        original_cells, diametrized_cells, feature1, feature2, neurite_types
    )

    if figname_prefix and not figname_prefix.endswith("_"):
        figname_prefix = figname_prefix + "_"
    figname_prefix = re.sub(r"[^\w\s-]", "_", figname_prefix)

    figure_name = f"{figname_prefix}cumulative_{prefix1}_{basename1}_{basename2}"
    fig.savefig(out_dir / f"{figure_name}{ext}", bbox_inches="tight")
    plt.close(fig)

    if individual:
        if not (out_dir / (figure_name + "_individual")).exists():
            os.mkdir(out_dir / (figure_name + "_individual"))

        for i, (original_cell, diametrized_cell) in enumerate(
            zip(original_cells, diametrized_cells)
        ):
            f, _ = plot_cumulative_distribution(
                [original_cell],
                [diametrized_cell],
                feature1,
                feature2,
                neurite_types,
                auto_limit=False,
            )
            fname = f"{figure_name}_{original_cell.name}{ext}"
            f.savefig(
                out_dir / (figure_name + "_individual") / (str(i) + "_" + fname),
                bbox_inches="tight",
            )
            plt.close(f)


def _load_morphologies(morph_path, mtypes_file="./neuronDB.xml"):
    """Load the morphologies from a directory, by mtypes or all at once."""
    morphologies_dict = utils.create_morphologies_dict(morph_path, mtypes_file=mtypes_file)
    return {
        mtype: [nm.load_morphology(i) for i in morphologies_dict[mtype]]
        for mtype in morphologies_dict
    }


def cumulative_analysis(
    original_path,
    diametrized_path,
    out_dir,
    individual=False,
    mtypes_file=None,
    neurite_types=None,
    ext=".png",
):
    """Make plots for cumulative distributions."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_original_cells = _load_morphologies(original_path, mtypes_file=mtypes_file)
    all_diametrized_cells = _load_morphologies(diametrized_path, mtypes_file=mtypes_file)
    for mtype in tqdm(all_original_cells):
        original_cells = all_original_cells[mtype]
        diametrized_cells = all_diametrized_cells[mtype]
        for feature1, feature2 in CUMULATIVE_FEATURE_PAIRS:
            make_cumulative_figures(
                original_cells,
                diametrized_cells,
                feature1,
                feature2,
                neurite_types,
                out_dir,
                individual=individual,
                figname_prefix=mtype,
                ext=ext,
            )


def get_features_all(object1, object2, flist, neurite_type):
    """Compute features from module mod."""

    def _get_feature(feat, obj):
        v = get(feat, obj, neurite_type=neurite_type)
        return v if isinstance(v, list) else [v]

    collect_all = []
    for feat in flist:
        feature_pop = []
        for obj in object1:
            feature_pop = feature_pop + _get_feature(feat, obj)
        feature_neu = []
        for obj in object2:
            feature_neu = feature_neu + _get_feature(feat, obj)

        collect_all.append([feature_pop, feature_neu])
    return collect_all


def transform2DataFrame(data, pop_names, flist):
    """Return a DataFrame in the appropriate format from a set of features."""
    values = []
    names = []
    feat = []
    # Loop through each feature
    for i, d1 in enumerate(data):
        m = np.mean(d1[0])
        st = np.std(d1[0])
        # Loop through populations
        for j, d2 in enumerate(d1):
            values = values + [(d3 - m) / st for d3 in d2]
            names = names + len(d2) * [pop_names[j]]
            feat = feat + len(d2) * [flist[i]]

    frame = pandas.DataFrame({"Data": names, "Values": values, "Morphological features": feat})

    return frame


def plot_violins(data, x="Morphological features", y="Values", hues="Data", ax=None):
    """Plot the split violins of all features."""
    import seaborn  # pylint: disable=import-outside-toplevel

    if ax is None:
        fig = plt.figure(figsize=(12, 6))
        ax = plt.gca()
    else:
        fig = None

    seaborn.violinplot(
        x=x,
        y=y,
        hue=hues,
        data=data,
        palette="muted",
        split=True,
        inner="quartile",
        ax=ax,
    )
    ax.tick_params(axis="x", rotation=20)
    return fig, ax


def violin_analysis(
    original_path,
    diametrized_path,
    out_dir,
    mtypes_file=None,
    max_cells=200,
    with_axon=False,
):
    """Plot violin distributions."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    orig_morphologies_dict = utils.create_morphologies_dict(original_path, mtypes_file=mtypes_file)
    diametrized_morphologies_dict = utils.create_morphologies_dict(
        diametrized_path, mtypes_file=mtypes_file
    )

    cells_data = [
        [orig_morphologies_dict[mtype], diametrized_morphologies_dict[mtype], mtype]
        for mtype in orig_morphologies_dict
    ]
    analyze_from_dict = partial(_analyze_from_dict, max_cells, with_axon=with_axon)

    pool = multiprocessing.Pool()  # pylint: disable=consider-using-with
    try:
        figs = list(
            tqdm(
                pool.imap_unordered(analyze_from_dict, cells_data),
                total=len(cells_data),
            )
        )
    finally:
        pool.close()
        pool.join()

    with PdfPages(out_dir / "morphometrics.pdf") as pdf:
        for mtype, fig in figs:
            if mtype is not None:
                pdf.savefig(fig)


def _analyze_from_dict(max_cells, cells, with_axon=False):
    cell_orig, cell_diametrized, mtype = cells
    cell_diametrized = cell_diametrized[:max_cells]
    cell_orig = cell_orig[:max_cells]
    original_cells = [nm.load_morphology(i) for i in cell_orig]
    diametrized_cells = [nm.load_morphology(i) for i in cell_diametrized]

    pop_names = ["Original cells of " + mtype, "Synthetised cells of " + mtype]
    data = get_features_all(
        original_cells,
        diametrized_cells,
        flist=VIOLIN_FEATURES_LIST,
        neurite_type=BASAL_DENDRITE,
    )

    data_frame = transform2DataFrame(data, pop_names, flist=VIOLIN_FEATURES_NAME)

    if with_axon:
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(15, 10))
    else:
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(15, 10))
    fig.subplots_adjust(hspace=0.5)
    plot_violins(data_frame.replace([np.inf, -np.inf], np.nan).dropna(), ax=axes[0])
    axes[0].set_ylim(-3, 5)
    axes[0].title.set_text("basal dendrites")

    if any(i.type == nm.NeuriteType.apical_dendrite for i in original_cells[0].neurites) and any(
        i.type == nm.NeuriteType.apical_dendrite for i in diametrized_cells[0].neurites
    ):
        data = get_features_all(
            original_cells,
            diametrized_cells,
            flist=VIOLIN_FEATURES_LIST,
            neurite_type=APICAL_DENDRITE,
        )

        data_frame = transform2DataFrame(data, pop_names, flist=VIOLIN_FEATURES_NAME)
        plot_violins(data_frame.replace([np.inf, -np.inf], np.nan).dropna(), ax=axes[1])
        axes[1].set_ylim(-3, 5)
        axes[1].title.set_text("apical dendrites")

    if with_axon:
        data = get_features_all(
            original_cells,
            diametrized_cells,
            flist=VIOLIN_FEATURES_LIST,
            neurite_type=AXON,
        )
        data_frame = transform2DataFrame(data, pop_names, flist=VIOLIN_FEATURES_NAME)
        plot_violins(data_frame.replace([np.inf, -np.inf], np.nan).dropna(), ax=axes[2])
        axes[2].set_ylim(-3, 5)

    return mtype, fig
