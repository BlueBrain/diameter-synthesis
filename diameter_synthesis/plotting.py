"""plotting functions"""
import json
import logging
import os

import pandas
import matplotlib
import matplotlib.pyplot as plt
import neurom as nm
import numpy as np

from neurom import APICAL_DENDRITE, BASAL_DENDRITE, get, iter_sections, viewer
from neurom.geom import bounding_box
from scipy import stats
from tqdm import tqdm

from diameter_synthesis import io, utils
from diameter_synthesis.distribution_fitting import evaluate_distribution
from diameter_synthesis.io import iter_morphology_filenames
from diameter_synthesis.types import STR_TO_TYPES

# pylint: disable=too-many-statements,too-many-locals,too-many-arguments

matplotlib.use("Agg")
L = logging.getLogger(__name__)
COLORS = {"basal": "r", "apical": "m", "axon": "b"}

CUMULATIVE_FEATURE_PAIRS = [
    ("section_radial_distances", "section_volumes"),
    ("section_radial_distances", "section_areas"),
    ("section_path_distances", "section_volumes"),
    ("section_path_distances", "section_areas"),
    ("section_branch_orders", "section_volumes"),
    ("section_branch_orders", "section_areas"),
]

VIOLIN_FEATURES_LIST = [
    "segment_radii",
    "section_areas",
    "section_volumes",
    "sibling_ratio",
    "diameter_power_relation",
]

# to include morphologies only features
VIOLIN_FEATURES_LIST += [
    "number_of_neurites",
    "number_of_sections_per_neurite",
    "number_of_terminations",
    "number_of_bifurcations",
    "section_lengths",
    "section_tortuosity",
    "section_radial_distances",
    "section_path_distances",
    "section_branch_orders",
    "remote_bifurcation_angles",
]

VIOLIN_FEATURES_NAME = [
    "Segment radii",
    "Section areas",
    "Section volumes",
    "Sibling ratios",
    "Diameter power relation",
]

# to include morphologies only features
VIOLIN_FEATURES_NAME += [
    "Number of neurites",
    "Number of sections",
    "Number of terminations",
    "Number of bifurcations",
    "Section lengths",
    "Section tortuosity",
    "Section radial distances",
    "Section path distances",
    "Section branch orders",
    "Remote bif angles",
]


def _compute_neurite_diff(
    neuron_orig, neuron_new, neuron_diff_pos, neuron_diff_neg, neurite_types
):
    """compute the differences between neurite diameters"""
    for neurite_type in neurite_types:
        neurites_orig = [
            neurite
            for neurite in neuron_orig.neurites
            if neurite.type == STR_TO_TYPES[neurite_type]
        ]
        neurites_new = [
            neurite
            for neurite in neuron_new.neurites
            if neurite.type == STR_TO_TYPES[neurite_type]
        ]
        neurites_diff_neg = [
            neurite
            for neurite in neuron_diff_neg.neurites
            if neurite.type == STR_TO_TYPES[neurite_type]
        ]
        neurites_diff_pos = [
            neurite
            for neurite in neuron_diff_pos.neurites
            if neurite.type == STR_TO_TYPES[neurite_type]
        ]

        for neurite_orig, neurite_new, neurite_diff_pos, neurite_diff_neg in zip(
            neurites_orig, neurites_new, neurites_diff_pos, neurites_diff_neg
        ):
            diam_orig = []
            diam_new = []
            for section_orig, section_new in zip(
                iter_sections(neurite_orig), iter_sections(neurite_new)
            ):
                diam_orig.append(utils._get_diameters(section_orig))
                diam_new.append(utils._get_diameters(section_new))

            for j, section in enumerate(iter_sections(neurite_diff_pos)):
                diff = diam_new[j] - diam_orig[j]
                diff_pos = diff.copy()
                diff_pos[diff_pos < 0] = 0
                utils._set_diameters(section, diff_pos)

            for j, section in enumerate(iter_sections(neurite_diff_neg)):
                diff = diam_new[j] - diam_orig[j]
                diff_neg = -diff.copy()
                diff_neg[diff_neg < 0] = 0
                utils._set_diameters(section, diff_neg)


def plot_diameter_diff(
    neuron_name, morph_path, neuron_new, neurite_types, folder, ext="png"
):
    """plot original morphology, new one and differences"""

    if not os.path.isdir(folder):
        os.mkdir(folder)

    neuron_orig = io.load_neuron(neuron_name, morph_path)
    neuron_diff_pos = io.load_neuron(neuron_name, morph_path)
    neuron_diff_neg = io.load_neuron(neuron_name, morph_path)
    _compute_neurite_diff(
        neuron_orig, neuron_new, neuron_diff_pos, neuron_diff_neg, neurite_types
    )

    bbox = bounding_box(neuron_orig)

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    viewer.plot_neuron(axs[0, 0], neuron_orig)
    axs[0, 0].set_title("Original neuron")
    axs[0, 0].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[0, 0].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[0, 0].set_aspect("equal")

    viewer.plot_neuron(axs[0, 1], neuron_new)
    axs[0, 1].set_title("New neuron")
    axs[0, 1].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[0, 1].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[0, 1].set_aspect("equal")

    viewer.plot_neuron(axs[1, 0], neuron_diff_pos)
    axs[1, 0].set_title("Positive diameter differences")
    axs[1, 0].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[1, 0].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[1, 0].set_aspect("equal")

    viewer.plot_neuron(axs[1, 1], neuron_diff_neg)
    axs[1, 1].set_title("Negative diameter differences")
    axs[1, 1].set_xlim([bbox[0, 0], bbox[1, 0]])
    axs[1, 1].set_ylim([bbox[0, 1], bbox[1, 1]])
    axs[1, 1].set_aspect("equal")

    fig.savefig(os.path.join(folder, neuron_new.name + "." + ext), dpi=500)
    plt.close()


def _plot_attribute_scatter(data, model, neurite_types, fig_name, figsize, ext):
    """plot scatter of parameters with attribute if any attributes"""
    plt.figure(figsize=figsize)
    save_plot = False
    max_val = -1e10
    min_val = 1e10
    for neurite_type in neurite_types:
        if model[neurite_type]["sequential"] == "asymmetry_threshold":
            save_plot = True
            tpes = np.asarray(data[neurite_type])[:, 1]
            values = np.asarray(data[neurite_type])[:, 0]
            plt.scatter(values, tpes, s=5, c=COLORS[neurite_type], alpha=0.5)
            plt.axhline(0.2, c="k")

            min_val = min(min_val, model[neurite_type]["params"]["min"])
            max_val = max(max_val, model[neurite_type]["params"]["max"])

    if os.path.basename(fig_name) == "sibling_ratios":
        plt.gca().set_xlim(0.0, 1.0)
    else:
        plt.gca().set_xlim(min_val * 0.5, max_val * 2.0)

    if save_plot:
        plt.savefig(fig_name + "_scatter." + ext)

    plt.close()


def plot_distribution_fit(
    data, model, neurite_types, fig_name="test", ext="png", figsize=(5, 4)
):
    """ Plot the data distribution and its fit """

    _plot_attribute_scatter(data, model, neurite_types, fig_name, figsize, ext)

    plt.figure()
    for neurite_type in neurite_types:

        min_val = model[neurite_type]["params"]["min"]
        max_val = model[neurite_type]["params"]["max"]
        if os.path.basename(fig_name) == "sibling_ratios":
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

    plt.savefig(fig_name + "." + ext, bbox_inches="tight")
    plt.close()


def _create_data(
    feature1, feature2, original_cells, diametrized_cells, step_size, neurite_types
):  # noqa, pylint: disable=too-many-locals,too-many-arguments
    def feature_data(cell, neurite_type):
        nm_neurite_type = STR_TO_TYPES[neurite_type]
        return [
            get(feat, cell, neurite_type=nm_neurite_type)
            for feat in (feature1, feature2)
        ]

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
            yield list(
                create_paired_features(original_cells, diametrized_cells, neurite_type)
            )

    assert len(original_cells) == len(diametrized_cells)

    n_cells = len(original_cells)
    iter_neurite_data = per_neurite_data(
        original_cells, diametrized_cells, neurite_types
    )
    for _, data_pairs in enumerate(iter_neurite_data):

        try:
            upper_bound = find_upper_bound(data_pairs)
        except BaseException:  # pylint: disable=broad-except
            L.exception("failed to find upper bound, most likely due to no data points")
            upper_bound = 200

        bin_centers, bins = create_bins(step_size, upper_bound)

        stats1 = np.empty((n_cells, len(bin_centers)), dtype=np.float)
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
    """
    Plot the cumulative distribution of feature2 with respect to
    the metric values determined via feature1

    Args:
        original_cells: list of NeuroM objects
        diametrized_cells: list of NeuroM objects
            The new cell with the changed diameters.
        neurite_types: string
            The list neurite types to be considered. e.g. ['basal', 'axon']
        step_size: float
            The step size of the cumulative histogram

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
        axes[0].plot(
            bin_centers, means, c=color, linestyle="-", lw=3, label="original cells"
        )

        if len(stats1) > 1:
            for st1 in stats1:
                axes[0].plot(
                    bin_centers, st1, c=color, linestyle="-", lw=0.5, alpha=0.2
                )

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
                axes[0].plot(
                    bin_centers, st2, c=color, linestyle="-", lw=0.3, alpha=0.5
                )

            sdevs = stats2.std(axis=0)
            axes[0].plot(bin_centers, means - sdevs, c=color, linestyle="--", lw=3)
            axes[0].plot(bin_centers, means + sdevs, c=color, linestyle="--", lw=3)

        axes[0].set_xlabel("path distance")
        axes[0].set_ylabel("cummulative section areas")

        axes[0].set_xlabel("path distance")
        axes[0].set_ylabel("cummulative section areas")

        stats1[stats1 == 0] = 1
        diffs = stats1 - stats2  # / stats1

        diff_means = diffs.mean(axis=0)

        color = "C2"
        if len(diffs) > 1:
            for dfs in diffs:
                axes[1].plot(
                    bin_centers, dfs, c=color, linestyle="-", lw=0.3, alpha=0.5
                )

            diff_sdevs = diffs.std(axis=0)
            axes[1].plot(
                bin_centers, diff_means - diff_sdevs, c=color, linestyle="--", lw=3
            )
            axes[1].plot(
                bin_centers, diff_means + diff_sdevs, c=color, linestyle="--", lw=3
            )

        axes[1].plot(bin_centers, diff_means, c=color, linestyle="-", lw=3)
        axes[1].axhline(0, ls="--", c="k")

        axes[1].set_xlabel("path distance")
        axes[1].set_ylabel("difference in cummulative section areas")

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
            "L2 error = "
            + str(np.around(np.linalg.norm(stats2[:, -1] - stats1[:, -1]), 1)),
            loc="left",
        )

        axes[2].set_xlabel("total surface area of original cells")
        axes[2].set_ylabel("total surface area of diametrized cells")

    return fig, axes


def _split_prefix(neurom_feature_name):

    name_list = neurom_feature_name.split("_")

    prefix = name_list[0]
    basename = "_".join(name_list[1:])

    return prefix, basename


def make_cumulative_figures(
    original_cells,
    diametrized_cells,
    feature1,
    feature2,
    config,
    out_dir,
    individual=False,
    figname_prefix="",
):

    """make plots for cumulative distributions for a pair of features"""
    prefix1, basename1 = _split_prefix(feature1)
    prefix2, basename2 = _split_prefix(feature2)

    assert prefix1 == prefix2

    fig, _ = plot_cumulative_distribution(
        original_cells, diametrized_cells, feature1, feature2, config["neurite_types"]
    )

    figure_name = figname_prefix + "cumulative_{}_{}_{}".format(
        prefix1, basename1, basename2
    )

    fig.savefig(os.path.join(out_dir, figure_name + ".svg"), bbox_inches="tight")
    plt.close(fig)

    if individual:
        if not os.path.isdir(os.path.join(out_dir, figure_name + "_individual")):
            os.mkdir(os.path.join(out_dir, figure_name + "_individual"))

        for i, (original_cell, diametrized_cell) in enumerate(
            zip(original_cells, diametrized_cells)
        ):
            f, _ = plot_cumulative_distribution(
                [original_cell],
                [diametrized_cell],
                feature1,
                feature2,
                config["neurite_types"],
                auto_limit=False,
            )
            fname = "{}_{}.svg".format(figure_name, original_cell.name)
            f.savefig(
                os.path.join(
                    out_dir, figure_name + "_individual/", str(i) + "_" + fname
                ),
                bbox_inches="tight",
            )
            plt.close(f)


def _load_cells(config):
    filenames_original = sorted(list(iter_morphology_filenames(config["morph_path"])))
    filenames_diametrized = sorted(
        list(iter_morphology_filenames(config["new_morph_path"]))
    )

    original_filepaths = (
        os.path.join(config["morph_path"], filename) for filename in filenames_original
    )
    diametrized_filepaths = (
        os.path.join(config["new_morph_path"], filename)
        for filename in filenames_diametrized
    )

    original_cells = list(map(nm.load_neuron, original_filepaths))
    diametrized_cells = list(map(nm.load_neuron, diametrized_filepaths))

    return original_cells, diametrized_cells


def _load_morphologies(
    morph_path, mtypes_sort="all", mtypes_file="./neuronDB.xml", ext=".asc", prefix="",
):
    """ Load the morphologies from a directory, by mtypes or all at once """
    name_dict = utils.create_morphologies_dict(
        morph_path,
        mtypes_sort=mtypes_sort,
        mtypes_file=mtypes_file,
        ext=ext,
        prefix=prefix,
    )
    return io.load_morphologies_from_dict(morph_path, name_dict)


def cumulative_analysis(config, out_dir, individual):
    """make plots for cumulative distributions"""
    with open(config, "r") as filename:
        config = json.load(filename)

    # out_dir += "_" + config["mtypes_sort"]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if len(config) > 1:
        L.warning(
            "multiple models provided, will only use the first in the list for analysis"
        )
    config = config[list(config.keys())[0]]

    all_original_cells = _load_morphologies(
        config["morph_path"], mtypes_sort=config["mtypes_sort"]
    )

    all_diametrized_cells = _load_morphologies(
        config["new_morph_path"], mtypes_sort=config["mtypes_sort"]
    )

    for mtype in tqdm(all_original_cells):
        original_cells = all_original_cells[mtype]
        diametrized_cells = all_diametrized_cells[mtype]
        try:
            for feature1, feature2 in CUMULATIVE_FEATURE_PAIRS:
                make_cumulative_figures(
                    original_cells,
                    diametrized_cells,
                    feature1,
                    feature2,
                    config,
                    out_dir,
                    individual=individual,
                    figname_prefix=mtype,
                )
        except BaseException:  # pylint: disable=broad-except
            pass


def get_features_all(object1, object2, flist, neurite_type):
    """Computes features from module mod"""
    collect_all = []

    for feat in flist:
        feature_pop = []
        for obj in object1:
            feature_pop = (
                feature_pop + get(feat, obj, neurite_type=neurite_type).tolist()
            )
        feature_neu = []
        for obj in object2:
            feature_neu = (
                feature_neu + get(feat, obj, neurite_type=neurite_type).tolist()
            )

        collect_all.append([feature_pop, feature_neu])

    return collect_all


def transform2DataFrame(data, pop_names, flist):
    """Returns a DataFrame in the appropriate
    format from a set of features"""
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

    frame = pandas.DataFrame(
        {"Data": names, "Values": values, "Morphological features": feat}
    )

    return frame


def plot_violins(data, x="Morphological features", y="Values", hues="Data", **kwargs):
    """Plots the split violins of all features"""
    import seaborn  # pylint: disable=import-outside-toplevel
    plt.figure(figsize=(12, 6))
    axs = seaborn.violinplot(
        x=x,
        y=y,
        hue=hues,
        data=data,
        palette="muted",
        split=True,
        inner="quartile",
        **kwargs
    )
    plt.xticks(rotation=20)
    plt.tight_layout(True)
    return axs


def violin_analysis(config, out_dir):
    """plot violin distributions"""
    with open(config, "r") as filename:
        config = json.load(filename)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if len(config) > 1:
        L.warning(
            "multiple models provided, will only use the first in the list for analysis"
        )
    config = config[list(config.keys())[0]]
    all_original_cells = _load_morphologies(
        config["morph_path"], mtypes_sort=config["mtypes_sort"]
    )

    all_diametrized_cells = _load_morphologies(
        config["new_morph_path"], mtypes_sort=config["mtypes_sort"]
    )

    for mtype in tqdm(all_original_cells):
        original_cells = all_original_cells[mtype]
        diametrized_cells = all_diametrized_cells[mtype]

        pop_names = ["Original cells of " + mtype, "Diametrized cells of " + mtype]
        data = get_features_all(
            original_cells,
            diametrized_cells,
            flist=VIOLIN_FEATURES_LIST,
            neurite_type=BASAL_DENDRITE,
        )

        try:
            data_frame = transform2DataFrame(
                data, pop_names, flist=VIOLIN_FEATURES_NAME
            )
            ax = plot_violins(data_frame)
            ax.set_ylim(-3, 5)
            plt.savefig(os.path.join(out_dir, "violin_basal_" + mtype + ".png"))
            plt.close()
        except BaseException:  # pylint: disable=broad-except
            pass

        try:
            data = get_features_all(
                original_cells,
                diametrized_cells,
                flist=VIOLIN_FEATURES_LIST,
                neurite_type=APICAL_DENDRITE,
            )
            data_frame = transform2DataFrame(
                data, pop_names, flist=VIOLIN_FEATURES_NAME
            )
            ax = plot_violins(data_frame)
            ax.set_ylim(-3, 5)
            plt.savefig(os.path.join(out_dir, "violin_apical_" + mtype + ".png"))
            plt.close()
        except BaseException:  # pylint: disable=broad-except
            pass
