""" plotting functions """
import collections
import json
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas
from neurom import APICAL_DENDRITE, BASAL_DENDRITE, get, iter_sections, viewer
from neurom.geom import bounding_box
from scipy import stats
import seaborn

import diameter_synthesis.utils as utils
from diameter_synthesis import io
from diameter_synthesis.distribution_fitting import (evaluate_distribution,
                                                     evaluate_spline)
from diameter_synthesis.io import iter_morphology_filenames, load_morphology
from diameter_synthesis.types import STR_TO_TYPES

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

feat_list = [
    "segment_radii",
    "section_areas",
    "section_volumes",
]

feat_names = [
    "Segment radii",
    "Section areas",
    "Section volumes",
]


######################
# plotting functions #
######################


def plot_fit_distribution_params(
    model, neurite_types, fig_name="test", ext=".png"
):  # pylint: disable=too-many-locals

    """ plot the fit of the parameter of the distributions """

    fig = plt.figure(figsize=(5, 8))

    fig.subplots_adjust(hspace=0.0)
    try:
        ax1 = fig.add_subplot(511)
        for neurite_type in neurite_types:
            tpes_model = list(model[neurite_type]["params"]["params_data"])

            var_x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
            ax1.plot(
                var_x,
                evaluate_spline(var_x, model[neurite_type]["params"]["a"]),
                c=COLORS[neurite_type],
            )

            var_as = np.array(
                [v["a"] for v in model[neurite_type]["params"]["params_data"].values()]
            )
            var_w = np.array(
                [
                    v["num_value"]
                    for v in model[neurite_type]["params"]["params_data"].values()
                ]
            )

            # prevent large values of a from bad fitting
            var_as[var_as > utils.A_MAX] = utils.A_MAX
            ax1.scatter(
                tpes_model,
                var_as,
                s=var_w,
                edgecolors=COLORS[neurite_type],
                facecolors="none",
            )

        ax1.set_xlabel("max path distance")
        ax1.set_ylabel("a")
        n_plot = 0
    except BaseException:  # pylint: disable=broad-except
        n_plot = 101

    ax2 = fig.add_subplot(512 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = list(model[neurite_type]["params"]["params_data"])

        var_x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax2.plot(
            var_x,
            evaluate_spline(var_x, model[neurite_type]["params"]["loc"]),
            c=COLORS[neurite_type],
        )

        locs = [v["loc"] for v in model[neurite_type]["params"]["params_data"].values()]
        var_w = np.array(
            [
                v["num_value"]
                for v in model[neurite_type]["params"]["params_data"].values()
            ]
        )

        ax2.scatter(
            tpes_model,
            locs,
            s=var_w,
            edgecolors=COLORS[neurite_type],
            facecolors="none",
        )
    ax2.set_ylabel("loc")

    ax3 = fig.add_subplot(513 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = list(model[neurite_type]["params"]["params_data"])

        var_x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax3.plot(
            var_x,
            evaluate_spline(var_x, model[neurite_type]["params"]["scale"]),
            c=COLORS[neurite_type],
        )

        scales = [
            v["scale"] for v in model[neurite_type]["params"]["params_data"].values()
        ]
        ax3.scatter(
            tpes_model,
            scales,
            s=var_w,
            edgecolors=COLORS[neurite_type],
            facecolors="none",
        )
    ax3.set_ylabel("scale")

    ax4 = fig.add_subplot(514 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = list(model[neurite_type]["params"]["params_data"])

        var_x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax4.plot(
            var_x,
            evaluate_spline(var_x, model[neurite_type]["params"]["min"]),
            c=COLORS[neurite_type],
        )

        mins = [v["min"] for v in model[neurite_type]["params"]["params_data"].values()]
        ax4.scatter(
            tpes_model,
            mins,
            s=var_w,
            edgecolors=COLORS[neurite_type],
            facecolors="none",
        )
    ax4.set_ylabel("min")

    ax5 = fig.add_subplot(515 - n_plot)
    for neurite_type in neurite_types:
        tpes_model = list(model[neurite_type]["params"]["params_data"])

        var_x = np.linspace(tpes_model[0], tpes_model[-1], 1000)
        ax5.plot(
            var_x,
            evaluate_spline(var_x, model[neurite_type]["params"]["max"]),
            c=COLORS[neurite_type],
        )

        maxs = [v["max"] for v in model[neurite_type]["params"]["params_data"].values()]
        ax5.scatter(
            tpes_model,
            maxs,
            s=var_w,
            edgecolors=COLORS[neurite_type],
            facecolors="none",
        )

    ax5.set_xlabel("max branching order")
    ax5.set_ylabel("max")

    fig.savefig(fig_name + ext, bbox_inches="tight")
    plt.close(fig)


def plot_distribution_fit(  # noqa, pylint: disable=too-many-locals,too-many-arguments,too-many-branches,too-many-statements
    data, model, neurite_types, fig_name="test", ext=".png", figsize=(5, 4), n_bins=10
):
    """ Plot the data distribution and its fit """

    plt.figure()
    save_plot = False
    for neurite_type in neurite_types:
        if model[neurite_type]["sequential"] == "asymmetry_threshold":
            tpes = np.asarray(data[neurite_type])[:, 1]  # collect the type of point
            values = np.asarray(data[neurite_type])[:, 0]  # collect the type of point

            plt.scatter(values, tpes, s=5, c=COLORS[neurite_type], alpha=0.5)
            plt.axhline(0.2, c="k")
            save_plot = True
    if save_plot:
        plt.savefig(fig_name + "_scatter" + ext)
    plt.close()

    fig = plt.figure(figsize=figsize)

    for neurite_type in neurite_types:
        # if the fits are done sequantially
        if (
            not isinstance(model[neurite_type]["sequential"], str)
            or model[neurite_type]["sequential"] == "asymmetry_threshold"
        ):

            min_val = model[neurite_type]["params"]["min"]
            max_val = model[neurite_type]["params"]["max"]

            if len(data[neurite_type]) > 0:
                if len(np.shape(data[neurite_type])) > 1:
                    plt.hist(
                        np.array(data[neurite_type])[:, 0],
                        bins=50,
                        log=False,
                        density=True,
                        histtype="bar",
                        lw=0.5,
                        color=COLORS[neurite_type],
                        alpha=0.5,
                        range=[min_val * 0.8, max_val * 1.2],
                    )
                else:
                    plt.hist(
                        data[neurite_type],
                        bins=50,
                        log=False,
                        density=True,
                        histtype="bar",
                        lw=0.5,
                        color=COLORS[neurite_type],
                        alpha=0.5,
                        range=[min_val * 0.8, max_val * 1.2],
                    )

                var_x = np.linspace(min_val, max_val, 1000)
                plt.plot(
                    var_x,
                    evaluate_distribution(
                        var_x,
                        model[neurite_type]["distribution"],
                        model[neurite_type]["params"],
                    ),
                    c=COLORS[neurite_type],
                    lw=3,
                    ls="--",
                    label=neurite_type,
                )

            plt.legend(loc="best")
            plt.gca().set_xlim(min_val * 0.8, max_val * 1.2)

        else:

            tpes = np.asarray(data[neurite_type])[:, 1]  # collect the type of point
            values = np.asarray(data[neurite_type])[:, 0]  # collect the data itself

            bins, _ = utils.set_bins(tpes, n_bins)

            min_val = 1e10
            max_val = -1e10

            for i in range(len(bins) - 1):
                values_tpe = values[
                    (tpes >= bins[i]) & (tpes < bins[i + 1])
                ]  # select the values by its type
                tpe_mean = (bins[i] + bins[i + 1]) / 2.0

                bottom_shift = tpe_mean

                height = (bins[i + 1] - bins[i]) / 2.0

                plt.axhline(bottom_shift, lw=0.2, c="k")
                try:  # try to plot, may not be a fit to plot
                    min_tpe = evaluate_spline(
                        tpe_mean, model[neurite_type]["params"]["min"]
                    )
                    max_tpe = evaluate_spline(
                        tpe_mean, model[neurite_type]["params"]["max"]
                    )

                    min_val = np.min([min_val, min_tpe])
                    max_val = np.max([max_val, max_tpe])
                    values_tpe = values_tpe[values_tpe < max_tpe * 1.2]
                    values_tpe = values_tpe[values_tpe > min_tpe * 0.8]

                    var_n, var_b = np.histogram(values_tpe, bins=20)
                    plt.bar(
                        var_b[1:],
                        height * np.array(var_n) / np.max(var_n),
                        width=var_b[1] - var_b[0],
                        bottom=bottom_shift,
                        color=COLORS[neurite_type],
                        alpha=0.5,
                    )

                    var_x = np.linspace(min_tpe, max_tpe, 1000)
                    params = {}
                    try:
                        params["a"] = evaluate_spline(
                            tpe_mean, model[neurite_type]["params"]["a"]
                        )
                    except BaseException:  # pylint: disable=broad-except
                        pass
                    params["loc"] = evaluate_spline(
                        tpe_mean, model[neurite_type]["params"]["loc"]
                    )
                    params["scale"] = evaluate_spline(
                        tpe_mean, model[neurite_type]["params"]["scale"]
                    )
                    pdf = evaluate_distribution(
                        var_x, model[neurite_type]["distribution"], params
                    )

                    plt.plot(
                        var_x,
                        bottom_shift + height * pdf / np.max(pdf),
                        c=COLORS[neurite_type],
                        lw=1,
                        ls="--",
                    )
                except Exception as exc:  # pylint: disable=broad-except
                    L.exception("plotting exception: %s", exc)
                    var_n, var_b = np.histogram(values_tpe, bins=20)
                    plt.bar(
                        var_b[:-1],
                        height * np.array(var_n) / np.max(var_n),
                        width=var_b[1] - var_b[0],
                        bottom=bottom_shift,
                        color=COLORS[neurite_type],
                        alpha=0.5,
                    )

    title_txt = "Fit parameters:\n"
    try:
        for neurite_type in neurite_types:
            params_title = {}
            for param in model[neurite_type]["params"]:
                if param != "params_data":
                    params_title[param] = model[neurite_type]["params"][param]
            title_txt += neurite_type + ": " + str(params_title)
            title_txt += "\n"
    except BaseException:  # pylint: disable=broad-except
        title_txt += " no parameter could be fitted."

    # plt.title(title_txt, fontsize = 8)

    plt.savefig(fig_name + ext, bbox_inches="tight")
    plt.close(fig)

    if (
        isinstance(
            # if the fits are done sequantially
            model[neurite_type]["sequential"],
            str,
        )
        and model[neurite_type]["sequential"] != "asymmetry_threshold"
    ):
        plot_fit_distribution_params(
            model, neurite_types, fig_name=fig_name + "_param_fit", ext=ext
        )


def plot_fit_param_boxes(  # pylint: disable=too-many-locals,too-many-arguments,too-many-branches
    model_params,
    model="M0",
    neurite_type="basal",
    figname="test",
    ext=".png",
    figsize=(6, 3),
):
    """ box plots for the fits of the model parameters """

    data = collections.OrderedDict()

    mtype = list(model_params)[0]
    for fit in model_params[mtype][model]:
        if fit != "trunk_diameter":
            for params in model_params[mtype][model][fit][neurite_type]["params"]:
                data[fit + "_" + params] = []
        else:
            for params in model_params[mtype][model][fit][neurite_type]["params"]:
                data[fit + "_" + params + "_0"] = []
                data[fit + "_" + params + "_1"] = []

    for mtype in model_params:
        for fit in model_params[mtype][model]:
            if fit != "trunk_diameter":
                for params in model_params[mtype][model][fit][neurite_type]["params"]:
                    data[fit + "_" + params].append(
                        model_params[mtype][model][fit][neurite_type]["params"][params]
                    )
            else:
                for params in model_params[mtype][model][fit][neurite_type]["params"]:
                    params_0 = model_params[mtype][model][fit][neurite_type]["params"][
                        params
                    ][0]
                    params_1 = model_params[mtype][model][fit][neurite_type]["params"][
                        params
                    ][1]
                    data[fit + "_" + params + "_0"].append(params_0)
                    data[fit + "_" + params + "_1"].append(params_1)

    plt.figure(figsize=figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data) + 1), list(data), rotation="vertical")
    # plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + ext, bbox_inches="tight")
    plt.close()

    data = collections.OrderedDict()

    bos = [
        "0.0",
        "1.0",
        "2.0",
        "3.0",
        "4.0",
        "5.0",
        "6.0",
        "7.0",
        "8.0",
    ]
    mtype = list(model_params)[0]
    for fit in model_params[mtype][model]:
        if fit == "trunk_diameter":
            for params in model_params[mtype][model][fit][neurite_type]["params_data"][
                "0.0"
            ]:
                for branch_order in bos:
                    data[branch_order + "_" + params] = []

    for mtype in model_params:
        fit_gen = (fit for fit in model_params[mtype][model] if fit == "trunk_diameter")
        for _ in fit_gen:
            params_tmp = model_params[mtype][model][fit][neurite_type]["params_data"]
            bo_gen = (
                bo
                for bo in model_params[mtype][model][fit][neurite_type]["params_data"]
                if bo in bos
            )
            for branch_order in bo_gen:
                bo_list = model_params[mtype][model][fit][neurite_type]["params_data"][
                    branch_order
                ]
                for params in bo_list:
                    data[branch_order + "_" + params].append(params_tmp)

    plt.figure(figsize=figsize)
    plt.boxplot(data.values())
    plt.xticks(np.arange(1, len(data) + 1), list(data), rotation="vertical")
    # plt.axis([0, len(data)+1, 0., 5])
    plt.savefig(figname + "_trunk" + ext, bbox_inches="tight")
    plt.close()


def plot_diameter_diff(  # pylint: disable=too-many-locals,too-many-arguments
    neuron_name, morph_path, neuron_new, model, neurite_types, folder, ext=".png"
):
    """ plot original morphology, new one and differences """

    if not os.path.isdir(folder):
        os.mkdir(folder)

    neuron_orig = io.load_neuron(neuron_name, None, morph_path)
    neuron_diff_pos = io.load_neuron(neuron_name, None, morph_path)
    neuron_diff_neg = io.load_neuron(neuron_name, None, morph_path)

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

    fig.savefig(
        folder + "/" + model + "_" + neuron_new.name + "_" + folder + "_" + model + ext,
        dpi=500,
    )
    plt.close("all")


def _compute_neurite_diff(
    neuron_orig, neuron_new, neuron_diff_pos, neuron_diff_neg, neurite_types
):  # pylint: disable=too-many-locals
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
                diam_orig.append(utils.get_diameters(section_orig))
                diam_new.append(utils.get_diameters(section_new))

            for j, section in enumerate(iter_sections(neurite_diff_pos)):
                diff = diam_new[j] - diam_orig[j]
                diff_pos = diff.copy()
                diff_pos[diff_pos < 0] = 0
                utils.set_diameters(section, diff_pos)

            for j, section in enumerate(iter_sections(neurite_diff_neg)):
                diff = diam_new[j] - diam_orig[j]
                diff_neg = -diff.copy()
                diff_neg[diff_neg < 0] = 0
                utils.set_diameters(section, diff_neg)


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
            upper_bound = 2000

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


def plot_cumulative_distribution(  # noqa, pylint: disable=too-many-statements,too-many-locals,too-many-arguments
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

        if len(stats1) > 1:  # don't plot std if there is only one curve
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
            # axes[1].fill_between(bin_centers, diff_means - diff_sdevs,
            # diff_means + diff_sdevs, color=color, alpha=0.2)
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

        # if len(diffs) > 1:
        #    for i in range(len(stats1)):
        #        axes[2].annotate(i, (stats1[i,-1], stats2[i,-1]))

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

        # if len(diffs) > 1:
        #    axes[2].plot(x,x-diff_means[-1], ls='--', c=color)
        #    axes[2].plot(x,x+diff_means[-1], ls='--', c=color)

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
):  # pylint: disable=too-many-locals

    """make plots for cumulative distributions for a pair of features"""
    prefix1, basename1 = _split_prefix(feature1)
    prefix2, basename2 = _split_prefix(feature2)

    assert prefix1 == prefix2

    fig, _ = plot_cumulative_distribution(
        original_cells, diametrized_cells, feature1, feature2, config["neurite_types"]
    )

    figure_name = "cumulative_{}_{}_{}".format(prefix1, basename1, basename2)

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

    original_cells = list(map(load_morphology, original_filepaths))
    diametrized_cells = list(map(load_morphology, diametrized_filepaths))

    return original_cells, diametrized_cells


def cumulative_analysis(config, out_dir, individual):
    """make plots for cumulative distributions"""
    with open(config, "r") as filename:
        config = json.load(filename)

    out_dir += "_" + config["mtypes_sort"]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if len(config["models"]) > 1:
        L.warning(
            "multiple models provided, will only use the first in the list for analysis"
        )

    original_cells, diametrized_cells = _load_cells(config)

    for feature1, feature2 in CUMULATIVE_FEATURE_PAIRS:
        make_cumulative_figures(
            original_cells,
            diametrized_cells,
            feature1,
            feature2,
            config,
            out_dir,
            individual=individual,
        )


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

    out_dir += "_" + config["mtypes_sort"]
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if len(config["models"]) > 1:
        L.warning(
            "multiple models provided, will only use the first in the list for analysis"
        )

    original_cells, diametrized_cells = _load_cells(config)

    pop_names = ["Original cells", "Diametrized cells"]
    data = get_features_all(
        original_cells, diametrized_cells, flist=feat_list, neurite_type=BASAL_DENDRITE
    )
    data_frame = transform2DataFrame(data, pop_names, flist=feat_list)
    ax = plot_violins(data_frame)
    ax.set_ylim(-3, 5)
    plt.savefig(os.path.join(out_dir, "violin_basal.png"))

    data = get_features_all(
        original_cells, diametrized_cells, flist=feat_list, neurite_type=APICAL_DENDRITE
    )
    data_frame = transform2DataFrame(data, pop_names, flist=feat_list)
    ax = plot_violins(data_frame)
    ax.set_ylim(-3, 5)
    plt.savefig(os.path.join(out_dir, "violin_apical.png"))
