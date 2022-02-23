"""Test the plotting module."""

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

import neurom as nm
import pytest
from diff_pdf_visually import pdfdiff
from matplotlib import pyplot as plt

from diameter_synthesis import plotting
from diameter_synthesis.exception import DiameterSynthesisError


def test_plot_diameter_diff(neuron_diametrized_path, tmpdir, expected_images_path):
    """Test the plot_diameter_diff function."""
    # Test with existing directory and with positive changes

    # Load neuron with NeuroM
    new_neuron = nm.load_morphology(neuron_diametrized_path)

    # Update the neuron
    for sec in new_neuron.sections:
        sec.points[:, nm.COLS.R] *= 20

    # Plot the figure
    plotting.plot_diameter_diff(
        neuron_diametrized_path, new_neuron, ["basal", "apical"], tmpdir, ext=".pdf"
    )

    # Check the figure
    assert pdfdiff(
        str(expected_images_path / "test_plot_diameter_diff_mult.pdf"),
        str(tmpdir / neuron_diametrized_path.with_suffix(".pdf").name),
    )

    # Test with new directory and with negative changes

    # Load neuron with NeuroM
    new_neuron = nm.load_morphology(neuron_diametrized_path)

    # Update the neuron
    for sec in new_neuron.sections:
        sec.points[:, nm.COLS.R] /= 20

    # Plot the figure
    plotting.plot_diameter_diff(
        neuron_diametrized_path, new_neuron, ["basal", "apical"], tmpdir / "new_dir", ext=".pdf"
    )

    # Check the figure
    assert pdfdiff(
        str(expected_images_path / "test_plot_diameter_diff_div.pdf"),
        str(tmpdir / "new_dir" / neuron_diametrized_path.with_suffix(".pdf").name),
    )


def test_plot_distribution_fit(model_params, model_data, tmpdir, expected_images_path):
    """Test the plot_distribution_fit function."""
    # Plot the figure
    plotting.plot_distribution_fit(
        model_data["diameter_power_relation"],
        model_params["generic"]["L5_TPC:A"]["diameter_power_relation"],
        ["basal", "apical"],
        fig_name=tmpdir / "test_plot_distribution_fit",
        ext=".pdf",
    )

    # Check the figures
    assert pdfdiff(
        str(expected_images_path / "test_plot_distribution_fit.pdf"),
        str(tmpdir / "test_plot_distribution_fit.pdf"),
    )
    assert pdfdiff(
        str(expected_images_path / "test_plot_distribution_fit_scatter.pdf"),
        str(tmpdir / "test_plot_distribution_fit_scatter.pdf"),
    )


def test_plot_cumulative_distribution(
    single_pop, single_pop_diametrized, tmpdir, expected_images_path
):
    """Test the plot_cumulative_distribution function."""
    # Plot the figure
    plotting.plot_cumulative_distribution(
        single_pop,
        single_pop_diametrized,
        "segment_radial_distances",
        "segment_volumes",
        ["basal", "apical"],
    )
    plt.savefig(tmpdir / "test_plot_cumulative_distribution.pdf", bbox_inches="tight")
    plt.close()

    # Check the figures
    assert pdfdiff(
        str(expected_images_path / "test_plot_cumulative_distribution.pdf"),
        str(tmpdir / "test_plot_cumulative_distribution.pdf"),
    )


def test_make_cumulative_figures(single_pop, single_pop_diametrized, tmpdir, expected_images_path):
    """Test the make_cumulative_figures function."""
    # Plot the figures
    plotting.make_cumulative_figures(
        single_pop,
        single_pop_diametrized,
        "segment_radial_distances",
        "segment_volumes",
        ["basal", "apical"],
        tmpdir,
        individual=True,
        figname_prefix="with_individual_",
        ext=".pdf",
    )

    # Check the figures
    images_path = expected_images_path / "test_make_cumulative_figures"
    assert pdfdiff(
        str(images_path / "cumulative_segment_radial_distances_volumes.pdf"),
        str(tmpdir / "with_individual_cumulative_segment_radial_distances_volumes.pdf"),
    )
    assert pdfdiff(
        str(
            images_path
            / "cumulative_segment_radial_distances_volumes_individual"
            / "0_cumulative_segment_radial_distances_volumes.pdf"
        ),
        str(
            tmpdir
            / "with_individual_cumulative_segment_radial_distances_volumes_individual"
            / (
                "0_with_individual_cumulative_segment_radial_distances_volumes_C030796A-P3_lite"
                ".h5.pdf"
            )
        ),
    )


def test_cumulative_analysis(
    single_pop_neurondb,
    tmpdir,
    expected_images_path,
    single_pop_data_dir,
    single_pop_diametrized_data_dir,
):
    """Test the cumulative_analysis function."""
    # Plot the figures
    plotting.cumulative_analysis(
        single_pop_data_dir,
        single_pop_diametrized_data_dir,
        tmpdir / "analysis",
        True,
        single_pop_data_dir / "neurondb.dat",
        ["basal", "apical"],
        ext=".pdf",
    )

    # Check the figures
    images_path = expected_images_path / "test_cumulative_analysis" / "analysis"
    assert pdfdiff(
        str(images_path / "L5_TPC_A_cumulative_section_path_distances_areas.pdf"),
        str(tmpdir / "analysis" / "L5_TPC_A_cumulative_section_path_distances_areas.pdf"),
    )
    assert pdfdiff(
        str(images_path / "L5_TPC_A_cumulative_section_path_distances_volumes.pdf"),
        str(tmpdir / "analysis" / "L5_TPC_A_cumulative_section_path_distances_volumes.pdf"),
    )
    assert pdfdiff(
        str(
            images_path
            / "L5_TPC_A_cumulative_section_path_distances_areas_individual"
            / "0_L5_TPC_A_cumulative_section_path_distances_areas_C030796A-P3_lite.pdf"
        ),
        str(
            tmpdir
            / "analysis"
            / "L5_TPC_A_cumulative_section_path_distances_areas_individual"
            / "0_L5_TPC_A_cumulative_section_path_distances_areas_C030796A-P3_lite.h5.pdf"
        ),
    )
    assert pdfdiff(
        str(
            images_path
            / "L5_TPC_A_cumulative_section_path_distances_volumes_individual"
            / "0_L5_TPC_A_cumulative_section_path_distances_volumes_C030796A-P3_lite.pdf"
        ),
        str(
            tmpdir
            / "analysis"
            / "L5_TPC_A_cumulative_section_path_distances_volumes_individual"
            / "0_L5_TPC_A_cumulative_section_path_distances_volumes_C030796A-P3_lite.h5.pdf"
        ),
    )

    # Test with inconsistent neurondb
    single_pop_neurondb.loc[0, "morphology"] = "UNKNOWN_MORPHOLOGY"
    single_pop_neurondb.to_csv(
        single_pop_data_dir / "neurondb.dat", sep=" ", header=False, index=False
    )
    with pytest.raises(DiameterSynthesisError):
        plotting.cumulative_analysis(
            single_pop_data_dir,
            single_pop_diametrized_data_dir,
            tmpdir / "analysis",
            True,
            single_pop_data_dir / "neurondb.dat",
            ["basal", "apical"],
            ext=".pdf",
        )


def test_violin_analysis(
    tmpdir, expected_images_path, single_pop_data_dir, single_pop_diametrized_data_dir
):
    """Test the violin_analysis function."""
    # Plot the figures
    plotting.violin_analysis(
        single_pop_data_dir,
        single_pop_diametrized_data_dir,
        tmpdir / "analysis",
        single_pop_data_dir / "neurondb.dat",
    )

    # Check the figures
    assert pdfdiff(
        str(expected_images_path / "test_violin_analysis.pdf"),
        str(tmpdir / "analysis" / "morphometrics.pdf"),
    )
