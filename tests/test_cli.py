"""Tests for the diameter_synthesis.cli module."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

import json
from pathlib import Path

import dictdiffer
import pytest
from click.testing import CliRunner
from numpy.testing import assert_almost_equal

from diameter_synthesis import cli
from diameter_synthesis import main


def test_cli(cli_runner):
    # pylint: disable=unused-argument
    """Test the CLI."""
    result = cli_runner.invoke(
        cli.main,
        [
            "--version",
        ],
    )
    assert result.exit_code == 0
    assert result.output.startswith("diameter-synthesis, version ")


@pytest.mark.parametrize("command", ["run_models", "run_diameters", "plot_diff", "run_analysis"])
def test_entry_points(script_runner, command):
    """Test the entry points."""
    ret = script_runner.run("diameter-synthesis", command)
    assert not ret.success
    assert f"Usage: diameter-synthesis {command}" in ret.stderr
    assert ret.stdout == ""


def test_run_models(tmpdir, single_pop_data_dir, single_pop_diametrized_data_dir, config):
    """Test the run_models entry point."""
    # Prepare inputs
    extract_models_params = config
    extract_models_params["mtypes_file"] = str(single_pop_data_dir / "neurondb.dat")
    extract_models_params["morph_path"] = str(single_pop_data_dir)
    extract_models_params["new_morph_path"] = str(single_pop_diametrized_data_dir)
    extract_models_params["models_params_file"] = str(tmpdir / "model_params_mtypes.json")
    extract_models_params["fig_folder"] = str(tmpdir / "model_figures")
    extract_models_params["n_cpu"] = 1
    extract_models_params["neurite_types"] = ["basal_dendrite", "apical_dendrite"]

    config_file = str(tmpdir / "diametrizer_params.json")
    with open(config_file, "w", encoding="utf-8") as json_file:
        json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

    # Run with CLI
    runner = CliRunner()
    runner.invoke(cli.main, ["run_models", config_file], catch_exceptions=False)

    # Check results
    with open(extract_models_params["models_params_file"], "r", encoding="utf-8") as json_file:
        res = json.load(json_file)

    assert list(res.keys()) == ["generic"]
    assert list(res["generic"].keys()) == ["L5_TPC:A"]

    # Check only diameter_power_relation entry

    assert dictdiffer.diff(
        res["generic"]["L5_TPC:A"]["diameter_power_relation"],
        {
            "apical_dendrite": {
                "distribution": "exponnorm",
                "params": {
                    "a": 0.3,
                    "loc": 10.328229,
                    "max": 21.302097,
                    "min": 3.745589,
                    "num_value": 4,
                    "scale": 7.391765,
                },
                "sequential": None,
            },
            "basal_dendrite": {
                "distribution": "exponnorm",
                "params": {
                    "a": 4.0,
                    "loc": 2.0812598,
                    "max": 13.171166,
                    "min": 2.431371,
                    "num_value": 10,
                    "scale": 0.815449,
                },
                "sequential": None,
            },
        },
        absolute_tolerance=1e-5,
    )


def test_run_diameters(tmpdir, single_pop_data_dir, config, model_params):
    """Test the run_diameters entry point."""
    # Prepare inputs
    res_path = Path(tmpdir / "new_morphologies")
    extract_models_params = config
    extract_models_params["mtypes_file"] = str(single_pop_data_dir / "neurondb.dat")
    extract_models_params["morph_path"] = str(single_pop_data_dir)
    extract_models_params["n_cpu"] = 1
    extract_models_params["neurite_types"] = ["basal_dendrite", "apical_dendrite"]
    extract_models_params["new_morph_path"] = str(res_path)

    config_file = str(tmpdir / "diametrizer_params.json")
    with open(config_file, "w", encoding="utf-8") as json_file:
        json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

    model_params_file = str(tmpdir / "model_params.json")
    with open(model_params_file, "w", encoding="utf-8") as json_file:
        json.dump({"generic": {"L5_TPC:A": model_params}}, json_file, sort_keys=True, indent=4)

    # Run with CLI
    runner = CliRunner()
    runner.invoke(
        cli.main, ["run_diameters", config_file, model_params_file], catch_exceptions=False
    )

    # Check results
    assert [i.name for i in res_path.iterdir()] == ["C030796A-P3_lite.h5"]


def test_run_diametrize_single_neuron(neuron):
    """Test diametrize single neuron."""
    n_root = len(neuron.root_sections)
    main.diametrize_single_neuron(neuron)
    assert len(neuron.root_sections) == n_root
    assert_almost_equal(neuron.root_sections[1].diameters, [1.627442, 1.6274352])


def test_run_diametrize_single_neuron_basal(neuron):
    """Test diametrize single neuron with only basal."""
    main.diametrize_single_neuron(neuron, neurite_types=["basal_dendrite"])
    assert_almost_equal(neuron.root_sections[1].diameters, [1.627442, 1.6274352])


def test_plot_diff(tmpdir, single_pop_data_dir, single_pop_diametrized_data_dir):
    """Test the plot_diff entry point."""
    res_path = Path(tmpdir / "figures")

    # Run with CLI
    runner = CliRunner()
    runner.invoke(
        cli.main,
        [
            "plot_diff",
            "--orig-path",
            str(single_pop_data_dir),
            "--diam-path",
            str(single_pop_diametrized_data_dir),
            "--out-dir",
            str(res_path),
            "--ext",
            ".pdf",
        ],
        catch_exceptions=False,
    )

    # Check results
    assert [i.name for i in (res_path / "diffs").iterdir()] == ["C030796A-P3_lite.pdf"]


def test_run_analysis(tmpdir, single_pop_data_dir, single_pop_diametrized_data_dir):
    """Test the run_analysis entry point."""
    res_path = Path(tmpdir / "figures")

    # Run with CLI
    runner = CliRunner()
    runner.invoke(
        cli.main,
        [
            "run_analysis",
            "--orig-path",
            str(single_pop_data_dir),
            "--diam-path",
            str(single_pop_diametrized_data_dir),
            "--mtypes-file",
            str(single_pop_data_dir / "neurondb.dat"),
            "--out-dir",
            str(res_path),
            "--cumulative",
            "--individual",
            "--violin",
            "--ext",
            ".pdf",
        ],
        catch_exceptions=False,
    )

    # Check results
    assert sorted(
        [i.relative_to(res_path).as_posix() for i in res_path.glob("**") if i != res_path]
    ) == [
        "apical_dendrite",
        "apical_dendrite/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "apical_dendrite/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
        "axon",
        "axon/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "axon/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
        "basal_dendrite",
        "basal_dendrite/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "basal_dendrite/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
    ]

    # Test with neither cumukative nor violin
    with pytest.raises(ValueError):
        runner.invoke(
            cli.main,
            [
                "run_analysis",
                "--orig-path",
                str(single_pop_data_dir),
                "--diam-path",
                str(single_pop_diametrized_data_dir),
                "--mtypes-file",
                str(single_pop_data_dir / "neurondb.dat"),
                "--out-dir",
                str(res_path),
                "--ext",
                ".pdf",
            ],
            catch_exceptions=False,
        )
