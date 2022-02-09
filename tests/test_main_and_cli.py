"""Test the cli module."""

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

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from diameter_synthesis import cli

from .testing_tools import compare_dicts


def test_run_models(tmpdir, single_pop_data_dir, single_pop_diametrized_data_dir, config):
    """Test the run_models entry point."""
    # Prepare inputs
    extract_models_params = config["generic"]
    extract_models_params["mtypes_file"] = str(single_pop_data_dir / "neurondb.dat")
    extract_models_params["morph_path"] = str(single_pop_data_dir)
    extract_models_params["new_morph_path"] = str(single_pop_diametrized_data_dir)
    extract_models_params["models_params_file"] = str(tmpdir / "model_params_mtypes.json")
    extract_models_params["fig_folder"] = str(tmpdir / "model_figures")
    extract_models_params["n_cpu"] = 1
    extract_models_params["neurite_types"] = ["basal", "apical"]

    config_file = str(tmpdir / "diametrizer_params.json")
    with open(config_file, "w", encoding="utf-8") as json_file:
        json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

    # Run with CLI
    runner = CliRunner()
    runner.invoke(cli.cli, ["run_models", config_file], catch_exceptions=False)

    # Check results
    with open(extract_models_params["models_params_file"], "r", encoding="utf-8") as json_file:
        res = json.load(json_file)

    assert list(res.keys()) == ["generic"]
    assert list(res["generic"].keys()) == ["L5_TPC:A"]

    # Check only diameter_power_relation entry
    compare_dicts(
        res["generic"]["L5_TPC:A"]["diameter_power_relation"],
        {
            "apical": {
                "distribution": "exponnorm",
                "params": {
                    "a": 0.3,
                    "loc": 10.324635885438347,
                    "max": 21.302095794677733,
                    "min": 3.745589017868042,
                    "num_value": 4,
                    "scale": 7.39024091228236,
                },
                "sequential": "asymmetry_threshold",
            },
            "basal": {
                "distribution": "exponnorm",
                "params": {
                    "a": 4.0,
                    "loc": 2.081371164504187,
                    "max": 13.171166133880611,
                    "min": 2.4313707590103153,
                    "num_value": 10,
                    "scale": 0.8153109268097628,
                },
                "sequential": "asymmetry_threshold",
            },
        },
        precision=6,
    )


def test_run_diameters(tmpdir, single_pop_data_dir, config, model_params_path):
    """Test the run_diameters entry point."""
    # Prepare inputs
    res_path = Path(tmpdir / "new_morphologies")
    extract_models_params = config["generic"]
    extract_models_params["mtypes_file"] = str(single_pop_data_dir / "neurondb.dat")
    extract_models_params["morph_path"] = str(single_pop_data_dir)
    extract_models_params["n_cpu"] = 1
    extract_models_params["neurite_types"] = ["basal", "apical"]
    extract_models_params["new_morph_path"] = str(res_path)

    config_file = str(tmpdir / "diametrizer_params.json")
    with open(config_file, "w", encoding="utf-8") as json_file:
        json.dump(extract_models_params, json_file, sort_keys=True, indent=4)

    # Run with CLI
    runner = CliRunner()
    runner.invoke(
        cli.cli, ["run_diameters", config_file, str(model_params_path)], catch_exceptions=False
    )

    # Check results
    assert [i.name for i in res_path.iterdir()] == ["C030796A-P3_lite.h5"]


def test_plot_diff(tmpdir, single_pop_data_dir, single_pop_diametrized_data_dir):
    """Test the plot_diff entry point."""
    res_path = Path(tmpdir / "figures")

    # Run with CLI
    runner = CliRunner()
    runner.invoke(
        cli.cli,
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
        cli.cli,
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
        "apical",
        "apical/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "apical/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
        "axon",
        "axon/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "axon/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
        "basal",
        "basal/L5_TPC_A_cumulative_section_path_distances_areas_individual",
        "basal/L5_TPC_A_cumulative_section_path_distances_volumes_individual",
    ]

    # Test with neither cumukative nor violin
    with pytest.raises(ValueError):
        runner.invoke(
            cli.cli,
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
