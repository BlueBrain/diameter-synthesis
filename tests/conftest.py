"""Prepare the tests."""

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

# pylint: disable=redefined-outer-name
import json
import shutil
from pathlib import Path

import morphio.mut
import neurom as nm
import numpy as np
import pandas as pd
import pytest
from morphio import PointLevel
from morphio import SectionType


@pytest.fixture(autouse=True)
def set_random():
    """Path to the test directory."""
    np.random.seed(0)


@pytest.fixture
def test_path():
    """Path to the test directory."""
    return Path(__file__).absolute().parent


@pytest.fixture
def test_data_path(test_path):
    """Path to the test data directory."""
    return test_path / "data"


@pytest.fixture
def expected_images_path(test_path):
    """Path to the expected images directory."""
    return test_path / "expected_images"


@pytest.fixture
def config_path(test_data_path):
    """Path to the config file."""
    return test_data_path / "config.json"


@pytest.fixture
def config(config_path):
    """The config used in tests."""
    with open(config_path, encoding="utf-8") as filename:
        return json.load(filename)


@pytest.fixture
def model_params_path(test_data_path):
    """The path to the model parameters."""
    return test_data_path / "model_params.json"


@pytest.fixture
def model_params(model_params_path):
    """The model parameters."""
    with open(model_params_path, encoding="utf-8") as filename:
        return json.load(filename)


@pytest.fixture
def model_data(test_data_path):
    """The model data."""
    with open(test_data_path / "model_data.json", "r", encoding="utf-8") as filename:
        return json.load(filename)


@pytest.fixture
def neuron_path(test_data_path):
    """The path to the tested morphology."""
    return test_data_path / "C030796A-P3_lite.h5"


@pytest.fixture
def neuron(neuron_path):
    """The tested morphology."""
    return morphio.mut.Morphology(neuron_path)


@pytest.fixture
def neuron_diametrized_path(test_data_path):
    """The path to the diametrized morphology."""
    return test_data_path / "C030796A-P3_lite_diametrized.h5"


@pytest.fixture
def neuron_diametrized(neuron_diametrized_path):
    """The diametrized morphology."""
    return morphio.mut.Morphology(neuron_diametrized_path)


@pytest.fixture
def single_pop(neuron_path):
    """A population with one morphology."""
    return [nm.load_morphology(neuron_path)]


@pytest.fixture
def single_pop_diametrized(neuron_diametrized_path):
    """A population with one diametrized morphology."""
    return [nm.load_morphology(neuron_diametrized_path)]


@pytest.fixture
def single_neurite(neuron_path):
    """A neurite."""
    neuron = nm.load_morphology(neuron_path)
    yield neuron.neurites[0]


@pytest.fixture
def single_pop_neurondb_dat_path(test_data_path):
    """The path to the neurondb.dat file."""
    return test_data_path / "neurondb.dat"


@pytest.fixture
def single_pop_neurondb(single_pop_neurondb_dat_path):
    """The DF of the singme population."""
    df = pd.read_csv(single_pop_neurondb_dat_path, header=None, sep=r"\s+")
    df.rename(columns={0: "morphology", 1: "layer", 2: "mtype"}, inplace=True)
    return df


@pytest.fixture
def single_pop_data_dir(tmpdir, single_pop_neurondb, neuron_path):
    """Prepare test data."""
    single_pop_dir = tmpdir / "single_pop"

    single_pop_dir.mkdir()
    single_pop_neurondb.to_csv(single_pop_dir / "neurondb.dat", sep=" ", header=False, index=False)
    shutil.copyfile(neuron_path, single_pop_dir / neuron_path.name)

    return single_pop_dir


@pytest.fixture
def single_pop_diametrized_data_dir(
    tmpdir, single_pop_neurondb, neuron_path, neuron_diametrized_path
):
    """Prepare test data with diametrized morphologies."""
    single_pop_diametrized_dir = tmpdir / "single_pop_diametrized"
    single_pop_diametrized_dir.mkdir()
    single_pop_neurondb.loc[0, "morphology"] += "_diametrized"
    single_pop_neurondb.to_csv(
        single_pop_diametrized_dir / "neurondb.dat", sep=" ", header=False, index=False
    )
    shutil.copyfile(neuron_diametrized_path, single_pop_diametrized_dir / neuron_path.name)

    return single_pop_diametrized_dir


@pytest.fixture
def empty_build_result():
    """The result of an empty build."""
    return (
        {
            "diameter_power_relation": {"apical": {}, "basal": {}},
            "sibling_ratios": {"apical": {}, "basal": {}},
            "tapers": {"apical": {}, "basal": {}},
            "terminal_diameters": {"apical": {}, "basal": {}},
            "trunk_diameters": {"apical": {}, "basal": {}},
        },
        {
            "diameter_power_relation": {"apical": [], "basal": []},
            "sibling_ratios": {"apical": [], "basal": []},
            "tapers": {"apical": [], "basal": []},
            "terminal_diameters": {"apical": [], "basal": []},
            "trunk_diameters": {"apical": [], "basal": []},
        },
    )


@pytest.fixture
def small_morph():
    """A small morphology used in tests."""
    N = 3
    dx = 100
    dy = 100
    Ndy = (N - 1) * dy

    nrn = morphio.mut.Morphology()

    # Define a basal with only one section
    x0_b1 = 0
    nrn.append_root_section(
        PointLevel([[x0_b1, i * dy, 0] for i in range(N)], [0.2] * N), SectionType(3)
    )

    # Define a basal with several sections
    x0_b2 = 0.1
    root = nrn.append_root_section(
        PointLevel([[x0_b2, i * dy, 0] for i in range(N)], [0.2] * N), SectionType(3)
    )
    root.append_section(
        PointLevel([[x0_b2, Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(3)
    )
    root.append_section(
        PointLevel([[x0_b2 + i * -0.1, Ndy + i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(3),
    )

    # Define the apical with several bifurcations
    x0_ap = 0.2
    apical_root = nrn.append_root_section(
        PointLevel([[x0_ap, i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )
    apical_root.append_section(
        PointLevel([[x0_ap + i * dx, Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )
    major = apical_root.append_section(
        PointLevel([[x0_ap, Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )

    major.append_section(
        PointLevel([[x0_ap + i * dx, 2 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(4),
    )
    major = major.append_section(
        PointLevel([[x0_ap, 2 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )

    major.append_section(
        PointLevel([[x0_ap + i * dx, 3 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(4),
    )
    major = major.append_section(
        PointLevel([[x0_ap, 3 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )

    major.append_section(
        PointLevel([[x0_ap + i * dx, 4 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(4),
    )
    major = major.append_section(
        PointLevel([[x0_ap, 4 * Ndy + i * dy, 0] for i in range(N)], [0.2] * N), SectionType(4)
    )

    # Define an axon with one section
    x0_ax1 = 0
    y0_ax1 = -0.1
    nrn.append_root_section(
        PointLevel([[x0_ax1, y0_ax1 - i * dy, 0] for i in range(N)], [0.2] * N), SectionType(2)
    )

    # Define an axon with one bifurcation
    x0_ax2 = 0.1
    y0_ax2 = -0.1
    root = nrn.append_root_section(
        PointLevel([[x0_ax2, y0_ax2 - i * dy, 0] for i in range(N)], [0.2] * N), SectionType(2)
    )
    root.append_section(
        PointLevel([[x0_ax2, y0_ax2 - Ndy - i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(3),
    )
    root.append_section(
        PointLevel([[x0_ax2 + i * -0.1, y0_ax2 - Ndy - i * dy, 0] for i in range(N)], [0.2] * N),
        SectionType(3),
    )

    # Define soma according to the root sections
    nrn.soma.points = [
        [x0_b1, 0, 0],
        [x0_b2, 0, 0],
        [x0_ap, 0, 0],
        [x0_ax2, y0_ax2, 0],
        [x0_ax1, y0_ax1, 0],
    ]
    nrn.soma.diameters = [0] * len(nrn.soma.points)

    return nrn
