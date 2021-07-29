import json
from pathlib import Path

import morphio.mut
import neurom as nm
import pandas as pd
import pytest
import shutil
from morphio import PointLevel
from morphio import SectionType


@pytest.fixture
def test_path():
    return Path(__file__).absolute().parent


@pytest.fixture
def test_data_path(test_path):
    return test_path / "data"


@pytest.fixture
def expected_images_path(test_path):
    return test_path / "expected_images"


@pytest.fixture
def config_path(test_data_path):
    return test_data_path / "config.json"


@pytest.fixture
def config(config_path):
    with open(config_path) as filename:
        return json.load(filename)


@pytest.fixture
def model_params_path(test_data_path):
    return test_data_path / "model_params.json"


@pytest.fixture
def model_params(model_params_path):
    with open(model_params_path) as filename:
        return json.load(filename)


@pytest.fixture
def model_data(test_data_path):
    with open(test_data_path / "model_data.json", "r") as filename:
        return json.load(filename)


@pytest.fixture
def neuron_path(test_data_path):
    return test_data_path / "C030796A-P3_lite.h5"


@pytest.fixture
def neuron(neuron_path):
    return morphio.mut.Morphology(neuron_path)


@pytest.fixture
def neuron_diametrized_path(test_data_path):
    return test_data_path / "C030796A-P3_lite_diametrized.h5"


@pytest.fixture
def neuron_diametrized(neuron_diametrized_path):
    return morphio.mut.Morphology(neuron_diametrized_path)


@pytest.fixture
def single_pop(neuron_path):
    return [nm.load_morphology(neuron_path)]


@pytest.fixture
def single_pop_diametrized(neuron_diametrized_path):
    return [nm.load_morphology(neuron_diametrized_path)]


@pytest.fixture
def single_neurite(neuron_path):
    neuron = nm.load_morphology(neuron_path)
    yield neuron.neurites[0]


@pytest.fixture
def single_pop_neurondb_dat_path(test_data_path):
    return test_data_path / "neurondb.dat"


@pytest.fixture
def single_pop_neurondb(single_pop_neurondb_dat_path):
    df = pd.read_csv(single_pop_neurondb_dat_path, header=None, sep=r"\s+")
    df.rename(columns={0: "morphology", 1: "layer", 2: "mtype"}, inplace=True)
    return df


@pytest.fixture
def single_pop_data_dir(tmpdir, single_pop_neurondb, neuron_path):
    # Prepare test data
    single_pop_dir = tmpdir / "single_pop"

    single_pop_dir.mkdir()
    single_pop_neurondb.to_csv(single_pop_dir / "neurondb.dat", sep=" ", header=False, index=False)
    shutil.copyfile(neuron_path, single_pop_dir / neuron_path.name)

    return single_pop_dir


@pytest.fixture
def single_pop_diametrized_data_dir(
    tmpdir, single_pop_neurondb, neuron_path, neuron_diametrized_path
):
    # Prepare test data
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
    N = 3
    dx = 100
    dy = 100
    Ndx = (N - 1) * dx
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
