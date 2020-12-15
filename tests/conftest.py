import json
from pathlib import Path

import morphio.mut
import neurom as nm
import pandas as pd
import pytest
import shutil


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
    return nm.load_neurons([neuron_path])


@pytest.fixture
def single_pop_diametrized(neuron_diametrized_path):
    return nm.load_neurons([neuron_diametrized_path])


@pytest.fixture
def single_neurite(single_pop):
    return single_pop.neurites[0]


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
