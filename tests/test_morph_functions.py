"""Test the morph_functions module."""

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

import pytest

from diameter_synthesis import morph_functions
from diameter_synthesis.exception import DiameterSynthesisError


def test_partition_asymmetry_length(neuron):
    """Test the partition_asymmetry_length function."""
    section = neuron.root_sections[0]
    res = morph_functions.partition_asymmetry_length(section)
    assert res == pytest.approx(157.6603)

    # Test with only one child section
    one_child_section = neuron.root_sections[1]
    neuron.delete_section(one_child_section.children[0], recursive=True)
    with pytest.raises(DiameterSynthesisError):
        morph_functions.partition_asymmetry_length(one_child_section)


def test_compute_sibling_ratios(single_neurite):
    """Test the compute_sibling_ratios function."""
    res = morph_functions.compute_sibling_ratios(single_neurite, "mean")
    assert res == pytest.approx([0.508197, 0.434402, 1, 1, 0.395023, 1, 1])

    res_bounds = morph_functions.compute_sibling_ratios(
        single_neurite, "mean", bounds=[0, 1 + 1e-5]
    )
    assert res_bounds == pytest.approx(res)


def test_compute_diameter_power_relation(single_neurite):
    """Test the compute_diameter_power_relation function."""
    res = morph_functions.compute_diameter_power_relation(single_neurite, "mean")
    assert res == pytest.approx(
        [6.964765, 4.492706, 5.942041, 2.534333, 7.254541, 2.822203, 16.933426]
    )

    res_bounds = morph_functions.compute_diameter_power_relation(
        single_neurite, "mean", bounds=[0, 100]
    )
    assert res_bounds == pytest.approx(res)


def test_terminal_diameters(single_neurite):
    """Test the terminal_diameters function."""
    res = morph_functions.terminal_diameters(single_neurite, "mean")
    assert res == pytest.approx([0.3, 0.3, 0.23, 0.23, 0.3, 0.3, 0.23, 0.23])

    res_bounds = morph_functions.terminal_diameters(single_neurite, "mean", bounds=[0, 100])
    assert res_bounds == pytest.approx(res)

    res_first = morph_functions.terminal_diameters(single_neurite, "first")
    assert res_first == pytest.approx(res)

    with pytest.raises(DiameterSynthesisError):
        morph_functions.terminal_diameters(single_neurite, "UNKNOWN")


def test_min_diameter(single_neurite):
    """Test the min_diameter function."""
    res = morph_functions.min_diameter(single_neurite)
    assert res == pytest.approx([0.23])

    res_bounds = morph_functions.min_diameter(single_neurite, bounds=[0, 100])
    assert res_bounds == pytest.approx(res)


def test_max_diameter(single_neurite):
    """Test the max_diameter function."""
    res = morph_functions.max_diameter(single_neurite)
    assert res == pytest.approx([1.84])

    res_bounds = morph_functions.max_diameter(single_neurite, bounds=[0, 100])
    assert res_bounds == pytest.approx(res)


def test_trunk_diameter(single_neurite):
    """Test the trunk_diameter function."""
    res = morph_functions.trunk_diameter(single_neurite)
    assert res == pytest.approx([1.84])

    res_bounds = morph_functions.trunk_diameter(single_neurite, bounds=[0, 100])
    assert res_bounds == pytest.approx(res)

    res_bounds = morph_functions.trunk_diameter(single_neurite, method="mean")
    assert res_bounds == pytest.approx(res)

    res_bounds = morph_functions.trunk_diameter(single_neurite, method="first")
    assert res_bounds == pytest.approx(res)


def test_taper(single_neurite):
    """Test the taper function."""
    params = {"min": -1e-9, "max": 1e-9}
    res = morph_functions.taper(single_neurite, params)
    assert res == pytest.approx([0.0] * 12)

    res_tot_length = morph_functions.taper(single_neurite, params, "tot_length")
    assert len(res_tot_length) == 1
    assert res_tot_length[0] == pytest.approx([0.0, 828.830532])

    # Test with None params
    with pytest.raises(DiameterSynthesisError):
        morph_functions.taper(single_neurite, None)


def test_get_additional_attribute(single_neurite):
    """Test the get_additional_attribute function."""
    res_asymmetry = morph_functions.get_additional_attribute("asymmetry", single_neurite)
    assert res_asymmetry == pytest.approx(
        [0.19022014, 0.02628145, 0.16578421, 0.16073337, 0.00042449, 0.00093826, 0.02613346],
        rel=1e-5,
    )

    res_asymmetry_threshold = morph_functions.get_additional_attribute(
        "asymmetry_threshold", single_neurite
    )
    assert res_asymmetry_threshold == pytest.approx(res_asymmetry)

    section = next(single_neurite.iter_sections())
    res_asymmetry_pair = morph_functions.get_additional_attribute("asymmetry_pair", section=section)
    assert len(res_asymmetry_pair) == 1
    assert res_asymmetry_pair[0] == pytest.approx([7, 7])

    res_tot_length = morph_functions.get_additional_attribute("tot_length", single_neurite)
    assert res_tot_length == pytest.approx([828.83])

    res_max_path = morph_functions.get_additional_attribute("max_path", single_neurite)
    assert res_max_path == pytest.approx([325.016])

    res_max_branch = morph_functions.get_additional_attribute("max_branch", single_neurite)
    assert res_max_branch == pytest.approx([3])

    res_root_strahler = morph_functions.get_additional_attribute("root_strahler", single_neurite)
    assert res_root_strahler == pytest.approx([4])

    res_sibling = morph_functions.get_additional_attribute("sibling", single_neurite)
    assert res_sibling == pytest.approx([0.508197, 0.434402, 1, 1, 0.395023, 1, 1])

    res_sibling_section = morph_functions.get_additional_attribute("sibling", section=section)
    assert res_sibling_section == pytest.approx(0.508197)

    # Test with unknown attribute name
    with pytest.raises(DiameterSynthesisError):
        morph_functions.get_additional_attribute("UNKNOWN", single_neurite)

    # Test with no neurite nor section
    with pytest.raises(DiameterSynthesisError):
        morph_functions.get_additional_attribute("asymmetry")

    # Test with both neurite and section
    for attr in [
        "asymmetry",
        "asymmetry_threshold",
        "asymmetry_pair",
        "tot_length",
        "max_path",
        "max_branch",
        "root_strahler",
        "sibling",
    ]:
        with pytest.raises(DiameterSynthesisError):
            morph_functions.get_additional_attribute(attr)


def test_add_additional_attributes(single_neurite):
    """Test the add_additional_attributes function."""
    res_max_branch = morph_functions.add_additional_attributes(
        [-1, 1, 10], single_neurite, "max_branch"
    )
    assert len(res_max_branch) == 1
    assert res_max_branch[0] == pytest.approx([-1, 3])

    res_max_branch_bounds = morph_functions.add_additional_attributes(
        [-1, 1, 10], single_neurite, "max_branch", bounds=[-1000, 1000]
    )
    assert len(res_max_branch_bounds) == 1
    assert res_max_branch_bounds[0] == pytest.approx([-1, 3])
