"""Some tools used in tests."""

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

import neurom as nm
import numpy as np
from numpy.testing import assert_allclose


def _copy_diameters(neuron_a, neuron_b):
    for section_a, section_b in zip(nm.iter_sections(neuron_a), nm.iter_sections(neuron_b)):
        section_a.diameters = section_b.diameters


def _compare_diameters(neuron_a, neuron_b, **kwargs):
    """Test if all diameters are the same."""
    errs = []
    for section_a, section_b in zip(neuron_a.sections, neuron_b.sections):
        try:
            assert_allclose(
                neuron_a.sections[section_a].diameters,
                neuron_b.sections[section_b].diameters,
                err_msg=f"Error for section {section_a}",
                **kwargs,
            )
        except AssertionError as e:
            errs.append(e)

    if errs:
        raise AssertionError("\n".join([i.args[0] for i in errs]))


def nested_round(obj, precision=6):
    """Round all floats (recursively) in a nested dictionary."""
    # pylint: disable=too-many-return-statements
    if precision is None:
        precision = 0
    if isinstance(obj, np.ndarray):
        return [nested_round(i, precision) for i in obj.tolist()]
    if isinstance(obj, np.floating):
        return round(float(obj), precision)
    if isinstance(obj, np.integer):
        return round(int(obj), precision)
    if isinstance(obj, float):
        return round(obj, precision)
    if isinstance(obj, dict):
        return dict((k, nested_round(v, precision)) for k, v in obj.items())
    if isinstance(obj, (list, tuple)):
        return [nested_round(i, precision) for i in obj]
    return obj


def compare_dicts(ref, test, precision=None):
    """Compare two dictionaries."""
    if precision is not None:
        ref = nested_round(ref, precision)
        test = nested_round(test, precision)
    return json.dumps(ref, sort_keys=True) == json.dumps(test, sort_keys=True)
