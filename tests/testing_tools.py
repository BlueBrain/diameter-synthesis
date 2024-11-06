"""Some tools used in tests."""

# Copyright (C) 2021-2024  Blue Brain Project, EPFL
#
# SPDX-License-Identifier: Apache-2.0

from numpy.testing import assert_allclose


def compare_diameters(neuron_a, neuron_b, **kwargs):
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
