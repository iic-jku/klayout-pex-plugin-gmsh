#
# --------------------------------------------------------------------------------
# SPDX-FileCopyrightText: 2026 Martin Jan Köhler and Harald Pretl
# Johannes Kepler University, Institute for Integrated Circuits.
#
# This file is part of klayout-pex-plugin-gmsh
# (see https://github.com/iic-jku/klayout-pex-plugin-gmsh).
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
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# SPDX-License-Identifier: GPL-3.0-or-later
# --------------------------------------------------------------------------------
#

"""
Scenes for the tests.

``two_nets.pex25d`` is hand-written: two nets (one with a via), a floating
conductor, and a film wrapping a film that covers the field.
``overlap_plates_100um_x_100um_m1_m2.pex25d`` is what ``kpex pex25d`` writes for
the KLayout-PEX test design of that name in ihp-sg13g2.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from klayout_pex import pex25d

DATA = Path(__file__).parent.parent / 'testdata'


def load_scene(name: str):
    return pex25d.resolve(pex25d.read(str(DATA / name)))


@pytest.fixture(scope='session')
def two_nets():
    return load_scene('two_nets.pex25d')


@pytest.fixture(scope='session')
def overlap_plates():
    return load_scene('overlap_plates_100um_x_100um_m1_m2.pex25d')
