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

"""The occupancy rule and the WRAPS growth, checked point by point."""

from __future__ import annotations

from itertools import combinations

import klayout.db as kdb
import pytest

from klayout_pex_plugins.gmsh.solids import SceneSolids, scene_solids


def body_at(solids: SceneSolids, x_um: float, y_um: float, z_um: float) -> str:
    x, y, z = (round(value / solids.grid_um) for value in (x_um, y_um, z_um))
    probe = kdb.Region(kdb.Box(x, y, x + 1, y + 1))
    names = [prism.body.name for prism in solids.prisms
             if prism.zlow <= z < prism.zhigh and not (prism.region & probe).is_empty()]
    assert len(names) == 1, names
    return names[0]


@pytest.mark.parametrize('scene', ['two_nets', 'overlap_plates'])
def test_prisms_fill_the_domain_exactly_once(scene: str, request: pytest.FixtureRequest):
    solids = scene_solids(request.getfixturevalue(scene), field_margin_um=8.0)
    volume = sum(prism.region.area() * (prism.zhigh - prism.zlow) for prism in solids.prisms)
    assert volume == solids.field.area() * (solids.zhigh - solids.zlow)
    for a, b in combinations(solids.prisms, 2):
        if a.zlow < b.zhigh and b.zlow < a.zhigh:
            assert (a.region & b.region).is_empty(), (a.body.name, b.body.name)


@pytest.mark.parametrize('point, body', [
    ((0.5, 0.25, 1.2), 'neta'),        # met1 of A
    ((0.3, 0.2, 1.7), 'neta'),         # via1 of A, same net as both metals
    ((2.5, 1.75, 2.2), 'F'),           # a floating conductor is its own body
    ((-1.5, -1.5, -0.2), 'subs'),      # the ground plane spans the field
    ((0.5, 0.25, 1.45), 'lint'),       # film over the metal ...
    ((0.5, 0.25, 1.55), 'cap'),        # ... and the film on that film
    ((1.05, 0.25, 1.2), 'lint'),       # beside the metal ...
    ((1.15, 0.25, 1.2), 'cap'),        # ... grown once more by the outer film
    ((1.25, 0.25, 1.2), 'ild'),        # the band fills what the films leave
    ((-1.5, -1.5, 1.02), 'lint'),      # on the field: lint from met1's bottom ...
    ((-1.5, -1.5, 1.1), 'cap'),        # ... and cap on top of lint, not inside it
    ((-1.5, -1.5, 0.5), 'fox'),
    ((-1.5, -1.5, 2.2), 'air'),        # the background above the last band
], ids=lambda value: value if isinstance(value, str) else None)
def test_occupancy(two_nets, point, body: str):
    assert body_at(scene_solids(two_nets, field_margin_um=8.0), *point) == body


def test_domain_comes_from_the_scene_or_the_margin(two_nets, overlap_plates):
    solids = scene_solids(two_nets, field_margin_um=8.0)       # DOMAIN_MARGIN X 2 Y 2 Z 1
    assert (solids.field.left, solids.field.right) == (-2000, 5000)
    assert (solids.zlow, solids.zhigh) == (-400, 3400)

    solids = scene_solids(overlap_plates, field_margin_um=4.0)   # no domain in the file
    assert (solids.field.left, solids.field.right) == (-40000, 1540000)
    assert solids.zhigh - 40000 == max(prism.zhigh for prism in solids.prisms
                                       if prism.body.name != 'air')
