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

"""Ground openings and typed ownership survive prism construction and real meshing."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from unittest.mock import patch

import klayout.db as kdb
import pytest

from klayout_pex import pex25d
from klayout_pex.plugin_api.v1 import SceneCapability
from klayout_pex_plugins.gmsh import GmshMeshOptions, GmshSceneExporter, write_mesh
from klayout_pex_plugins.gmsh.solids import BodyKind, _SolidsBuilder, scene_solids
from .test_exporter import meshed_volumes, physical_groups
from .test_solids import body_at


def test_cutout_union_and_holes_preserve_reference_islands_and_coatings(ground_cutout: Any):
    before = ground_cutout.SerializeToString(deterministic=True)
    solids = scene_solids(ground_cutout, 8)
    ground = [p for p in solids.prisms if p.body.name == 'VSS']
    # Independent rectangle arithmetic: 36 - (16 - 4 + 2) um².
    assert sum(p.region.area()*(p.zhigh-p.zlow) for p in ground) == 22_000_000*400
    for x, y in [(-2.5, -2.5), (0, 0)]:
        assert body_at(solids, x, y, -.2) == 'VSS'
        assert body_at(solids, x, y, .05) == 'oxide'
    assert body_at(solids, -1.5, 0.25, -.2) == 'air'
    assert body_at(solids, -1.5, 0.25, .05) == 'signal'
    assert body_at(solids, 2.5, 0, -.2) == 'air'
    assert ground_cutout.SerializeToString(deterministic=True) == before


def test_reference_net_merges_contacts_but_profile_name_remains_independent(ground_cutout: Any):
    ground_cutout.conductors[0].net = 'VSS'
    ground_cutout.ground_plane.ClearField('cutouts')
    solids = scene_solids(ground_cutout, 8)
    assert {b.name for b in solids.bodies if b.kind == BodyKind.CONDUCTOR} == {'VSS'}
    assert body_at(solids, -1.5, .25, .05) == 'VSS'
    assert body_at(solids, 0, 0, .05) == 'oxide'


def test_floating_name_collision_keeps_two_electrical_bodies(ground_cutout: Any):
    conductor = ground_cutout.conductors[0]
    conductor.name, conductor.net, conductor.floating = 'VSS', '', True
    solids = scene_solids(ground_cutout, 8)
    assert {b.name for b in solids.bodies if b.kind == BodyKind.CONDUCTOR} == {'VSS', 'FLOATING:VSS'}
    assert body_at(solids, -1.5, .25, .05) == 'FLOATING:VSS'


def test_no_ground_creates_no_reference_body(ground_cutout: Any):
    ground_cutout.ClearField('ground_plane')
    ground_cutout.ClearField('dielectrics')
    solids = scene_solids(ground_cutout, 8)
    assert {b.name for b in solids.bodies if b.kind == BodyKind.CONDUCTOR} == {'signal'}


def test_complete_reference_removal_fails_before_mesh_output(ground_cutout: Any, tmp_path: Path):
    ground_cutout.ground_plane.ClearField('cutouts')
    polygon = ground_cutout.ground_plane.cutouts.add().polygon
    for x, y in [(-4000, -4000), (4000, -4000), (4000, 4000), (-4000, 4000)]:
        polygon.outer.points.add(x=x, y=y)
    with pytest.raises(pex25d.ExportError, match='complete reference body'):
        write_mesh(ground_cutout, GmshMeshOptions(), str(tmp_path/'output'))
    assert not (tmp_path/'output').exists()


@pytest.mark.parametrize('axis', ['x', 'y'])
@pytest.mark.parametrize('edge,accepted', [(2147483646, True), (2147483647, False),
                                          (2147483648, False), (-2147483647, True),
                                          (-2147483648, False),
                                          (-2147483649, False)])
def test_native_field_boundaries_fail_before_box_construction(ground_cutout: Any, axis: str,
                                                             edge: int, accepted: bool):
    ground_cutout.ClearField('ground_plane')
    ground_cutout.ClearField('dielectrics')
    box = ground_cutout.domain.box
    setattr(box.lower_left, axis, edge-10 if edge > 0 else edge)
    setattr(box.upper_right, axis, edge if edge > 0 else edge+10)
    # Avoid invoking any native boolean with a rejected endpoint. The public
    # entry point must reject the field before native box construction wraps it.
    with patch.object(_SolidsBuilder, 'build', return_value=None), patch(
            'klayout_pex_plugins.gmsh.solids.kdb.Box', wraps=kdb.Box) as native_box:
        if accepted:
            scene_solids(ground_cutout, 8)
            native_box.assert_called_once()
        else:
            with pytest.raises(pex25d.ExportError, match='ground geometry cannot be represented'):
                scene_solids(ground_cutout, 8)
            native_box.assert_not_called()


@pytest.mark.parametrize('floating', [False, True])
def test_public_export_mesh_has_cutout_volume_and_distinct_terminal_groups(
        ground_cutout: Any, tmp_path: Path, floating: bool):
    if floating:
        conductor = ground_cutout.conductors[0]
        conductor.name, conductor.net, conductor.floating = 'VSS', '', True
    before = ground_cutout.SerializeToString(deterministic=True)
    paths = pex25d.export_with_backend(GmshSceneExporter(), ground_cutout, 'gmsh', str(tmp_path),
                                      prefix='cut_', options={'mesh_size_edge_um': .2,
                                                             'mesh_size_max_um': 2.0})
    groups = physical_groups(paths[0])
    assert {name for dim, name in groups if dim == 2} == {
        'VSS', 'FLOATING:VSS' if floating else 'signal', 'outer_boundary'}
    assert all(count > 0 for count in groups.values())
    # Domain is 36×1.4 um³; subtract ground, oxide and the .5×.5×.2 um contact.
    assert meshed_volumes(paths[0]) == pytest.approx({'oxide': 2.2, 'air': 39.35}, rel=1e-9)
    assert paths[:2] == [str(tmp_path/'cut_mesh.msh'), str(tmp_path/'cut_mesh.json')]
    if floating:
        assert paths[2:] == [str(tmp_path/'cut_body_names.json')]
        mapping = json.loads(Path(paths[2]).read_text(encoding='utf-8'))
        assert mapping['ground_net'] == 'VSS'
        assert {'kind': 'floating-body', 'name': 'VSS', 'solver_name': 'FLOATING:VSS'} in mapping['bodies']
    else:
        assert len(paths) == 2
    assert ground_cutout.SerializeToString(deterministic=True) == before


def test_capabilities_are_available_without_constructing_exporter():
    with patch.object(GmshSceneExporter, '__init__', side_effect=AssertionError('constructed')):
        assert GmshSceneExporter.supported_capabilities() == frozenset({
            SceneCapability.GROUND_CUTOUTS, SceneCapability.ELECTRICAL_IDENTITY})


@pytest.mark.parametrize('direct', [False, True])
def test_accounting_is_rejected_before_output(tmp_path: Path, direct: bool, ground_cutout: Any):
    from klayout_pex_protobuf.kpex.pex25d.pex25d_capabilities_pb2 import REQUIRED_CAPABILITY_CAPACITANCE_ACCOUNTING
    original = ground_cutout
    original.required_capabilities.append(REQUIRED_CAPABILITY_CAPACITANCE_ACCOUNTING)
    before = original.SerializeToString(deterministic=True)
    directory = tmp_path/'unsupported'
    with pytest.raises(pex25d.ExportError, match='capacitance-accounting'):
        if direct:
            GmshSceneExporter().export(original, output_dir_path=str(directory))
        else:
            pex25d.export_with_backend(GmshSceneExporter(), original, 'gmsh', str(directory))
    assert not directory.exists()
    assert original.SerializeToString(deterministic=True) == before


def test_mesh_helper_rejects_accounting_before_geometry(tmp_path: Path, ground_cutout: Any):
    from klayout_pex_protobuf.kpex.pex25d.pex25d_capabilities_pb2 import REQUIRED_CAPABILITY_CAPACITANCE_ACCOUNTING
    ground_cutout.required_capabilities.append(REQUIRED_CAPABILITY_CAPACITANCE_ACCOUNTING)
    with patch('klayout_pex_plugins.gmsh.solids.scene_solids',
               side_effect=AssertionError('geometry ran before capability validation')) as build:
        with pytest.raises(pex25d.ExportError, match='capacitance-accounting'):
            write_mesh(ground_cutout, GmshMeshOptions(), str(tmp_path/'unsupported'))
        build.assert_not_called()
    assert not (tmp_path/'unsupported').exists()
