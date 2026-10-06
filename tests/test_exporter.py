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

"""The mesh, its physical groups, and the exporter as KLayout-PEX sees it."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import gmsh
import pytest

from klayout_pex import pex25d
from klayout_pex.pex25d.pex25d_cli import Pex25DCLI
from klayout_pex_plugins.gmsh import GmshMeshOptions, write_mesh
from klayout_pex_plugins.gmsh.solids import BodyKind, scene_solids

from .conftest import DATA


COARSE = {'mesh_size_edge_um': 5.0, 'mesh_size_max_um': 40.0}
"""Enough for checking groups and volumes, and fast."""


def physical_groups(msh_path: str):
    """{(dim, name): element count} as read back from the file."""
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.open(msh_path)
        groups = {}
        for dim, tag in gmsh.model.getPhysicalGroups():
            count = 0
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                _, element_tags, _ = gmsh.model.mesh.getElements(dim, entity)
                count += sum(len(tags) for tags in element_tags)
            groups[(dim, gmsh.model.getPhysicalName(dim, tag))] = count
        return groups
    finally:
        gmsh.finalize()


def meshed_volumes(msh_path: str) -> Dict[str, float]:
    """Volume of each dielectric group as read back from the file, in µm³."""
    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.open(msh_path)
        volumes = {}
        for dim, tag in gmsh.model.getPhysicalGroups(3):
            volume = 0.0
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                _, element_tags, _ = gmsh.model.mesh.getElements(dim, entity)
                for tags in element_tags:
                    volume += sum(gmsh.model.mesh.getElementQualities(tags, 'volume'))
            volumes[gmsh.model.getPhysicalName(dim, tag)] = volume
        return volumes
    finally:
        gmsh.finalize()


@pytest.mark.parametrize('scene', ['two_nets', 'overlap_plates'])
def test_meshed_volume_of_each_dielectric_matches_its_prisms(
        scene: str, request: pytest.FixtureRequest, tmp_path: Path):
    scene = request.getfixturevalue(scene)
    solids = scene_solids(scene, field_margin_um=8.0)
    expected: Dict[str, float] = {}
    for prism in solids.prisms:
        if prism.body.kind == BodyKind.DIELECTRIC:
            expected[prism.body.name] = expected.get(prism.body.name, 0.0) \
                + prism.region.area() * (prism.zhigh - prism.zlow) * solids.grid_um ** 3
    result = write_mesh(scene, GmshMeshOptions(**COARSE), str(tmp_path))
    assert meshed_volumes(result.msh_path) == pytest.approx(expected, rel=1e-9)


def test_default_sizes_come_from_the_thinnest_conductor_layer(two_nets, tmp_path: Path):
    result = write_mesh(two_nets, GmshMeshOptions(), str(tmp_path))
    assert result.mesh_size_edge_um == pytest.approx(0.4)      # met1 and met2
    assert result.mesh_size_max_um == pytest.approx(8.0)
    document = json.loads(Path(result.groups_path).read_text())
    assert (document['mesh_size_edge_um'], document['mesh_size_max_um']) == \
        (result.mesh_size_edge_um, result.mesh_size_max_um)


def test_finer_edges_mean_more_elements(two_nets, tmp_path: Path):
    coarse = write_mesh(two_nets, GmshMeshOptions(), str(tmp_path / 'coarse'))
    fine = write_mesh(two_nets, GmshMeshOptions(mesh_size_edge_um=0.1), str(tmp_path / 'fine'))
    assert fine.tetrahedron_count > 2 * coarse.tetrahedron_count


def test_mesh_has_dielectric_volumes_and_terminal_surfaces(two_nets, tmp_path: Path):
    result = write_mesh(two_nets, GmshMeshOptions(), str(tmp_path))

    assert [group.name for group in result.dielectrics] == ['fox', 'lint', 'cap', 'ild', 'air']
    assert [group.permittivity for group in result.dielectrics] == [3.9, 7.3, 5.0, 4.1, 1.0]
    assert [group.name for group in result.terminals] == ['subs', 'neta', 'netb', 'F']
    assert [group.tag for group in result.terminals] == [1, 2, 3, 4]
    assert result.outer_boundary.tag == 5

    groups = physical_groups(result.msh_path)
    expected = {(3, group.name) for group in result.dielectrics} \
        | {(2, group.name) for group in result.terminals} | {(2, 'outer_boundary')}
    assert set(groups) == expected
    assert all(count > 0 for count in groups.values())


def test_groups_document_names_every_tag(two_nets, tmp_path: Path):
    result = write_mesh(two_nets, GmshMeshOptions(), str(tmp_path), prefix='tiny_')
    assert result.paths == [str(tmp_path / 'tiny_mesh.msh'), str(tmp_path / 'tiny_mesh.json')]
    document = json.loads(Path(result.groups_path).read_text())
    assert document['mesh'] == 'tiny_mesh.msh'
    assert document['length_unit_m'] == 1e-6
    assert document['terminals'][1] == {'tag': 2, 'name': 'neta'}
    assert document['dielectrics'][0] == {'tag': 1, 'name': 'fox', 'permittivity': 3.9}
    assert document['outer_boundary'] == {'tag': 5, 'name': 'outer_boundary'}


def test_touching_nets_are_rejected(tmp_path: Path):
    text = (DATA / 'two_nets.pex25d').read_text().replace(
        'OUTER 2.0 0.0 3.0 0.0 3.0 1.0 2.0 1.0', 'OUTER 1.0 0.0 3.0 0.0 3.0 1.0 1.0 1.0')
    scene = pex25d.resolve(pex25d.read_text(text.encode()))
    with pytest.raises(pex25d.ExportError, match="'neta' and 'netb' touch"):
        write_mesh(scene, GmshMeshOptions(), str(tmp_path))


def test_registered_as_exporter_plugin():
    info = pex25d.exporter_registry().get_info('gmsh')
    assert info.distribution == 'klayout-pex-plugin-gmsh'
    assert info.target == 'klayout_pex_plugins.gmsh:create_exporter'


def test_public_api_export(overlap_plates, tmp_path: Path):
    written = pex25d.export(overlap_plates, 'klayout-pex-plugin-gmsh:gmsh', str(tmp_path),
                            options={'field_margin_um': 4.0, **COARSE})
    assert written == [str(tmp_path / 'mesh.msh'), str(tmp_path / 'mesh.json')]
    groups = physical_groups(written[0])
    assert {name for dim, name in groups if dim == 2} == {'subs', 'LOWER', 'UPPER',
                                                          'outer_boundary'}


@pytest.mark.parametrize('options, message', [
    ({'delaunay_b': 1.0}, "unexpected keyword argument 'delaunay_b'"),
    ({'mesh_size_max_um': 'fine'}, "'mesh_size_max_um' must be float or None"),
    ({'verbose': 1}, "'verbose' must be bool"),
    ({'field_margin_um': -1.0}, "must not be negative"),
], ids=['unknown', 'type', 'bool', 'value'])
def test_invalid_options_are_export_errors(two_nets, tmp_path: Path, options, message: str):
    with pytest.raises(pex25d.ExportError, match=message):
        pex25d.export(two_nets, 'gmsh', str(tmp_path), options=options)


def test_cli(tmp_path: Path):
    with pytest.raises(SystemExit) as completion:
        Pex25DCLI().main(['pex25d', 'export', str(DATA / 'two_nets.pex25d'),
                          '--to', 'gmsh', '--out_dir', str(tmp_path)])
    assert completion.value.code == pex25d.ExitCode.OK
    assert sorted(path.name for path in tmp_path.iterdir()) == ['mesh.json', 'mesh.msh']
