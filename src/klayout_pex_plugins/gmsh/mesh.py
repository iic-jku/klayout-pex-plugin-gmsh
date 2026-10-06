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
Disjoint prisms → gmsh volume mesh for electrostatics.

The prisms are extruded with OpenCASCADE and fragmented, so that touching solids
share their faces and the mesh is conformal. Conductors are then removed: only
the dielectrics are meshed, and each net is the set of faces where it touches a
dielectric. Physical groups carry the names:

* one volume group per dielectric, the background included
* one surface group per terminal: each net, each floating conductor, the ground plane
* one surface group ``outer_boundary`` for the faces of the domain box

Tags count from 1 within each kind, so they can serve directly as Palace
attributes and Elmer body / boundary indices. Coordinates are in µm.
"""

from __future__ import annotations

import contextlib
from dataclasses import asdict, dataclass, field
from importlib import metadata
import json
import math
import os
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional

from klayout_pex.log import debug, info
from klayout_pex.plugin_api.v1 import ExportError, ExporterUnavailable

if TYPE_CHECKING:
    from klayout_pex_protobuf.kpex.pex25d.pex25d_scene_pb2 import PEX25DScene

    from .options import GmshMeshOptions
    from .solids import Body, SceneSolids

OUTER_BOUNDARY = 'outer_boundary'

LENGTH_UNIT_M = 1e-6
"""Mesh coordinates are in µm."""

MAX_SIZE_PER_EDGE_SIZE = 20.0
"""Default largest element, relative to the size at conductor edges."""

SIZE_GROWTH = 0.5
"""How fast elements grow away from conductor edges: size increase per distance."""

MAX_EDGE_SAMPLES = 2000


@dataclass(frozen=True)
class PhysicalGroup:
    dim: int
    tag: int
    name: str
    permittivity: Optional[float] = None
    """Relative permittivity of a dielectric volume."""


@dataclass
class MeshResult:
    """What was written, and which physical group is what."""

    msh_path: str
    groups_path: str
    dielectrics: List[PhysicalGroup] = field(default_factory=list)
    terminals: List[PhysicalGroup] = field(default_factory=list)
    outer_boundary: Optional[PhysicalGroup] = None
    mesh_size_edge_um: Optional[float] = None
    """Element size used at conductor edges."""

    mesh_size_max_um: Optional[float] = None
    """Largest element size used."""

    node_count: int = 0
    tetrahedron_count: int = 0
    terminal_faces: List[int] = field(default_factory=list, repr=False)
    """The model's faces of all terminals, for refining the mesh near conductors."""

    @property
    def paths(self) -> List[str]:
        return [self.msh_path, self.groups_path]


def write_mesh(scene: PEX25DScene,
               settings: GmshMeshOptions,
               output_dir_path: str,
               prefix: str = '') -> MeshResult:
    """
    Mesh ``scene`` and write ``<prefix>mesh.msh`` (gmsh 2.2, ASCII), plus
    ``<prefix>mesh.json`` naming each physical group.
    """
    try:
        import gmsh
    except ImportError as exc:
        raise ExporterUnavailable(f"The gmsh Python module is not installed: {exc}") from exc
    from .solids import scene_solids

    solids = scene_solids(scene, settings.field_margin_um)

    os.makedirs(output_dir_path, exist_ok=True)
    result = MeshResult(msh_path=os.path.join(output_dir_path, f"{prefix}mesh.msh"),
                        groups_path=os.path.join(output_dir_path, f"{prefix}mesh.json"))
    with _session(gmsh, settings.verbose):
        _build(gmsh, solids, result)
        _generate(gmsh, settings, solids, result)
        gmsh.write(result.msh_path)

    with open(result.groups_path, 'w', encoding='utf-8') as file:
        json.dump(_groups_document(result), file, indent=2)
        file.write('\n')
    info(f"Meshed {len(result.dielectrics)} dielectric(s) and {len(result.terminals)} "
         f"terminal(s): {result.node_count} nodes, {result.tetrahedron_count} tetrahedra, "
         f"{result.mesh_size_edge_um} µm at conductor edges, up to {result.mesh_size_max_um} µm")
    return result


@contextlib.contextmanager
def _session(gmsh: Any, verbose: bool) -> Iterator[None]:
    """A model of our own; gmsh itself is finalized only if it was not running."""
    owner = not gmsh.isInitialized()
    if owner:
        gmsh.initialize(readConfigFiles=False)
    gmsh.option.setNumber('General.Terminal', 1 if verbose else 0)
    gmsh.model.add('klayout-pex-plugin-gmsh')
    try:
        yield
    finally:
        gmsh.model.remove()
        if owner:
            gmsh.finalize()


def _build(gmsh: Any, solids: SceneSolids, result: MeshResult):
    from .solids import BodyKind

    occ = gmsh.model.occ
    body_of_volume: Dict[int, Body] = {}
    for prism in solids.prisms:
        z = solids.um(prism.zlow)
        height = solids.um(prism.zhigh - prism.zlow)
        for polygon in prism.region.each():
            for volume in _extrude(occ, polygon, z, height, solids.grid_um):
                body_of_volume[volume] = prism.body
    debug(f"{len(body_of_volume)} volumes before fragmenting")

    inputs = [(3, volume) for volume in body_of_volume]
    if len(inputs) > 1:
        # The prisms are disjoint, so each fragment comes from exactly one of them.
        _, children = occ.fragment(inputs[:1], inputs[1:])
        fragments: Dict[int, Body] = {}
        for (_, parent), parent_children in zip(inputs, children):
            for _, child in parent_children:
                if fragments.setdefault(child, body_of_volume[parent]) != body_of_volume[parent]:
                    raise ExportError(f"Solids of '{fragments[child].name}' and "
                                      f"'{body_of_volume[parent].name}' overlap")
        body_of_volume = fragments
    occ.synchronize()

    terminal_faces: Dict[Body, List[int]] = {}
    outer_faces: List[int] = []
    for _, face in gmsh.model.getEntities(2):
        volumes, _ = gmsh.model.getAdjacencies(2, face)
        bodies = [body_of_volume[volume] for volume in volumes]
        conductors = [body for body in bodies if body.kind == BodyKind.CONDUCTOR]
        if len(bodies) == 1 and not conductors:
            outer_faces.append(face)
        elif len(bodies) == 2 and len(conductors) == 1:
            terminal_faces.setdefault(conductors[0], []).append(face)
        elif len(conductors) == 2 and conductors[0] != conductors[1]:
            raise ExportError(f"Conductors '{conductors[0].name}' and "
                              f"'{conductors[1].name}' touch, which shorts them")

    # Only the dielectrics are meshed. Removing the conductors keeps their faces
    # as the dielectrics' boundaries; faces left without a volume are dropped.
    occ.remove([(3, volume) for volume, body in body_of_volume.items()
                if body.kind == BodyKind.CONDUCTOR])
    occ.synchronize()
    dangling = [(2, face) for _, face in gmsh.model.getEntities(2)
                if len(gmsh.model.getAdjacencies(2, face)[0]) == 0]
    if dangling:
        occ.remove(dangling, recursive=True)
        occ.synchronize()

    for body in solids.bodies:
        volumes = [volume for volume, owner in body_of_volume.items()
                   if owner == body and body.kind == BodyKind.DIELECTRIC]
        if volumes:
            result.dielectrics.append(_physical_group(
                gmsh, 3, volumes, len(result.dielectrics) + 1, body.name, body.permittivity))
    for body in solids.bodies:
        if terminal_faces.get(body):
            result.terminals.append(_physical_group(
                gmsh, 2, terminal_faces[body], len(result.terminals) + 1, body.name))
            result.terminal_faces.extend(terminal_faces[body])
    if any(group.name == OUTER_BOUNDARY for group in result.terminals):
        raise ExportError(f"A net is named '{OUTER_BOUNDARY}', the name of the "
                          f"domain boundary group")
    if outer_faces:
        result.outer_boundary = _physical_group(
            gmsh, 2, outer_faces, len(result.terminals) + 1, OUTER_BOUNDARY)


def _extrude(occ: Any, polygon: Any, z: float, height: float, scale: float) -> List[int]:
    loops = [_curve_loop(occ, polygon.each_point_hull(), z, scale)]
    loops.extend(_curve_loop(occ, polygon.each_point_hole(hole), z, scale)
                 for hole in range(polygon.holes()))
    surface = occ.addPlaneSurface(loops)
    return [tag for dim, tag in occ.extrude([(2, surface)], 0, 0, height) if dim == 3]


def _curve_loop(occ: Any, points: Any, z: float, scale: float) -> int:
    # gmsh wants the holes wound like the hull, KLayout winds them the other
    # way; opposite windings make OCC add a hole's area instead of cutting it.
    points = [(point.x, point.y) for point in points]
    twice_area = sum(x0 * y1 - x1 * y0
                     for (x0, y0), (x1, y1) in zip(points, points[1:] + points[:1]))
    if twice_area < 0:
        points.reverse()
    tags = [occ.addPoint(x * scale, y * scale, z) for x, y in points]
    lines = [occ.addLine(start, end) for start, end in zip(tags, tags[1:] + tags[:1])]
    return occ.addCurveLoop(lines)


def _physical_group(gmsh: Any, dim: int, entities: List[int], tag: int, name: str,
                    permittivity: Optional[float] = None) -> PhysicalGroup:
    gmsh.model.addPhysicalGroup(dim, entities, tag)
    gmsh.model.setPhysicalName(dim, tag, name)
    return PhysicalGroup(dim, tag, name, permittivity)


def _generate(gmsh: Any, settings: GmshMeshOptions, solids: SceneSolids, result: MeshResult):
    gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
    gmsh.option.setNumber('Mesh.Binary', 0)
    gmsh.option.setNumber('Mesh.SaveAll', 0)

    edge_size = settings.mesh_size_edge_um
    if edge_size is None and solids.conductor_thickness:
        edge_size = solids.um(solids.conductor_thickness)
    max_size = settings.mesh_size_max_um
    if max_size is None and edge_size is not None:
        max_size = edge_size * MAX_SIZE_PER_EDGE_SIZE
    if edge_size is not None and result.terminal_faces:
        _refine_at_conductor_edges(gmsh, result.terminal_faces, edge_size, max_size)
        result.mesh_size_edge_um = edge_size
    result.mesh_size_max_um = max_size

    if settings.mesh_size_min_um is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMin', settings.mesh_size_min_um)
    if max_size is not None:
        gmsh.option.setNumber('Mesh.MeshSizeMax', max_size)
    gmsh.model.mesh.generate(3)

    result.node_count = len(gmsh.model.mesh.getNodes()[0])
    element_types, element_tags, _ = gmsh.model.mesh.getElements(3)
    result.tetrahedron_count = sum(len(tags) for element_type, tags
                                   in zip(element_types, element_tags) if element_type == 4)


def _refine_at_conductor_edges(gmsh: Any, faces: List[int], edge_size: float,
                               max_size: float):
    """
    The field concentrates at conductor edges: elements are ``edge_size`` there and
    grow with the distance from the nearest edge, up to ``max_size``.
    """
    curves = sorted({curve for _, curve in gmsh.model.getBoundary(
        [(2, face) for face in faces], combined=False, oriented=False)})
    longest = max(gmsh.model.occ.getMass(1, curve) for curve in curves)

    fields = gmsh.model.mesh.field
    distance = fields.add('Distance')
    fields.setNumbers(distance, 'CurvesList', curves)
    # Sample the longest edge at the element size, so the distance is not bumpy.
    fields.setNumber(distance, 'Sampling', min(MAX_EDGE_SAMPLES, max(
        2, math.ceil(longest / edge_size) + 1)))
    threshold = fields.add('Threshold')
    fields.setNumber(threshold, 'InField', distance)
    fields.setNumber(threshold, 'SizeMin', edge_size)
    fields.setNumber(threshold, 'SizeMax', max_size)
    fields.setNumber(threshold, 'DistMin', 0.0)
    fields.setNumber(threshold, 'DistMax', (max_size - edge_size) / SIZE_GROWTH)
    fields.setAsBackgroundMesh(threshold)

    # The field alone decides; sizes from points or boundaries would override it.
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints', 0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature', 0)


def _groups_document(result: MeshResult) -> Dict[str, Any]:
    def entry(group: PhysicalGroup) -> Dict[str, Any]:
        document = asdict(group)
        document.pop('dim')
        if document['permittivity'] is None:
            document.pop('permittivity')
        return document

    return {
        'generator': f"klayout-pex-plugin-gmsh {metadata.version('klayout-pex-plugin-gmsh')}",
        'mesh': os.path.basename(result.msh_path),
        'length_unit_m': LENGTH_UNIT_M,
        'mesh_size_edge_um': result.mesh_size_edge_um,
        'mesh_size_max_um': result.mesh_size_max_um,
        'dielectrics': [entry(group) for group in result.dielectrics],
        'terminals': [entry(group) for group in result.terminals],
        'outer_boundary': entry(result.outer_boundary) if result.outer_boundary else None,
    }
