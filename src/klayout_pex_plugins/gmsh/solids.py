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
Resolved PEX25D scene → disjoint 2.5D prisms, one body per net and per dielectric.

Every PEX25D solid is a 2D region extruded over a z range, so the whole scene is a
set of prisms. Where several solids claim a point, the occupancy rule decides:
conductors and the ground plane first, then the dielectrics by ascending wrap
depth (the order a resolved scene lists them in), the background last. Slicing
the domain at every z where a solid starts or ends makes that rule a sequence of
2D booleans per slab. The result covers the domain exactly once, which is what a
volume mesher needs.

Conformal films grow incrementally along their WRAPS chain, as the format
specifies: a film's solid is the solid it wraps, grown sideways and upward, plus
a slab over the field where the wrapped solid is absent. A film that wraps a film
covering the field therefore lies on top of it there.

Coordinates stay in integer grid units; :attr:`SceneSolids.grid_um` converts.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from enum import Enum
from fractions import Fraction
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

import klayout.db as kdb

from klayout_pex import pex25d
from klayout_pex.plugin_api.v1 import ExportError

if TYPE_CHECKING:
    from klayout_pex_protobuf.kpex.pex25d.pex25d_scene_pb2 import PEX25DScene


class BodyKind(Enum):
    CONDUCTOR = 'conductor'
    DIELECTRIC = 'dielectric'


@dataclass(frozen=True)
class Body:
    """One net (or floating conductor, or the ground plane), or one dielectric."""

    name: str
    kind: BodyKind
    permittivity: Optional[float] = None
    """Relative permittivity; dielectrics only."""


@dataclass
class Prism:
    """A region of one body extruded over [zlow, zhigh], in grid units."""

    body: Body
    region: kdb.Region
    zlow: int
    zhigh: int


@dataclass
class SceneSolids:
    """The domain split into disjoint prisms; their union is ``field`` × [zlow, zhigh]."""

    grid_um: float
    field: kdb.Box
    zlow: int
    zhigh: int
    bodies: List[Body]
    """Conductors first, then dielectrics in occupancy order, the background last."""

    prisms: List[Prism] = dataclass_field(default_factory=list)

    conductor_thickness: Optional[int] = None
    """The thinnest layer that carries conductor geometry, the ground plane aside."""

    def um(self, grid_units: int) -> float:
        return grid_units * self.grid_um


Claim = Tuple[kdb.Region, int, int]


def length_unit_um(units) -> Fraction:
    unit = pex25d.proto.units.Units
    factors = {
        unit.LENGTH_UNIT_UM: Fraction(1),
        unit.LENGTH_UNIT_NM: Fraction(1, 1000),
        unit.LENGTH_UNIT_M: Fraction(1_000_000),
    }
    if units.length not in factors:
        raise ExportError("The scene declares no length unit")
    return factors[units.length]


def scene_solids(scene: PEX25DScene, field_margin_um: float) -> SceneSolids:
    """Split the domain of ``scene`` into disjoint prisms by the occupancy rule."""
    return _SolidsBuilder(scene, field_margin_um).build()


class _SolidsBuilder:
    def __init__(self, scene: PEX25DScene, field_margin_um: float):
        self.scene = scene
        if not scene.units.grid_denominator:
            raise ExportError("The scene declares no grid")
        grid = Fraction(scene.units.grid_numerator, scene.units.grid_denominator)
        self.grid_um = float(grid * length_unit_um(scene.units))
        self.field_margin = int(round(field_margin_um / self.grid_um))

        self.layers = {layer.name: layer for layer in scene.layers}
        self.dielectrics = {d.name: d for d in scene.dielectrics}
        self.field = self._field_box()
        self._solids: Dict[str, List[Claim]] = {}

    # ---------------------------------------------------------------- field

    def _field_box(self) -> kdb.Box:
        """
        The lateral domain: the scene's domain box, or the finite geometry grown by
        the field margin, as the FasterCap exporter does.
        """
        domain = self.scene.domain
        if self.scene.HasField('domain') and domain.origin != domain.ORIGIN_UNSPECIFIED:
            box = domain.box
            return kdb.Box(box.lower_left.x, box.lower_left.y,
                           box.upper_right.x, box.upper_right.y)
        bounds = domain.geometry_bounds
        box = kdb.Box(bounds.lower_left.x, bounds.lower_left.y,
                      bounds.upper_right.x, bounds.upper_right.y)
        if box.empty() or box.width() == 0 or box.height() == 0:
            raise ExportError("The scene has neither a domain nor finite geometry, so "
                              "there is nothing to bound the unbounded materials with")
        return box.enlarged(self.field_margin, self.field_margin)

    @property
    def field_region(self) -> kdb.Region:
        return kdb.Region(self.field)

    # --------------------------------------------------------------- claims

    def build(self) -> SceneSolids:
        claims: List[Tuple[Body, List[Claim]]] = []
        claims.extend(self._conductor_claims())
        for dielectric in self.scene.dielectrics:
            body = Body(dielectric.name, BodyKind.DIELECTRIC, dielectric.permittivity)
            claims.append((body, [(region & self.field_region, zlow, zhigh)
                                  for region, zlow, zhigh in self._solid(dielectric.name)]))

        zlow, zhigh = self._z_extent(claims)
        if not self.scene.HasField('background'):
            raise ExportError("The scene declares no background dielectric")
        background = self.scene.background
        claims.append((Body(background.name, BodyKind.DIELECTRIC, background.permittivity),
                       [(self.field_region, zlow, zhigh)]))

        solids = SceneSolids(grid_um=self.grid_um, field=self.field, zlow=zlow, zhigh=zhigh,
                             bodies=[body for body, _ in claims])
        thicknesses = [self.layers[conductor_region.layer].zhigh
                       - self.layers[conductor_region.layer].zlow
                       for conductor in self.scene.conductors
                       for conductor_region in conductor.regions]
        solids.conductor_thickness = min(thicknesses, default=None)
        solids.prisms = self._disjoint_prisms(claims, zlow, zhigh)
        return solids

    def _conductor_claims(self) -> List[Tuple[Body, List[Claim]]]:
        by_name: Dict[str, List[Claim]] = {}
        if self.scene.HasField('ground_plane'):
            plane = self.scene.ground_plane
            by_name[plane.name] = [(self.field_region, plane.zlow, plane.zhigh)]
        for conductor in self.scene.conductors:
            # A FLOATING conductor belongs to no net and is its own body.
            name = conductor.name if conductor.floating else conductor.net
            for conductor_region in conductor.regions:
                layer = self.layers.get(conductor_region.layer)
                if layer is None:
                    raise ExportError(f"Conductor '{conductor.name}' has geometry on "
                                      f"'{conductor_region.layer}', which the scene "
                                      f"does not declare")
                by_name.setdefault(name, []).append(
                    (region_of(conductor_region), layer.zlow, layer.zhigh))
        return [(Body(name, BodyKind.CONDUCTOR), claims) for name, claims in by_name.items()]

    def _z_extent(self, claims: List[Tuple[Body, List[Claim]]]) -> Tuple[int, int]:
        """
        The whole stack, extended upward to the domain's top, or by the field margin
        when the scene has no domain, so that the background lies above the stack.
        """
        zs = [z for _, body_claims in claims for _, zlow, zhigh in body_claims
              for z in (zlow, zhigh)]
        if not zs:
            raise ExportError("The scene has no solids")
        zlow, zhigh = min(zs), max(zs)
        domain = self.scene.domain
        if self.scene.HasField('domain') and domain.origin != domain.ORIGIN_UNSPECIFIED:
            return zlow, max(zhigh, domain.box.upper_right.z)
        return zlow, zhigh + self.field_margin

    # ------------------------------------------------------------ profiles

    def _solid(self, name: str) -> List[Claim]:
        """The solid of a profile as overlapping prisms, before occupancy."""
        if name not in self._solids:
            self._solids[name] = self._compute_solid(name)
        return self._solids[name]

    def _compute_solid(self, name: str) -> List[Claim]:
        if name in self.layers:
            layer = self.layers[name]
            shapes = self._shapes_of_layer(name)
            return [] if shapes.is_empty() else [(shapes, layer.zlow, layer.zhigh)]

        plane = self.scene.ground_plane
        if self.scene.HasField('ground_plane') and name == plane.name:
            return [(self.field_region, plane.zlow, plane.zhigh)]

        dielectric = self.dielectrics.get(name)
        if dielectric is None:
            raise ExportError(f"'{name}' is not a profile of the scene")

        kinds = pex25d.proto.dielectric
        if dielectric.kind == kinds.DIELECTRIC_KIND_SIMPLE:
            # A band spans the whole plane; the films it encloses win by depth.
            if dielectric.zhigh <= dielectric.zlow:
                return []
            return [(self.field_region, dielectric.zlow, dielectric.zhigh)]
        if dielectric.kind != kinds.DIELECTRIC_KIND_CONFORMAL:
            raise ExportError(f"Dielectric '{name}' has no kind")

        wrapped = self._solid(dielectric.wraps)
        lateral = dielectric.thickness_beside_wrapped
        solid = [(region.sized(lateral) if lateral else region, zlow,
                  zhigh + dielectric.thickness_over_wrapped)
                 for region, zlow, zhigh in wrapped]
        if dielectric.thickness_on_field > 0:
            footprint = kdb.Region()
            for region, _, _ in wrapped:
                footprint += region
            bottom = self._bottom(dielectric.wraps)
            solid.append((self.field_region - footprint, bottom,
                          bottom + dielectric.thickness_on_field))
        return solid

    def _bottom(self, name: str) -> int:
        """The bottom face of a profile, which a film on the field is measured from."""
        if name in self.layers:
            return self.layers[name].zlow
        if name in self.dielectrics:
            return self.dielectrics[name].zlow
        return self.scene.ground_plane.zlow

    def _shapes_of_layer(self, layer_name: str) -> kdb.Region:
        """Every conductor's geometry on one layer, unioned across nets."""
        region = kdb.Region()
        for conductor in self.scene.conductors:
            for conductor_region in conductor.regions:
                if conductor_region.layer == layer_name:
                    region += region_of(conductor_region)
        return region.merged()

    # ----------------------------------------------------------- occupancy

    def _disjoint_prisms(self, claims: List[Tuple[Body, List[Claim]]],
                         zlow: int, zhigh: int) -> List[Prism]:
        """
        Apply the occupancy rule slab by slab, then merge each body's slabs back
        into taller prisms wherever its region does not change.
        """
        zs = sorted({zlow, zhigh} | {z for _, body_claims in claims
                                     for _, z0, z1 in body_claims if z0 < z1
                                     for z in (z0, z1) if zlow < z < zhigh})
        open_runs: Dict[Body, Prism] = {}
        prisms: List[Prism] = []
        for z0, z1 in zip(zs, zs[1:]):
            claimed = kdb.Region()
            conductors = kdb.Region()
            for body, body_claims in claims:
                region = kdb.Region()
                for claim_region, claim_zlow, claim_zhigh in body_claims:
                    if claim_zlow <= z0 and claim_zhigh >= z1:
                        region += claim_region
                region.merge()
                if body.kind == BodyKind.CONDUCTOR:
                    short = region & conductors
                    if not short.is_empty():
                        raise ExportError(
                            f"Conductor '{body.name}' overlaps another net "
                            f"at z {z0}…{z1} (grid units)")
                    conductors += region
                region -= claimed
                claimed += region

                run = open_runs.get(body)
                if run is not None and run.zhigh == z0 and (run.region ^ region).is_empty():
                    run.zhigh = z1
                    continue
                if run is not None:
                    prisms.append(open_runs.pop(body))
                if not region.is_empty():
                    open_runs[body] = Prism(body, region, z0, z1)
        prisms.extend(open_runs.values())
        return prisms


def region_of(conductor_region) -> kdb.Region:
    """Build a region from one conductor's shapes on one layer, in grid units."""
    region = kdb.Region()
    for box in conductor_region.boxes:
        region.insert(kdb.Box(box.lower_left.x, box.lower_left.y,
                              box.upper_right.x, box.upper_right.y))
    for polygon in conductor_region.polygons:
        shape = kdb.Polygon([kdb.Point(p.x, p.y) for p in polygon.outer.points])
        for hole in polygon.holes:
            shape.insert_hole([kdb.Point(p.x, p.y) for p in hole.points])
        region.insert(shape)
    return region.merged()
