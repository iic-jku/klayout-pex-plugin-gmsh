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
gmsh volume meshes of PEX25D scenes: the KLayout-PEX exporter ``gmsh``, and the
mesh builder the FEM exporters (Elmer, Palace) share.

gmsh and KLayout load on export, so plugin discovery works without them.
"""

from .exporter import GmshSceneExporter
from .mesh import LENGTH_UNIT_M, OUTER_BOUNDARY, MeshResult, PhysicalGroup, write_mesh
from .options import GmshMeshOptions, parse_options


def create_exporter() -> GmshSceneExporter:
    return GmshSceneExporter()


__all__ = [
    'GmshMeshOptions',
    'GmshSceneExporter',
    'LENGTH_UNIT_M',
    'MeshResult',
    'OUTER_BOUNDARY',
    'PhysicalGroup',
    'create_exporter',
    'parse_options',
    'write_mesh',
]
