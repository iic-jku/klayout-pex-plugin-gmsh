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

"""The exporter ``gmsh``: the mesh on its own, for any FEM tool or for looking at."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, ClassVar, List, Mapping, Optional

from klayout_pex.plugin_api.v1 import PEX25DSceneExporter

from .mesh import write_mesh
from .options import GmshMeshOptions, parse_options

if TYPE_CHECKING:
    from klayout_pex_protobuf.kpex.pex25d.pex25d_scene_pb2 import PEX25DScene


class GmshSceneExporter(PEX25DSceneExporter):
    """Write a gmsh volume mesh of a resolved PEX25D scene."""

    name: ClassVar[str] = 'gmsh'
    default_prefix: ClassVar[str] = ''

    def export(self,
               scene: PEX25DScene,
               *,
               output_dir_path: str,
               prefix: str = '',
               options: Optional[Mapping[str, Any]] = None) -> List[str]:
        """Return the mesh first, then the JSON naming its physical groups."""
        settings = parse_options(GmshMeshOptions, options, self.name)
        result = write_mesh(scene, settings, output_dir_path, prefix or self.default_prefix)
        return result.paths
