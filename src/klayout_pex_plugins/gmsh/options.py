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

"""Mesh options, shared by the exporters built on this mesh (Elmer, Palace)."""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any, Mapping, Optional, Type, TypeVar

from klayout_pex.plugin_api.v1 import ExportError

_TYPES = {
    'bool': (bool,),
    'int': (int,),
    'float': (int, float),
    'str': (str,),
}


@dataclass
class GmshMeshOptions:
    """Domain and mesh settings."""

    field_margin_um: float = 8.0
    """
    How far the domain reaches beyond the finite geometry, sideways and above the
    stack. Ignored when the scene carries a domain, which says it outright. Same
    name and default as the FasterCap exporters' option.
    """

    mesh_size_edge_um: Optional[float] = None
    """
    Element size at conductor edges, where the field concentrates; elements grow
    away from them. Unset: the thickness of the thinnest conductor layer.
    """

    mesh_size_max_um: Optional[float] = None
    """Largest element size. Unset: 20 times the size at conductor edges."""

    mesh_size_min_um: Optional[float] = None
    """Smallest element size; unset leaves it to gmsh."""

    verbose: bool = False
    """Show gmsh's own output."""

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            optional = field.type.startswith('Optional[')
            type_name = field.type[len('Optional['):-1] if optional else field.type
            if value is None and optional:
                continue
            allowed = _TYPES[type_name]
            if (isinstance(value, bool) and bool not in allowed) \
                    or not isinstance(value, allowed):
                raise TypeError(f"'{field.name}' must be {type_name}"
                                f"{' or None' if optional else ''}, not {value!r}")
        if self.field_margin_um < 0:
            raise ValueError("'field_margin_um' must not be negative")
        for name in ('mesh_size_edge_um', 'mesh_size_max_um', 'mesh_size_min_um'):
            if getattr(self, name) is not None and getattr(self, name) <= 0:
                raise ValueError(f"'{name}' must be positive")


Options = TypeVar('Options', bound=GmshMeshOptions)


def parse_options(cls: Type[Options], options: Optional[Mapping[str, Any]],
                  exporter_name: str) -> Options:
    """Options as ``cls``; unknown names and wrong types are an ExportError."""
    try:
        return cls(**(options or {}))
    except (TypeError, ValueError) as exc:
        raise ExportError(f"Invalid options for '{exporter_name}': {exc}") from exc
