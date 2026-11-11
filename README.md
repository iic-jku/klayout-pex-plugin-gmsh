# klayout-pex-plugin-gmsh

gmsh volume meshes of [KLayout-PEX](https://github.com/iic-jku/klayout-pex) PEX25D scenes, for FEM field solvers:
- the exporter `gmsh`
- the mesh builder that [klayout-pex-plugin-elmer](https://github.com/iic-jku/klayout-pex-plugin-elmer) and [klayout-pex-plugin-palace](https://github.com/iic-jku/klayout-pex-plugin-palace) share

## Usage

```bash
pip install klayout-pex-plugin-gmsh
pex25d export cell.pex25d --to gmsh --out_dir mesh
```

Output:
- `mesh.msh`: gmsh 2.2 ASCII, coordinates in µm, with physical groups
  - one volume per dielectric, the background included
  - one surface per net, floating conductor and ground plane: where it touches a dielectric
  - `outer_boundary`: the faces of the domain box
- `mesh.json`: tag, name and permittivity of every group

Conductors are not meshed, as electrostatics only needs their surfaces.

## Options

Pass them with `pex25d export --option NAME=VALUE` (klayout-pex 0.6.3 or later) or `options={...}` in Python. `--field_margin` sets `field_margin_um`.

| Option | Default | Meaning |
| --- | --- | --- |
| `field_margin_um` | 8.0 | How far the domain reaches beyond the geometry, sideways and above the stack; ignored when the scene has a domain |
| `mesh_size_edge_um` | thinnest conductor layer | Element size at conductor edges, where the field concentrates; elements grow away from them |
| `mesh_size_max_um` | 20 × edge size | Largest element |
| `mesh_size_min_um` | gmsh | Smallest element |
| `verbose` | `False` | Show gmsh's output |

`mesh.json` records the sizes used.

How the default edge size does on the ihp-sg13g2 test design `overlap_plates_100um_x_100um_m1_m2` (LOWER–UPPER, Elmer):

| Edge size | Tetrahedra | C |
| --- | --- | --- |
| 0.84 µm | 54k | 182.5 fF |
| 0.42 µm (default: Metal1) | 94k | 181.3 fF |
| 0.21 µm | ~180k | 180.7 fF |

## How the mesh is built

- **Occupancy:** the domain is sliced at every z where a solid starts or ends. Per slab, conductors and the ground plane claim their regions first, then the dielectrics by wrap depth, then the background. The resulting prisms fill the domain exactly once.
- **Conformal films:** grow along their WRAPS chain, as PEX25D specifies. A film wrapping a film that covers the field lies on top of it there.
- **Conformal mesh:** the prisms are extruded with OpenCASCADE and fragmented, so touching solids share faces.

## Requirements

- klayout-pex 0.6.x (plugin API v1 is provisional)
- gmsh ≥ 4.12. IIC-OSIC-TOOLS provides 4.12.1 from APT, and PyPI has no Linux aarch64 wheel, so the code sticks to the 4.12 API.

## Development

```bash
poetry install
poetry run pytest
```

## Acknowledgements

This project is funded by the JKU/SAL [IWS Lab](https://research.jku.at/de/projects/jku-lit-sal-intelligent-wireless-systems-lab-iws-lab/), a collaboration of [Johannes Kepler University](https://jku.at) and [Silicon Austria Labs](https://silicon-austria-labs.com).

This project is further funded by the European project [ODE4EC-AMS](https://ode4ec.eu/projects/ams/), which has received funding from the Chips Joint Undertaking (Chips JU) under Grant Agreement No 101252962.

Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the Chips Joint Undertaking. Neither the European Union nor the granting authority can be held responsible for them.

<p align="center">
  <a href="https://iic.jku.at" target="_blank">
    <img src="https://github.com/iic-jku/klayout-pex-website/raw/main/figures/funding/iic-jku.svg" alt="Johannes Kepler University: Institute for Integrated Circuits and Quantum Computing" width="300"/>
  </a>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://silicon-austria-labs.com" target="_blank">
    <img src="https://github.com/iic-jku/klayout-pex-website/raw/main/figures/funding/silicon-austria-labs-logo.svg" alt="Silicon Austria Labs" width="300"/>
  </a>
</p>

<p align="center">
  <a href="https://ode4ec.eu/projects/ams/" target="_blank">
    <img src="https://github.com/iic-jku/klayout-pex-website/raw/main/figures/funding/eu-cofunded.svg" alt="Co-funded by the European Union" width="300"/>
  </a>
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <a href="https://www.chips-ju.europa.eu" target="_blank">
    <img src="https://github.com/iic-jku/klayout-pex-website/raw/main/figures/funding/chips-ju.png" alt="Chips Joint Undertaking" width="300"/>
  </a>
</p>
