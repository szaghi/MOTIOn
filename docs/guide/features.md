---
title: Features
---

# Features

## High-level API — XH5F

`xh5f_file_object` is the primary user-facing type. It coordinates an HDF5 data file and an XDMF metadata file through a single, unified interface. MPI parallelization is handled transparently: each rank writes its own HDF5 file while the master rank assembles the final XDMF descriptor via `MPI_GATHERV`.

### Grid topologies

| Block type constant | Description |
|---------------------|-------------|
| `XH5F_BLOCK_CARTESIAN_UNIFORM` | Uniform spacing — described by origin `emin(3)` and step `dxyz(3)` |
| `XH5F_BLOCK_CARTESIAN` | Variable spacing — described by per-axis node arrays `x`, `y`, `z` |
| `XH5F_BLOCK_CURVILINEAR` | General curvilinear mesh — described by full `nodes(3,0:ni,0:nj,0:nk)` array |

### Data fields

All field operations are overloaded for ranks and kinds:

| Field shape | Leading dimension of `nd` |
|-------------|--------------------------|
| Scalar 0D (grid-centered) | — (no `nd` argument) |
| Scalar 3D | `nd = [ni, nj, nk]` |
| Vector 3D | `nd = [3, ni, nj, nk]` |
| Tensor6 3D | `nd = [6, ni, nj, nk]` |
| Tensor 3D | `nd = [9, ni, nj, nk]` |
| Matrix 3D | `nd = [m, ni, nj, nk]` |

Supported Fortran kinds: `R8P`, `R4P`, `I8P`, `I4P`, `I2P`, `I1P` (via [PENF](https://github.com/szaghi/PENF)).

### Field centering

| Constant | Description |
|----------|-------------|
| `XDMF_ATTR_CENTER_CELL` | Cell-centered (supported) |
| `XDMF_ATTR_CENTER_GRID` | Grid-centered / 0D (supported) |
| `XDMF_ATTR_CENTER_NODE` | Node-centered (not yet tested) |
| `XDMF_ATTR_CENTER_EDGE` | Edge-centered (not yet tested) |
| `XDMF_ATTR_CENTER_FACE` | Face-centered (not yet tested) |

### DataItem format

| Constant | Description |
|----------|-------------|
| `XDMF_DATAITEM_NUMBER_FORMAT_HDF` | Data written to HDF5 file, path referenced in XDMF |
| `XDMF_DATAITEM_NUMBER_FORMAT_XML` | Data inlined in the XDMF file (suitable for small scalars) |

---

## Low-level API — HDF5

`hdf5_file_object` provides direct access to HDF5 file and dataset operations.

- Open / close HDF5 files with `open_file` / `close_file`
- Open / close simple dataspaces with `open_dspace` / `close_dspace`
- Save datasets with `save_dataset` (overloaded for ranks 0–7 and all PENF kinds)
- Load datasets with `load_dataset` (same overloads)
- Introspect: `does_dataset_exist`, `get_dataset_dims`

---

## Low-level API — XDMF

`xdmf_file_object` generates XDMF XML via [FoXy](https://github.com/szaghi/FoXy).

### Grid types

| Constant | Description |
|----------|-------------|
| `XDMF_GRID_TYPE_UNIFORM` | Single grid |
| `XDMF_GRID_TYPE_COLLECTION` | Synchronous collection |
| `XDMF_GRID_TYPE_COLLECTION_ASYNC` | Asynchronous collection (MPI use case) |
| `XDMF_GRID_TYPE_TREE` | Tree hierarchy |
| `XDMF_GRID_TYPE_SUBSET` | Subset of another grid |

### 3-D structured topologies

| Constant | Fortran / XDMF name |
|----------|---------------------|
| `XDMF_TOPOLOGY_TYPE_3DCORECTMESH` | Cartesian uniform |
| `XDMF_TOPOLOGY_TYPE_3DRECTMESH` | Cartesian (rectilinear) |
| `XDMF_TOPOLOGY_TYPE_3DSMESH` | Curvilinear |

---

## Compiler Support

| Compiler | Status |
|----------|--------|
| GNU gfortran ≥ 14.2.0 | Supported |
| NVIDIA nvfortran 25.3 | Supported |
| Intel Fortran | Not tested |
| IBM XL Fortran | Not tested |

## Design Principles

- **Pure Fortran** — no C extensions beyond the HDF5 Fortran interface
- **OOP** — all functionality exposed as type-bound procedures
- **KISS** — the high-level XH5F API intentionally exposes a minimal surface
- **Agnostic kinds** — the same method names work for all real/integer precisions
- **Free & Open Source** — multi-licensed for both FOSS and commercial use
