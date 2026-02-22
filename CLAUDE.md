# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MOTIOn** (Modular HPC Optimized Toolkit for IO in Fortran) is a pure Fortran library providing an agnostic, high-level API for efficient parallel I/O in HPC applications. It wraps HDF5 and XDMF file formats for use in MPI-parallel CFD solvers and structured/unstructured grid codes.

## Build System

The project uses **FoBiS.py** (Fortran Building System). All build modes are defined in the `fobos` file at the project root.

```bash
# Build shared library (GNU)
FoBiS.py build -mode shared-gnu

# Build static library (GNU)
FoBiS.py build -mode static-gnu

# Build test executables (GNU)
FoBiS.py build -mode tests-gnu

# Debug variants
FoBiS.py build -mode tests-gnu-debug
FoBiS.py build -mode tests-nvf          # NVIDIA Fortran compiler
FoBiS.py build -mode tests-nvf-debug

# Documentation (requires FORD)
FoBiS.py rule -ex makedoc

# Clean build artifacts
FoBiS.py rule -ex clean
```

Compiled test executables go to `exe/`, object files to `obj/`, `.mod` files to `mod/`.

## Running Tests

Tests require MPI. After building with `tests-gnu`:

```bash
mpirun -np 2 ./exe/motion_write_xdmf_file_test
```

There is no automated test runner — tests are run manually as MPI executables.

## Dependencies

**External:** HDF5 1.14.6 with Fortran and parallel I/O support. Pre-built libraries are included:
- `lib/hdf5/1.14.6/gnu/14.2.0/` — for GNU 14.2+
- `lib/hdf5/1.14.6/nvf/25.3/` — for nvfortran 25.3

**Internal git submodules** (in `src/third_party/`):
- **PENF** — portable Fortran kind parameters (`R8P`, `R4P`, `I8P`, `I4P`, etc.)
- **StringiFor** — string handling
- **FoXy** — XML generation (used for XDMF output)
- **BeFoR64** — Base64 encoding
- **FACE** — ANSI terminal colors

To update submodules: `./scripts/update_submodules.sh`

## Code Architecture

The library is organized in a three-layer hierarchy in `src/lib/`:

### Layer 1 — Base (`motion_file_base_object.F90`)
Abstract type `file_base_object` that handles MPI initialization and exposes `mpi_rank`/`mpi_procs`. All file types extend this.

### Layer 2 — Format-specific objects
- **`motion_hdf5_file_object.F90`** — Type `hdf5_file_object` + `HDF5_PARAMETERS`. Low-level HDF5 file, dataset, and attribute handling. Maps Fortran kinds (R8P…I1P) to HDF5 native types.
- **`motion_xdmf_file_object.F90`** — Type `xdmf_file_object` + `XDMF_PARAMETERS`. XML/XDMF metadata generation via FoXy. Handles grid types, attribute centers, and data item formats.

### Layer 3 — High-level unified API (`motion_xh5f_file_object.F90`)
Type `xh5f_file_object` + `XH5F_PARAMETERS`. Combines HDF5 and XDMF into a single interface. Supports block types: `cartesian`, `cartesian_uniform`, `curvilinear`. This is the primary user-facing API.

### Entry point (`motion.F90`)
Thin wrapper that re-exports all public types and parameter constants from the four modules above.

### Template includes (`.INC` files)
- `motion_hdf5_load_save_dataset_agnostic.INC` — Expands into overloaded procedures for each Fortran kind.
- `motion_xh5f_load_save_block_field_agnostic.INC` — Same pattern for block field I/O.

These are `#include`d inside modules to generate kind-specific overloads without code duplication.

## Preprocessor Flags

The standard build uses:
```
-D_ASCII_SUPPORTED -D_UCS4_SUPPORTED -D_R16P
```
`-D_R16P` enables 128-bit real support via PENF. The NVF templates do not define these (compiler differences).

## Compiler Support

| Compiler     | Status      |
|-------------|-------------|
| gfortran 14.2+ | Supported |
| nvfortran 25.3 | Supported |
| Intel/IBM     | Not tested |

The debug GNU mode uses `-std=f2018 -fimplicit-none -fbacktrace -finit-real=nan` among other strict flags — useful for catching regressions.
