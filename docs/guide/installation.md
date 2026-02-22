---
title: Installation
---

# Installation

## Prerequisites

| Requirement | Minimum version |
|-------------|----------------|
| GNU gfortran | 14.2.0 |
| NVIDIA nvfortran | 25.3 |
| MPI implementation | any standard |
| HDF5 (Fortran + parallel) | 1.14.6 |

Pre-built HDF5 1.14.6 libraries are included in the repository under `lib/hdf5/1.14.6/`, so no separate HDF5 build is required for the default GNU and NVIDIA toolchains.

## Download

MOTIOn uses git submodules for its third-party dependencies. Clone recursively:

```bash
git clone https://github.com/szaghi/MOTIOn --recursive
cd MOTIOn
```

If you already have a non-recursive clone:

```bash
./scripts/update_submodules.sh
# or equivalently:
git submodule update --init --recursive
```

### Third-party dependencies

The submodules live under `src/third_party/`:

| Library | Purpose |
|---------|---------|
| [PENF](https://github.com/szaghi/PENF) | Portable numeric kind parameters (`R8P`, `I4P`, …) |
| [StringiFor](https://github.com/szaghi/StringiFor) | String utilities |
| [FoXy](https://github.com/szaghi/FoXy) | Fortran XML generation (used by the XDMF layer) |
| [BeFoR64](https://github.com/szaghi/BeFoR64) | Base64 encoding |
| [FACE](https://github.com/szaghi/FACE) | ANSI terminal colors |

## Build with FoBiS.py

[FoBiS.py](https://github.com/szaghi/FoBiS) is the build system. Install it with pip if not already present:

```bash
pip install FoBiS.py
```

### List all build modes

```bash
FoBiS.py build -lmodes
```

### Build the library

```bash
# Shared library — GNU gfortran (release)
FoBiS.py build -mode shared-gnu

# Static library — GNU gfortran (release)
FoBiS.py build -mode static-gnu

# Shared library — GNU gfortran (debug)
FoBiS.py build -mode shared-gnu-debug

# Static library — GNU gfortran (debug)
FoBiS.py build -mode static-gnu-debug
```

Outputs are placed in `./shared/` or `./static/` respectively.

### Build and run the test suite

```bash
# Build tests (GNU, release)
FoBiS.py build -mode tests-gnu

# Build tests (GNU, debug — strict flags, NaN-initialised reals)
FoBiS.py build -mode tests-gnu-debug

# Build tests (NVIDIA Fortran, release)
FoBiS.py build -mode tests-nvf

# Run all tests
./scripts/run_tests.sh          # uses 2 MPI ranks by default
NP=4 ./scripts/run_tests.sh    # override rank count
```

Compiled test executables are placed in `./exe/`.

### NVIDIA Fortran

Switch to the `nvf` modes, which use the pre-built HDF5 under `lib/hdf5/1.14.6/nvf/25.3/`:

```bash
FoBiS.py build -mode tests-nvf
FoBiS.py build -mode tests-nvf-debug
```

### Documentation

```bash
FoBiS.py rule -ex makedoc   # build FORD HTML docs → doc/html/
FoBiS.py rule -ex deldoc    # delete generated docs
```

### Coverage analysis

```bash
FoBiS.py rule -ex makecoverage          # build + run + gcov
FoBiS.py rule -ex makecoverage-analysis # also save markdown reports
```

### Building HDF5 from source

If you need a custom HDF5 build (different compiler version, different options), use the provided helper:

```bash
./scripts/hdf5_build.sh -get -build -sdk GNU
# Options: -sdk GNU | NVF
#          -hdf5 <path>   override HDF5 source path
#          -lsrc <path>   override library source path
```
