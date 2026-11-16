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

## One-line install

```bash
curl -fsSL https://github.com/szaghi/MOTIOn/releases/latest/download/install.sh | bash
```

It installs the latest release into `~/.local` with the build system you choose among those MOTIOn ships
(`--build fobis|fpm`: a menu asks when you do not say), one install per compiler (`lib/motion/<compiler>-<version>/`,
`include/motion/<compiler>-<version>/`), and makes it usable through pkg-config (`motion`, with
`PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig`). The release carries the Fortran dependencies (PENF, StringiFor, FoXy,
BeFoR64, FACE); MPI and a parallel HDF5 with the Fortran API are yours: the fobis build compiles against the HDF5 under
`lib/hdf5/`, the fpm build against the one pkg-config finds (`hdf5_fortran`). Build your program with the MPI wrapper
and both packages:

```bash
mpifort my_program.F90 $(pkg-config --cflags --libs motion hdf5_fortran) -o my_program
```

Pass options after `bash -s --`, e.g. `... | bash -s -- --build fpm --prefix /opt/motion`; `--help` lists them all,
`--uninstall` removes the install. The installer is FoBiS's
([One-line installer](https://szaghi.github.io/FoBiS/cookbook/installer)).

## Download

Clone the repository, then fetch the third-party dependencies via FoBiS.py:

```bash
git clone https://github.com/szaghi/MOTIOn
cd MOTIOn
FoBiS.py fetch
```

`FoBiS.py fetch` reads `src/third_party/.deps_config.ini` and downloads the required libraries into `src/third_party/`.

### Third-party dependencies

The dependencies are placed under `src/third_party/`:

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

# Run all tests (serial)
./scripts/run_tests.sh

# Run all tests via MPI with N ranks
./scripts/run_tests.sh --np 2
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
