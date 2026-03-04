# MOTIOn

>#### Modular HPC Optimized Toolkit for IO in Fortran
>a pure Fortran 2003+ library providing a simple, agnostic API for parallel HDF5/XDMF I/O in MPI-based HPC applications.

[![CI](https://github.com/szaghi/MOTIOn/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/MOTIOn/actions)
[![Coverage](https://img.shields.io/codecov/c/github/szaghi/MOTIOn.svg)](https://app.codecov.io/gh/szaghi/MOTIOn)
[![GitHub tag](https://img.shields.io/github/tag/szaghi/MOTIOn.svg)](https://github.com/szaghi/MOTIOn/releases)
[![License](https://img.shields.io/badge/license-GPLv3%20%7C%20BSD%20%7C%20MIT-blue.svg)](#copyrights)

| 💾 **HDF5+XDMF output**<br>One call writes binary HDF5 data and XML metadata; MPI master/worker branching handled transparently | 🧱 **Three grid topologies**<br>Cartesian uniform, Cartesian, and Curvilinear — declare a block type and MOTIOn handles the geometry | 🔄 **Kind-generic field I/O**<br>`save_block_field`/`load_block_field` overloaded for all Fortran kinds (R8P…I1P) and shapes (3-D field, vector, tensor, matrix) | ⚙️ **Low-level control**<br>Access `hdf5_file_object` and `xdmf_file_object` directly for custom dataspaces, attributes, and XML topologies |
|:---:|:---:|:---:|:---:|
| 🏗️ **Three-layer OOP**<br>`file_base_object` → format objects → `xh5f_file_object`; extend any layer for custom formats | 🖥️ **MPI-parallel**<br>Designed for distributed CFD solvers; each rank writes its own HDF5 file while the XDMF master collects metadata | 📦 **Pre-built HDF5**<br>GNU 14.2+ and nvfortran 25.3 HDF5 1.14.6 libraries included — no separate HDF5 build required | 🔤 **Pure Fortran 2003+**<br>Standard-compliant, no external dependencies beyond MPI and HDF5 |

>#### [Documentation](https://szaghi.github.io/MOTIOn/)
> For full documentation (guide, API reference, examples, etc...) see the [MOTIOn website](https://szaghi.github.io/MOTIOn/).

---

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](https://szaghi.github.io/MOTIOn/guide/contributing) page.

## Copyrights

This project is distributed under a multi-licensing system:

- **FOSS projects**: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html)
- **Closed source / commercial**: [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause), [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause), or [MIT](http://opensource.org/licenses/MIT)

> Anyone interested in using, developing, or contributing to this project is welcome — pick the license that best fits your needs.

---

## Quick start

Write a parallel HDF5/XDMF output from a collection of Cartesian uniform blocks:

```fortran
use motion
implicit none
type(xh5f_file_object) :: xh5f

call xh5f%open_file(filename_hdf5='out_rank'//trim(strz(myrank,2))//'.h5', &
                    filename_xdmf='out_procs'//trim(strz(nprocs,2))//'.xdmf')
call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
do b = 1, nblocks
  call xh5f%open_block(block_type=XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                       block_name='block_'//trim(strz(b,2)),                    &
                       nijk=nijk, emin=emin(:,b), dxyz=dxyz, time=time)
  call xh5f%save_block_field(xdmf_field_name='pressure',                              &
                             nijk=nijk,                                               &
                             field=pressure(1:nijk(1),1:nijk(2),1:nijk(3),b),         &
                             field_center=XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,      &
                             field_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
  call xh5f%close_block
enddo
call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%close_file
```

See the full example in [`src/tests/motion_write_xdmf_file_test.F90`](src/tests/motion_write_xdmf_file_test.F90).

---

## Install

### FoBiS

**Standalone** — clone and fetch dependencies, then build:

```bash
git clone https://github.com/szaghi/MOTIOn && cd MOTIOn
FoBiS.py fetch
FoBiS.py build -mode tests-gnu && ./scripts/run_tests.sh --np 2
```

**As a project dependency** — declare MOTIOn in your `fobos` and run `fetch`:

```ini
[dependencies]
deps_dir = src/third_party
MOTIOn = https://github.com/szaghi/MOTIOn
```

```bash
FoBiS.py fetch           # fetch and build
FoBiS.py fetch --update  # re-fetch and rebuild
```

| Mode | Command |
|------|---------|
| Shared library | `FoBiS.py build -mode shared-gnu` |
| Static library | `FoBiS.py build -mode static-gnu` |
| Tests (debug)  | `FoBiS.py build -mode tests-gnu-debug` |
| NVIDIA Fortran | `FoBiS.py build -mode tests-nvf` |
