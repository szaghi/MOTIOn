# MOTIOn

**Modular HPC Optimized Toolkit for IO in Fortran** — a pure Fortran 2003+ library providing a simple, agnostic API for parallel HDF5/XDMF I/O in MPI-based HPC applications.

[![GitHub tag](https://img.shields.io/github/tag/szaghi/MOTIOn.svg)](https://github.com/szaghi/MOTIOn/releases)
[![License](https://img.shields.io/badge/license-GPLv3%20%7C%20BSD%20%7C%20MIT-blue.svg)](#copyrights)
[![Compiler](https://img.shields.io/badge/GNU-v14.2+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/nvfortran-25.3-brightgreen.svg)]()

---

## Features

- High-level `xh5f_file_object` API — one call writes HDF5 data and XDMF metadata; MPI master/worker branching handled transparently
- Three grid topologies out of the box: Cartesian uniform, Cartesian, and Curvilinear
- `save_block_field` / `load_block_field` overloaded for all Fortran kinds (R8P…I1P) and shapes (scalar, 3-D field, vector, tensor, tensor6, matrix)
- Low-level `hdf5_file_object` and `xdmf_file_object` for fine-grained control over dataspaces, attributes, geometries, and topologies
- OOP three-layer architecture — `file_base_object` → format objects → `xh5f_file_object`
- Pure Fortran 2003+ standard compliant

**[Documentation](https://szaghi.github.io/MOTIOn/)** | **[API Reference](https://szaghi.github.io/MOTIOn/api/)**

---

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](https://szaghi.github.io/MOTIOn/guide/contributing) page.

## Copyrights

This project is distributed under a multi-licensing system:

- **FOSS projects**: [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html)
- **Closed source / commercial**: [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause), [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause), or [MIT](http://opensource.org/licenses/MIT)

> Anyone interested in using, developing, or contributing to MOTIOn is welcome — pick the license that best fits your needs.

---

## Quick start

Write a parallel HDF5/XDMF output from a collection of Cartesian uniform blocks:

```fortran
call xh5f%open_file(filename_hdf5='simple-mpi_'//trim(strz(myrank,2))//'.h5', &
                    filename_xdmf='simple-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf')
call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%open_grid(grid_name='mpi_'//trim(strz(myrank,2)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
do b=1, domain%nb_proc
   call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                        block_name = 'block_'//trim(strz(mynb(1)-1+b,2)),          &
                        nijk       = nijk,                                         &
                        emin       = domain%emin(:,b),                             &
                        dxyz       = domain%dxyz,                                  &
                        time       = domain%time)
   call xh5f%save_block_field(xdmf_field_name = 'Time',                                &
                              field           = domain%time,                           &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID, &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
   do v=1, domain%nv
      call xh5f%save_block_field(xdmf_field_name = field_name(v)%chars(),                           &
                                 nijk            = nijk,                                            &
                                 field           = field(v,1:nijk(1),1:nijk(2),1:nijk(3),b),        &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,           &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//field_name(v)%chars())
   enddo
   call xh5f%close_block
enddo
call xh5f%close_grid
call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%close_file
```

See the full example in [`src/tests/motion_write_xdmf_file_test.F90`](src/tests/motion_write_xdmf_file_test.F90).

---

## Install

```sh
git clone https://github.com/szaghi/MOTIOn --recursive
cd MOTIOn
FoBiS.py build -mode tests-gnu
mpirun -np 2 ./exe/motion_write_xdmf_file_test
```

| Mode | Command |
|------|---------|
| Shared library | `FoBiS.py build -mode shared-gnu` |
| Static library | `FoBiS.py build -mode static-gnu` |
| Tests (debug)  | `FoBiS.py build -mode tests-gnu-debug` |
| NVIDIA Fortran | `FoBiS.py build -mode tests-nvf` |

Pre-built HDF5 1.14.6 libraries for GNU and nvfortran are included in `lib/`.
