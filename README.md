<a name="top"></a>

# MOTIOn [![GitHub tag](https://img.shields.io/github/tag/szaghi/MOTIOn.svg)](https://github.com/szaghi/MOTIOn/releases)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3%20,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

<!-- [![Status](https://img.shields.io/badge/status-stable-brightgreen.svg)]() -->
<!-- [![CI Status](https://github.com/szaghi/MOTIOn/actions/workflows/ci.yml/badge.svg)](https://github.com/szaghi/MOTIOn/actions) -->
<!-- [![Coverage Status](https://img.shields.io/codecov/c/github/szaghi/MOTIOn.svg)](https://app.codecov.io/gh/szaghi/MOTIOn) -->

### MOTIOn, Modular (HPC) Optimized Toolkit (for) IO (in fortra)n
A KISS, modular library for handling IO in HPC scenario for modern Fortran projects.

+ MOTIOn is a pure Fortran (KISS) library for handling efficient IO in HPC scenario for modern Fortran projects;
+ MOTIOn is Fortran 2003+ standard compliant;
+ MOTIOn is a Free, Open Source Project.

#### Table of Contents

- [What is MOTIOn?](#what-is-motion?)
- [Main features](#main-features)
- [Copyrights](#copyrights)
- [Documentation](#documentation)
	- [A Taste of MOTIOn](#a-taste-of-motion)

#### Issues

[![GitHub issues](https://img.shields.io/github/issues/szaghi/MOTIOn.svg)]()

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-v14.2.0+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/INTEL%20XL-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/NVIDIA%20XL-not%20tested-yellow.svg)]()

## What is MOTIOn?

MOTIOn, aka Modular (HPC) Optimized Toolkit (for) IO (in fortran)n, is a modular library for handling Input/Output operations
in HPC scenario for modern Fortran projects. MOTIOn aims to provide a simple, agnostic API to easy handle IO files tailored
to large, parallel applications where MPI is usually used. In particular, the main goal is to let final user deal with only
a very reduced set of operations while MOTIOn deals, in background, with the complexity of generating efficient, parallel IO,
e.g. exploiting HPC library like HDF5 with XDMF XML descriptor. The main focus is on CFD solver where a simulation domain is
discretized (in any form) and where some governing laws are integrated (common PDE problems): in such a scenario a set of
numerical grids (the number of which, in HPC applications, can be huge) and their integrated fields are processed, in parallel,
by a finite number of MPI processes that have the necessity to perform IO as efficient as possible; MOTIOn provides a convenient
API to IO grids/fields exploiting, in background, cutting edge libraries like HDF5/XDMF with very low effort.

Go to [Top](#top)

## Main features

Currently MOTIOn has the following features:

+ [ ] HDF5/XDMF:
   + [ ] Output:
      + [x] Cartesian uniform grids;
      + [ ] Cartesian grids;
      + [ ] Curvilinear grids;
      + [x] Scalar 3D fields;
      + [x] Scalar 0D (grid centered) fields;
      + [ ] Tensor 3D fields;
      + [ ] Matrix 3D fields;
   + [ ] Input:

To be completed.

Go to [Top](#top)

## Copyrights

MOTIOn is an open source project, it is distributed under a multi-licensing system:

+ for FOSS projects:
  - [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html);
+ for closed source/commercial projects:
  - [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause);
  - [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause);
  - [MIT](http://opensource.org/licenses/MIT).

Anyone is interest to use, to develop or to contribute to MOTIOn is welcome, feel free to select the license that best matches your soul!

More details can be found on [wiki](https://github.com/szaghi/MOTIOn/wiki/Copyrights).

Go to [Top](#top)

## Documentation

Besides this README file the MOTIOn documentation is contained into its own [wiki](https://github.com/szaghi/MOTIOn/wiki).
Detailed documentation of the API is contained into the [GitHub Pages](http://szaghi.github.io/MOTIOn/index.html) that can also
be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

### A Taste of MOTIOn

Let us consider a very simple numerical domain discretized by a collection of cartesian, uniform grids with the fields of pressure,
density and temperature integrated at cell center. The generation of output files in HDF5/XDMF form can be conveniently done by means
of the motion class `xh5f_file_object` by something like the following:

```fortran
call xh5f%open_file(filename_hdf5='simple-mpi_'//trim(strz(myrank,2))//'.h5', &
                    filename_xdmf='simple-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf')
call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%open_grid(grid_name='mpi_'//trim(strz(myrank,2)), grid_type=XDMF_GRID_TYPE_COLLECTION)
do b=1, domain%nb_proc
   call xh5f%open_block(block_type = XH5F_BLOCK_CARTESIAN_UNIFORM,        &
                        block_name = 'block_'//trim(strz(mynb(1)-1+b,2)), & ! global block numeration
                        nijk       = nijk,                                &
                        emin       = domain%emin(:,b),                    &
                        dxyz       = domain%dxyz,                         &
                        time       = domain%time)
   call xh5f%save_block_field(xdmf_field_name = 'Time',                &
                              field           = domain%time,           &
                              field_center    = XDMF_ATTR_CENTER_GRID, &
                              field_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
   do v=1, domain%nv
      call xh5f%save_block_field(xdmf_field_name = field_name(v)%chars(),                   &
                                 nijk            = nijk,                                    &
                                 field           = field(v,1:nijk(1),1:nijk(2),1:nijk(3),b),&
                                 field_center    = XDMF_ATTR_CENTER_CELL,                   &
                                 field_format    = XDMF_DATAITEM_NUMBER_FORMAT_HDF,         &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//field_name(v)%chars())
   enddo
   call xh5f%close_block
enddo
call xh5f%close_grid
call xh5f%close_grid(grid_type=XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%close_file
```

In the above (incomplete) example (see the full source
[here](src/tests/motion_write_xdmf_file_test.F90)) a set of HDF5 files are
created, one for each MPI processes used, as well as one XDMF file describing
the data contained into all HDF5 files. Each MPI process generates its own HDF5
file (in parallel and asyncronously) containing the data of the grids/fields
assigned to it. Also, each MPI process generates its own XDMF file structure,
but only the master process (myrank=0) creates a real XDMF file into witch it
writes the description of its data and also the descritpion of the data of the
other processes: as a matter of facts, MOTIOn objects (files classes) are able
to automatically understand which is the master process that must gather the
informations from other processes, committing a MPI_GATHERV only when all
processes have created their own XDMF structure. The creation of HDF5 files and
XDMF structures happens asyncronously, and the gathering of data for the final
XDMF file happens only on the last instructions.

The above example also shows some MOTIOn features:
+ a set (quite large) of global, named constants (parameters) are provided, e.g. `XDMF_ATTR_CENTER_GRID`,
  `XDMF_DATAITEM_NUMBER_FORMAT_XML`, ecc...; these parameters simplify the handling of HDF5/XDMF syntax;
+ the MPI *branching* between master process and other processes is automatically performed in background without bothering the
  end user;
+ the `xh5f_file_object` class provides an high level API with a simple, reduced set of methods, e.g. to save a scalar 0D field and
  a 3D one the same `save_block_field` method is used, without the necessity to use specialized API.

To be completed.

Go to [Top](#top)
