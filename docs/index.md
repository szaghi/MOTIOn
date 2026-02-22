---
layout: home

hero:
  name: MOTIOn
  text: Modular HPC IO for Fortran
  tagline: A pure Fortran 2003+ library providing a simple, agnostic API for parallel HDF5/XDMF IO in MPI-based HPC applications.
  actions:
    - theme: brand
      text: Guide
      link: /guide/
    - theme: alt
      text: API Reference
      link: /api/
    - theme: alt
      text: View on GitHub
      link: https://github.com/szaghi/MOTIOn

features:
  - icon: 🚀
    title: High-level XH5F API
    details: Single unified object that coordinates HDF5 data files and XDMF metadata in one call. MPI master/worker branching is handled transparently — no explicit rank checks in your code.
  - icon: 🗺️
    title: Three grid topologies
    details: Cartesian uniform (origin + spacing), Cartesian (per-axis node arrays), and Curvilinear (full 3-D node coordinate arrays) out of the box.
  - icon: 🔢
    title: All Fortran kinds
    details: save_block_field and load_block_field are overloaded for R8P, R4P, I8P, I4P, I2P, and I1P — scalars, 3-D fields, vectors, tensors, tensor6, and matrix shapes.
  - icon: 🔧
    title: Low-level HDF5 and XDMF APIs
    details: hdf5_file_object and xdmf_file_object give fine-grained control over dataspaces, attributes, geometries, and topologies when the high-level API is not enough.
  - icon: 🧪
    title: OOP / TDD Designed
    details: Three clean layers — file_base_object, format-specific objects (hdf5_file_object, xdmf_file_object), and the high-level xh5f_file_object — each tested individually.
  - icon: 🆓
    title: Free & Open Source
    details: Multi-licensed — GPLv3 for FOSS projects, BSD 2/3-Clause or MIT for commercial use. Fortran 2003+ standard compliant.
---

## Quick start

Let us consider a very simple numerical domain discretized by a collection of cartesian, uniform grids with the fields of pressure,
density and temperature integrated at cell center. The generation of output files in HDF5/XDMF form can be conveniently done by means
of the motion class `xh5f_file_object` by something like the following:

```fortran
call xh5f%open_file(filename_hdf5='simple-mpi_'//trim(strz(myrank,2))//'.h5', &
                    filename_xdmf='simple-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf')
call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%open_grid(grid_name='mpi_'//trim(strz(myrank,2)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
do b=1, domain%nb_proc
   call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                        block_name = 'block_'//trim(strz(mynb(1)-1+b,2)),          & ! global block numeration
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

In the above (incomplete) example (see the full source
[here](https://github.com/szaghi/MOTIOn/blob/main/src/tests/motion_write_xdmf_file_test.F90)) a set of HDF5 files are
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

## Authors

- Stefano Zaghi — [@szaghi](https://github.com/szaghi)

Contributions are welcome — see the [Contributing](/guide/contributing) page.

## Copyrights

FOSSIL is distributed under a multi-licensing system:

| Use case | License |
|----------|---------|
| FOSS projects | [GPL v3](http://www.gnu.org/licenses/gpl-3.0.html) |
| Closed source / commercial | [BSD 2-Clause](http://opensource.org/licenses/BSD-2-Clause) |
| Closed source / commercial | [BSD 3-Clause](http://opensource.org/licenses/BSD-3-Clause) |
| Closed source / commercial | [MIT](http://opensource.org/licenses/MIT) |

> Anyone interested in using, developing, or contributing to FOSSIL is welcome — pick the license that best fits your needs.
