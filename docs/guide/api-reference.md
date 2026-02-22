---
title: API Reference
---

# API Reference

The `motion` module re-exports all public types and parameter objects. Use a single `use` statement:

```fortran
use motion
```

---

## `file_base_object`

Abstract base type extended by all file handlers. Initializes MPI and exposes:

| Member | Type | Description |
|--------|------|-------------|
| `mpi_rank` | `integer(I4P)` | Rank of the calling process |
| `mpi_procs` | `integer(I4P)` | Total number of MPI processes |

Constants object: `FILE_PARAMETERS`

| Constant | Value | Description |
|----------|-------|-------------|
| `FILE_ACTION_OVERWRITE` | `'overwrite'` | Open file for writing, truncate if exists |
| `FILE_ACTION_READONLY` | `'readonly'` | Open file for reading |

---

## `xh5f_file_object`

High-level combined HDF5 + XDMF handler. Extends `file_base_object`.

Constants object: `XH5F_PARAMETERS`

| Constant | Value |
|----------|-------|
| `XH5F_BLOCK_CARTESIAN_UNIFORM` | `'cartesian_uniform'` |
| `XH5F_BLOCK_CARTESIAN` | `'cartesian'` |
| `XH5F_BLOCK_CURVILINEAR` | `'curvilinear'` |

### File lifecycle

```fortran
call xh5f%open_file(filename_hdf5='out.h5', filename_xdmf='out.xdmf', &
                    act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
call xh5f%close_file
```

### Grid / block lifecycle

```fortran
call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%open_grid(grid_name='mpi_00', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)

! Cartesian uniform block
call xh5f%open_block(block_type=XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                     block_name='block_01', nijk=nijk, emin=emin, dxyz=dxyz, time=t)

! Cartesian block (variable spacing)
call xh5f%open_block(block_type=XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN, &
                     block_name='block_01', nijk=nijk, x=x, y=y, z=z, time=t)

! Curvilinear block
call xh5f%open_block(block_type=XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR, &
                     block_name='block_01', nijk=nijk, nodes=nodes, time=t)

call xh5f%close_block
call xh5f%close_grid
call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
```

### `save_block_field`

Overloaded for scalar 0D (`field` is a scalar), 3-D arrays, and 4-D arrays (vector/tensor). Each variant accepts all PENF kinds.

```fortran
! Grid-centered scalar (0D), stored inline in XML
call xh5f%save_block_field(xdmf_field_name='Time', field=t, &
                           field_center=XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID, &
                           field_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)

! Cell-centered scalar 3D, stored in HDF5
call xh5f%save_block_field(xdmf_field_name='pressure', nd=nijk, &
                           field=pressure(1:ni,1:nj,1:nk), &
                           field_center=XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL, &
                           field_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                           hdf5_field_name='block_01-pressure')

! Cell-centered vector 3D
call xh5f%save_block_field(xdmf_field_name='velocity', nd=[3_HSIZE_T,ni,nj,nk], &
                           field=vel(1:3,1:ni,1:nj,1:nk), &
                           field_center=XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL, &
                           field_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                           hdf5_field_name='block_01-velocity')
```

### `load_block_field`

Mirror of `save_block_field` for reading. Same overloads.

---

## `hdf5_file_object`

Low-level HDF5 handler. Extends `file_base_object`.

Constants object: `HDF5_PARAMETERS`

| Constant | Value |
|----------|-------|
| `HDF5_DATASPACE_TYPE_SIMPLE` | `'simple'` |

### File and dataspace

```fortran
call hdf5%open_file(filename='out.h5', act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
call hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nijk)
! ... save/load datasets ...
call hdf5%close_dspace
call hdf5%close_file
```

### Dataset I/O

```fortran
! Save (overloaded: ranks 0-7, kinds R8P R4P I8P I4P I2P I1P)
call hdf5%save_dataset(dset_name='pressure', nd=nijk, dset=pressure(1:ni,1:nj,1:nk))

! Load
call hdf5%load_dataset(dset_name='pressure', dset=pressure)

! Scalar (0D)
call hdf5%save_dataset(dset_name='time', dset=t)
call hdf5%load_dataset(dset_name='time', dset=t)
```

### Introspection

```fortran
logical :: exists
integer(HSIZE_T), allocatable :: dims(:)

exists = hdf5%does_dataset_exist(dset_name='pressure')
call hdf5%get_dataset_dims(dset_name='pressure', nd=dims)
```

---

## `xdmf_file_object`

Low-level XDMF XML generator. Extends `file_base_object`.

Constants object: `XDMF_PARAMETERS` — contains all XDMF string constants (grid types, topology types, attribute centers, attribute types, geometry types, dataitem formats/types/precisions, time types).

### File lifecycle

```fortran
call xdmf%open_file(filename='out.xdmf', act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
call xdmf%open_domain_tag
! ... write tags ...
call xdmf%close_domain_tag
call xdmf%close_file
```

### Grid / attribute tags

```fortran
call xdmf%open_grid_tag(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xdmf%open_grid_tag(grid_name='block_01')
call xdmf%write_topology_tag(topology_type=XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                              topology_dimensions='16 16 16')
call xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_ODXYZ)
call xdmf%write_dataitem_tag(content='0.0 0.0 0.0', item_dimensions='3', &
                              number_type=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT, &
                              number_precision=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                              number_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
call xdmf%close_geometry_tag
call xdmf%write_time_tag(time_value='0.56')
call xdmf%open_attribute_tag(attribute_name='pressure', &
                              attribute_center=XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL, &
                              attribute_type=XDMF_PARAMETERS%XDMF_ATTR_TYPE_SCALAR)
call xdmf%write_dataitem_tag(content='out.h5:/block_01-pressure', &
                              item_dimensions='16 16 16', &
                              number_type=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT, &
                              number_precision=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                              number_format=XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
call xdmf%close_attribute_tag
call xdmf%close_grid_tag
call xdmf%close_grid_tag(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
```

### Reading an XDMF file

```fortran
call xdmf%open_file(filename='out.xdmf', act=FILE_PARAMETERS%FILE_ACTION_READONLY)
call xdmf%parse_file
call xdmf%close_file
```
