---
title: Usage
---

# Usage

The examples below are adapted from `src/tests/motion_write_xdmf_file_test.F90`.

## Quick start — XH5F (high-level API)

Import the library and declare a handler:

```fortran
use motion
implicit none
type(xh5f_file_object) :: xh5f
```

### Write a Cartesian uniform domain

The following writes one HDF5 file per MPI rank plus a single XDMF descriptor assembled by the master rank. All MPI coordination is internal — no explicit `if (myrank == 0)` branches needed.

```fortran
integer(HSIZE_T) :: nijk(3)
real(R8P)        :: emin(3), dxyz(3), time
integer(I4P)     :: b, myrank

nijk  = [16_HSIZE_T, 16_HSIZE_T, 16_HSIZE_T]
dxyz  = [1._R8P, 1._R8P, 1._R8P]
emin  = [0._R8P, 0._R8P, 0._R8P]
time  = 0.56_R8P

call xh5f%open_file(filename_hdf5='sim-mpi_'//trim(strz(myrank,2))//'.h5', &
                    filename_xdmf='sim-mpi_procs_02.xdmf',                 &
                    act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)

call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%open_grid(grid_name='mpi_'//trim(strz(myrank,2)), &
                    grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)

do b = 1, nb_proc  ! blocks owned by this rank
   call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM, &
                        block_name = 'block_'//trim(strz(global_b,2)),              &
                        nijk       = nijk,                                          &
                        emin       = emin,                                          &
                        dxyz       = dxyz,                                          &
                        time       = time)

   ! Grid-centered scalar stored inline (small value → XML)
   call xh5f%save_block_field(xdmf_field_name = 'Time',                                &
                              field           = time,                                  &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID, &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)

   ! Cell-centered scalar 3D stored in HDF5
   call xh5f%save_block_field(xdmf_field_name = 'pressure',                              &
                              nd              = nijk,                                    &
                              field           = pressure(1:nijk(1),1:nijk(2),1:nijk(3)), &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,   &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                              hdf5_field_name = 'block_01-pressure')

   call xh5f%close_block
enddo

call xh5f%close_grid
call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xh5f%close_file
```

### Vector, tensor6, tensor, and matrix fields

The leading dimension of `nd` selects the XDMF attribute type automatically:

```fortran
! Vector (3 components)
call xh5f%save_block_field(xdmf_field_name = 'velocity',                                       &
                           nd              = [3_HSIZE_T, nijk(1), nijk(2), nijk(3)],           &
                           field           = vel(1:3, 1:nijk(1), 1:nijk(2), 1:nijk(3)),        &
                           field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,            &
                           field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,  &
                           hdf5_field_name = 'block_01-velocity')

! Symmetric tensor (6 components)
call xh5f%save_block_field(xdmf_field_name = 'stress6',                                        &
                           nd              = [6_HSIZE_T, nijk(1), nijk(2), nijk(3)],           &
                           field           = sig6(1:6, 1:nijk(1), 1:nijk(2), 1:nijk(3)),       &
                           field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,            &
                           field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,  &
                           hdf5_field_name = 'block_01-stress6')

! Full tensor (9 components)
call xh5f%save_block_field(xdmf_field_name = 'stress',                                         &
                           nd              = [9_HSIZE_T, nijk(1), nijk(2), nijk(3)],           &
                           field           = sig(1:9, 1:nijk(1), 1:nijk(2), 1:nijk(3)),        &
                           field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,            &
                           field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,  &
                           hdf5_field_name = 'block_01-stress')

! General matrix (m components)
call xh5f%save_block_field(xdmf_field_name = 'matrix',                                         &
                           nd              = [12_HSIZE_T, nijk(1), nijk(2), nijk(3)],          &
                           field           = mat(1:12, 1:nijk(1), 1:nijk(2), 1:nijk(3)),       &
                           field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,            &
                           field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,  &
                           hdf5_field_name = 'block_01-matrix')
```

### Other Fortran kinds

The same method works for `R4P`, `I8P`, `I4P`, `I2P`, `I1P`:

```fortran
real(R4P) :: pressure_r4p(nijk(1), nijk(2), nijk(3))

call xh5f%save_block_field(xdmf_field_name = 'pressure-R4P',                             &
                           nd              = nijk,                                       &
                           field           = pressure_r4p,                               &
                           field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,     &
                           field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                           hdf5_field_name = 'block_01-pressure-R4P')
```

### Cartesian and curvilinear blocks

```fortran
! Cartesian (variable spacing) — pass per-axis node coordinate arrays
call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN, &
                     block_name = 'block_01',                            &
                     nijk       = nijk,                                  &
                     x          = x(0:nijk(1)),                          &
                     y          = y(0:nijk(2)),                          &
                     z          = z(0:nijk(3)),                          &
                     time       = time)

! Curvilinear — pass full node coordinate array [3, 0:ni, 0:nj, 0:nk]
call xh5f%open_block(block_type = XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR, &
                     block_name = 'block_01',                              &
                     nijk       = nijk,                                    &
                     nodes      = nodes(1:3, 0:nijk(1), 0:nijk(2), 0:nijk(3)), &
                     time       = time)
```

---

## Low-level API — HDF5 + XDMF separately

When you need full control over the XDMF structure, use `hdf5_file_object` and `xdmf_file_object` independently:

```fortran
type(hdf5_file_object) :: hdf5
type(xdmf_file_object) :: xdmf

call hdf5%open_file(filename='out.h5', act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
call xdmf%open_file(filename='out.xdmf', act=FILE_PARAMETERS%FILE_ACTION_OVERWRITE)

call xdmf%open_domain_tag
call xdmf%open_grid_tag(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xdmf%open_grid_tag(grid_name='mpi_00', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)

do b = 1, nb_proc
   call hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nijk)

   call xdmf%open_grid_tag(grid_name='block_'//trim(strz(b,2)))
   call xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_ODXYZ)
   call xdmf%write_dataitem_tag(content          = trim(str([emin(3),emin(2),emin(1)], separator=' ')), &
                                item_dimensions  = '3',                                                &
                                number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
                                number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,   &
                                number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
   call xdmf%write_dataitem_tag(content          = trim(str([dxyz(3),dxyz(2),dxyz(1)], separator=' ')), &
                                item_dimensions  = '3',                                                  &
                                number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,      &
                                number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,     &
                                number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
   call xdmf%close_geometry_tag
   call xdmf%write_topology_tag(topology_type       = XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                                topology_dimensions = trim(str([nijk(3),nijk(2),nijk(1)], separator=' ')))
   call xdmf%write_time_tag(time_value=trim(str(time)))

   ! Save data to HDF5 and reference it from XDMF
   call hdf5%save_dataset(dset_name='block_01-pressure', nd=nijk, dset=pressure)
   call xdmf%open_attribute_tag(attribute_name   = 'pressure',                              &
                                attribute_center = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,  &
                                attribute_type   = XDMF_PARAMETERS%XDMF_ATTR_TYPE_SCALAR)
   call xdmf%write_dataitem_tag(content          = hdf5%filename//':block_01-pressure',       &
                                item_dimensions  = trim(str([nijk(3),nijk(2),nijk(1)], separator=' ')), &
                                number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,     &
                                number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,    &
                                number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
   call xdmf%close_attribute_tag
   call hdf5%close_dspace
   call xdmf%close_grid_tag
enddo

call xdmf%close_grid_tag
call xdmf%close_grid_tag(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
call xdmf%close_domain_tag

call hdf5%close_file
call xdmf%close_file
```

---

## Reading back HDF5 data

```fortran
type(hdf5_file_object)        :: hdf5
real(R8P)                     :: time
integer(HSIZE_T), allocatable :: dims(:)

call hdf5%open_file(filename='out.h5', act=FILE_PARAMETERS%FILE_ACTION_READONLY)

! Load a scalar
call hdf5%load_dataset(dset_name='block_01-time', dset=time)

! Check existence
if (hdf5%does_dataset_exist(dset_name='block_01-pressure')) then
   ! Query dimensions before allocating
   call hdf5%get_dataset_dims(dset_name='block_01-pressure', nd=dims)
   ! ... allocate and load ...
endif

call hdf5%close_file
```
