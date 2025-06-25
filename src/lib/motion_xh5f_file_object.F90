!< MOTIOn, XDMF/HDF5 file object class.
module motion_xh5f_file_object
!< MOTIOn, XDMF/HDF5 file object class.

use motion_file_abst_object
use motion_hdf5_file_object
use motion_xdmf_file_object
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf
use stringifor
use hdf5
use mpi

implicit none
private
public :: xh5f_file_object
public :: XH5F_PARAMETERS

type :: xh5f_parameters_object
   !< Global named constants (paramters) class (container) of XH5F syntax.
   character(9)  :: XH5F_BLOCK_CARTESIAN         = 'cartesian'
   character(17) :: XH5F_BLOCK_CARTESIAN_UNIFORM = 'cartesian_uniform'
   character(11) :: XH5F_BLOCK_CURVILINEAR       = 'curvilinear'
endtype xh5f_parameters_object
type(xh5f_parameters_object), parameter :: XH5F_PARAMETERS=xh5f_parameters_object() !< List of XH5F named constants.

type, extends(file_base_object) :: xh5f_file_object
   !< XDMF/HDF5 file object class.
   type(hdf5_file_object)    :: hdf5             !< HDF5 file handler.
   type(xdmf_file_object)    :: xdmf             !< XDMF file handler.
   character(:), allocatable :: act              !< Action on files (readonly, overwrite, new...).
   logical                   :: with_xdmf=.true. !< Sentinel to exclude XDMF file.
   contains
      ! file methods
      procedure, pass(self) :: close_file !< Close XH5F file.
      procedure, pass(self) :: open_file  !< Open XH5F file.
      ! data methods
      procedure, pass(self) :: close_block                !< Close block.
      procedure, pass(self) :: close_grid                 !< Close grid.
      procedure, pass(self) :: open_grid                  !< Open grid.
      procedure, pass(self) :: open_block                 !< Open block.
      procedure, pass(self) :: save_block_field_xdmf_tags !< Save field in block, only XDMF tags.
      generic               :: load_block_field =>      &
                               load_block_field_4D_R8P    !< Load field in block.
      generic               :: save_block_field =>      &
                               save_block_field_4D_R8P, &
                               save_block_field_3D_R8P, &
                               save_block_field_0D_R8P, &
                               save_block_field_4D_R4P, &
                               save_block_field_3D_R4P, &
                               save_block_field_0D_R4P, &
                               save_block_field_4D_I8P, &
                               save_block_field_3D_I8P, &
                               save_block_field_0D_I8P, &
                               save_block_field_4D_I4P, &
                               save_block_field_3D_I4P, &
                               save_block_field_0D_I4P, &
                               save_block_field_4D_I2P, &
                               save_block_field_3D_I2P, &
                               save_block_field_0D_I2P, &
                               save_block_field_4D_I1P, &
                               save_block_field_3D_I1P, &
                               save_block_field_0D_I1P    !< Save field in block.
      ! private methods
      procedure, pass(self), private :: load_block_field_4D_R8P !< Load field in block, kind R8P, rank 4D.
      procedure, pass(self), private :: save_block_field_4D_R8P !< Save field in block, kind R8P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_R8P !< Save field in block, kind R8P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_R8P !< Save field in block, kind R8P, rank 0D.
      procedure, pass(self), private :: save_block_field_4D_R4P !< Save field in block, kind R4P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_R4P !< Save field in block, kind R4P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_R4P !< Save field in block, kind R4P, rank 0D.
      procedure, pass(self), private :: save_block_field_4D_I8P !< Save field in block, kind I8P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_I8P !< Save field in block, kind I8P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_I8P !< Save field in block, kind I8P, rank 0D.
      procedure, pass(self), private :: save_block_field_4D_I4P !< Save field in block, kind I4P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_I4P !< Save field in block, kind I4P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_I4P !< Save field in block, kind I4P, rank 0D.
      procedure, pass(self), private :: save_block_field_4D_I2P !< Save field in block, kind I2P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_I2P !< Save field in block, kind I2P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_I2P !< Save field in block, kind I2P, rank 0D.
      procedure, pass(self), private :: save_block_field_4D_I1P !< Save field in block, kind I1P, rank 4D.
      procedure, pass(self), private :: save_block_field_3D_I1P !< Save field in block, kind I1P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_I1P !< Save field in block, kind I1P, rank 0D.
endtype xh5f_file_object

contains
   ! files methods
   subroutine close_file(self)
   !< Close XH5F file.
   class(xh5f_file_object), intent(inout) :: self !< File handler.

   call self%hdf5%close_file
   if (self%with_xdmf) then
      call self%xdmf%close_domain_tag
      call self%xdmf%close_file
   endif
   endsubroutine close_file

   subroutine open_file(self, filename_hdf5, filename_xdmf, act)
   !< Open XH5F file.
   !< @NOTE MPI init must be invoked before this routine.
   class(xh5f_file_object), intent(inout)        :: self               !< File handler.
   character(*),            intent(in)           :: filename_hdf5      !< File name of HDF5 file.
   character(*),            intent(in), optional :: filename_xdmf      !< File name of XDMF file.
   character(*),            intent(in), optional :: act                !< File action ['readonly, overwrite'...].
   character(:), allocatable                     :: blocks_group_name_ !< Name of blocks group, local var.

   ! reset file handler
   select type(self)
   type is(xh5f_file_object)
      self = xh5f_file_object()
   endselect
   call self%file_base_object%initialize
   self%act = FILE_PARAMETERS%FILE_ACTION_READONLY ; if (present(act)) self%act = trim(adjustl(act))
   call self%hdf5%open_file(filename=filename_hdf5, act=act)
   if (present(filename_xdmf)) then
      self%with_xdmf = .true.
      call self%xdmf%open_file(filename=filename_xdmf, act=act)
      call self%xdmf%open_domain_tag
   else
      self%with_xdmf = .false.
   endif
   endsubroutine open_file

   ! data methods
   subroutine close_block(self)
   !< Close block.
   class(xh5f_file_object), intent(inout) :: self !< File handler.

   if (self%with_xdmf) call self%xdmf%close_grid_tag
   endsubroutine close_block

   subroutine close_grid(self, grid_type)
   !< Close grid.
   class(xh5f_file_object), intent(inout)        :: self      !< File handler.
   character(*),            intent(in), optional :: grid_type !< Grid type.

   if (self%with_xdmf) call self%xdmf%close_grid_tag(grid_type=grid_type)
   endsubroutine close_grid

   subroutine open_grid(self, grid_name, grid_type, grid_collection_type, grid_section)
   !< Open grid.
   class(xh5f_file_object), intent(inout)        :: self                 !< File handler.
   character(*),            intent(in), optional :: grid_name            !< Grid name.
   character(*),            intent(in), optional :: grid_type            !< Grid type.
   character(*),            intent(in), optional :: grid_collection_type !< Grid collection type.
   character(*),            intent(in), optional :: grid_section         !< Grid section.

   if (self%with_xdmf) then
      call self%xdmf%open_grid_tag(grid_name=grid_name, grid_type=grid_type, &
                                   grid_collection_type=grid_collection_type, grid_section=grid_section)
   endif
   endsubroutine open_grid

   subroutine open_block(self, block_type, block_name, nijk, emin, dxyz, x, y, z, nodes, time)
   !< Open block.
   class(xh5f_file_object), intent(inout)        :: self           !< File handler.
   character(*),            intent(in)           :: block_type     !< Block type.
   character(*),            intent(in), optional :: block_name     !< Block name.
   integer(HSIZE_T),        intent(in), optional :: nijk(3)        !< Cells number.
   real(R8P),               intent(in), optional :: emin(3)        !< Block minimum extents.
   real(R8P),               intent(in), optional :: dxyz(3)        !< Space steps for cartesian unform grid.
   real(R8P),               intent(in), optional :: x(:),y(:),z(:) !< Nodes coordinates for cartesian grid [1:nijk(ijk)+1].
   real(R8P),               intent(in), optional :: nodes(:,:,:,:) !< Nodes coordinates for curvilinear grid [3,nijk+1].
   real(R8P),               intent(in), optional :: time           !< Current time.
   integer(HSIZE_T), allocatable                 :: nd(:)          !< Dataspece dimensions.
   type(string)                                  :: block_name_    !< Block name, local var.

   block_name_ = block_name_%tempname(prefix='block-') ; if (present(block_name)) block_name_ = trim(adjustl(block_name))
   select case(trim(adjustl(block_type)))
   case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN)
      if ((.not.present(nijk)).or.(.not.present(x)).or.(.not.present(y)).or.(.not.present(z))) then
         write(stderr, '(A)') 'error: opening XH5F block of type "'//trim(adjustl(block_type))//&
                              '" needs "nijk", "x", "y", "z" dummy arguments'
         call MPI_FINALIZE(self%error)
         stop
      endif
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[nijk(1)+1])
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-x', nd=nijk(1)+1, dset=x)
      call self%hdf5%close_dspace
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[nijk(2)+1])
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-y', nd=nijk(2)+1, dset=y)
      call self%hdf5%close_dspace
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[nijk(3)+1])
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-z', nd=nijk(3)+1, dset=z)
      call self%hdf5%close_dspace
      if (self%with_xdmf) then
         call self%xdmf%open_grid_tag(grid_name=block_name_%chars())
         call self%xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_VXVYVZ)
         if (present(time)) call self%xdmf%write_time_tag(time_value=trim(str(time)))
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-x', &
                                           item_dimensions  = trim(adjustl(str(nijk(1)+1))),                      &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,   &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-y', &
                                           item_dimensions  = trim(adjustl(str(nijk(2)+1))),                      &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,   &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-z', &
                                           item_dimensions  = trim(adjustl(str(nijk(3)+1))),                      &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,   &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%close_geometry_tag
         call self%xdmf%write_topology_tag(topology_type=XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DRECTMESH, &
                                           topology_dimensions=trim(str([nijk(3)+1,nijk(2)+1,nijk(1)+1],separator=' ')))
      endif
   case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM)
      if ((.not.present(nijk)).or.(.not.present(emin)).or.(.not.present(dxyz))) then
         write(stderr, '(A)') 'error: opening XH5F block of type "'//trim(adjustl(block_type))//&
                              '" needs "nijk", "emin", "dxyz" dummy arguments'
         call MPI_FINALIZE(self%error)
         stop
      endif
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[3_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-origin', nd=3_HSIZE_T, dset=[emin(3),emin(2),emin(1)])
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-dxdydz', nd=3_HSIZE_T, dset=[dxyz(3),dxyz(2),dxyz(1)])
      call self%hdf5%close_dspace
      if (self%with_xdmf) then
         call self%xdmf%open_grid_tag(grid_name=block_name_%chars())
         call self%xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_ODXYZ)
         if (present(time)) call self%xdmf%write_time_tag(time_value=trim(str(time)))
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-origin', &
                                           item_dimensions  = '3',                                                     &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,         &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,        &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-dxdydz', &
                                           item_dimensions  = '3',                                                     &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,         &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,        &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%close_geometry_tag
         call self%xdmf%write_topology_tag(topology_type=XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                                           topology_dimensions=trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')))
      endif
   case(XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR)
      if ((.not.present(nijk)).or.(.not.present(nodes))) then
         write(stderr, '(A)') 'error: opening XH5F block of type "'//trim(adjustl(block_type))//&
                              '" needs "nijk", "nodes" dummy arguments'
         call MPI_FINALIZE(self%error)
         stop
      endif
      nd = [3_HSIZE_T,nijk(1)+1,nijk(2)+1,nijk(3)+1]
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=block_name_%chars()//'-nodes', nd=nd, dset=nodes)
      call self%hdf5%close_dspace
      if (self%with_xdmf) then
         call self%xdmf%open_grid_tag(grid_name=block_name_%chars())
         call self%xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_XYZ)
         if (present(time)) call self%xdmf%write_time_tag(time_value=trim(str(time)))
         call self%xdmf%write_dataitem_tag(content          = self%hdf5%filename//':'//block_name_%chars()//'-nodes',      &
                                           item_dimensions  = trim(adjustl(str([nd(4),nd(3),nd(2),nd(1)],separator=' '))), &
                                           number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,             &
                                           number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,            &
                                           number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call self%xdmf%close_geometry_tag
         call self%xdmf%write_topology_tag(topology_type=XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DSMESH, &
                                           topology_dimensions=trim(str([nd(4),nd(3),nd(2)],separator=' ')))
      endif
   case default
      write(stderr, '(A)') 'error: block of type "'//trim(adjustl(block_type))//'" is unknown'
      call MPI_FINALIZE(self%error)
      stop
   endselect
   endsubroutine open_block

   ! private methods
   ! load block field
   ! R8P
   subroutine load_block_field_4D_R8P(self, nd, field, xdmf_field_name, hdf5_field_name)
   !< Load field in block, kind R8P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self            !< File handler.
   integer(HSIZE_T),        intent(in)           :: nd(:)           !< Dataspace datasets dimensions.
   real(R8P),               intent(inout)        :: field(:,:,:,:)  !< Field.
   character(*),            intent(in), optional :: xdmf_field_name !< Field name in XDMF file.
   character(*),            intent(in), optional :: hdf5_field_name !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format    !< Field format, local var.
   character(:), allocatable                     :: field_name      !< Field name.

   select case(field_format)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%load_dataset(dset_name=field_name, nd=nd, dset=field)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
   case default
      write(stderr, '(A)') 'error: load block field error, field_format "'//field_format//'" unknown'
      call MPI_FINALIZE(self%error)
      stop
   endselect
   endsubroutine load_block_field_4D_R8P

   ! save block field
   subroutine save_block_field_xdmf_tags(self,number_type,number_precision,dataitem_content, &
                                         xdmf_field_name,field_format,nd,field_center)
   !< Save field in block, only XDMF tags.
   class(xh5f_file_object),   intent(inout)        :: self             !< File handler.
   character(*),              intent(in)           :: number_type      !< Number type.
   character(*),              intent(in)           :: number_precision !< Number precision.
   character(*),              intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   character(*),              intent(in)           :: dataitem_content !< Field content.
   character(*),              intent(in)           :: field_format     !< Field format, HDF, XML, Binary.
   integer(HSIZE_T),          intent(in), optional :: nd(:)            !< Dataspace datasets dimensions.
   character(*),              intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(:), allocatable                       :: field_center_    !< Field center, local var.
   character(:), allocatable                       :: item_dimensions  !< Field dimensions.
   character(:), allocatable                       :: attribute_type   !< Field type (scalar, vector, tensor, matrix).
   integer(I4P)                                    :: i                !< Counter.

   if (self%with_xdmf) then
      field_center_ = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL ; if (present(field_center)) field_center_ = trim(adjustl(field_center))
      item_dimensions = '1' ; if (present(nd)) item_dimensions = trim(str([(nd(i),i=size(nd),1,-1)],separator=' '))
      attribute_type = XDMF_PARAMETERS%XDMF_ATTR_TYPE_SCALAR
      if (present(nd)) then
         if (size(nd)>3) then
            select case(nd(1))
            case(2:3)
               attribute_type = XDMF_PARAMETERS%XDMF_ATTR_TYPE_VECTOR
            case(6)
               attribute_type = XDMF_PARAMETERS%XDMF_ATTR_TYPE_TENSOR6
            case(4:5, 7:9)
               attribute_type = XDMF_PARAMETERS%XDMF_ATTR_TYPE_TENSOR
            case(10:)
               attribute_type = XDMF_PARAMETERS%XDMF_ATTR_TYPE_MATRIX
            endselect
         endif
      endif
      call self%xdmf%open_attribute_tag(attribute_name   = trim(adjustl(xdmf_field_name)), &
                                        attribute_center = field_center_,                  &
                                        attribute_type   = attribute_type)
      call self%xdmf%write_dataitem_tag(content          = trim(adjustl(dataitem_content)), &
                                        item_dimensions  = item_dimensions,                 &
                                        number_type      = trim(adjustl(number_type)),      &
                                        number_precision = trim(adjustl(number_precision)), &
                                        number_format    = trim(adjustl(field_format)))
      call self%xdmf%close_attribute_tag
   endif
   endsubroutine save_block_field_xdmf_tags

   ! R8P
   subroutine save_block_field_4D_R8P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R8P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   real(R8P),               intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_R8P

   subroutine save_block_field_3D_R8P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R8P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   real(R8P),               intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_R8P

   subroutine save_block_field_0D_R8P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R8P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   real(R8P),               intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_R8P

   ! R4P
   subroutine save_block_field_4D_R4P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R4P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   real(R4P),               intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_R4P

   subroutine save_block_field_3D_R4P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R4P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   real(R4P),               intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_R4P

   subroutine save_block_field_0D_R4P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R4P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   real(R4P),               intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,  &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_R4P

   ! I8P
   subroutine save_block_field_4D_I8P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I8P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_I8P

   subroutine save_block_field_3D_I8P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I8P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_I8P

   subroutine save_block_field_0D_I8P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I8P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(I8P),            intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_I8P

   ! I4P
   subroutine save_block_field_4D_I4P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I4P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_I4P

   subroutine save_block_field_3D_I4P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I4P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_I4P

   subroutine save_block_field_0D_I4P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I4P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(I4P),            intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_4, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_I4P

   ! I2P
   subroutine save_block_field_4D_I2P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I2P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_2, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_I2P

   subroutine save_block_field_3D_I2P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I2P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_2, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_I2P

   subroutine save_block_field_0D_I2P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I2P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(I2P),            intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_2, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_I2P

   ! I1P
   subroutine save_block_field_4D_I1P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I1P, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)           :: field(:,:,:,:)   !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_1, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_4D_I1P

   subroutine save_block_field_3D_I1P(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I1P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)           :: field(:,:,:)     !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nd)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nd=nd, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[product(nd)]),separator=' '))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_1, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        nd               = nd,                                               &
                                        field_center     = field_center)
   endsubroutine save_block_field_3D_I1P

   subroutine save_block_field_0D_I1P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind I1P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(I1P),            intent(in)           :: field            !< Field.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_format_ = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF
   if (present(field_format)) field_format_ = trim(adjustl(field_format))
   hdf5_field_name_ = trim(adjustl(xdmf_field_name))
   if (present(hdf5_field_name)) hdf5_field_name_ = trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=[1_HSIZE_T])
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      call self%hdf5%close_dspace
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   endselect
   call self%save_block_field_xdmf_tags(number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_INT,    &
                                        number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_1, &
                                        dataitem_content = dataitem_content,                                 &
                                        xdmf_field_name  = xdmf_field_name,                                  &
                                        field_format     = field_format,                                     &
                                        field_center     = field_center)
   endsubroutine save_block_field_0D_I1P
endmodule motion_xh5f_file_object
