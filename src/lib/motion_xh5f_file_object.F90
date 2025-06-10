!< MOTIOn, XDMF/HDF5 file object class.
module motion_xh5f_file_object
!< MOTIOn, XDMF/HDF5 file object class.

use motion_hdf5_file_object
use motion_xdmf_file_object
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf
use stringifor
use mpi

implicit none
private
public :: xh5f_file_object
public :: XH5F_BLOCK_CARTESIAN
public :: XH5F_BLOCK_CARTESIAN_UNIFORM
public :: XH5F_BLOCK_CURVILINEAR

character(9),  parameter :: XH5F_BLOCK_CARTESIAN         = 'cartesian'
character(17), parameter :: XH5F_BLOCK_CARTESIAN_UNIFORM = 'cartesian_uniform'
character(11), parameter :: XH5F_BLOCK_CURVILINEAR       = 'curvilinear'

type :: xh5f_file_object
   !< XDMF/HDF5 file object class.
   type(hdf5_file_object) :: hdf5               !< HDF5 file handler.
   type(xdmf_file_object) :: xdmf               !< XDMF file handler.
   integer(I4P)           :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)           :: myrank=0_I4P       !< MPI ID process.
   integer(I4P)           :: error=0_I4P        !< Error status.
   contains
      ! file methods
      procedure, pass(self) :: close_file !< Close XH5F file.
      procedure, pass(self) :: open_file  !< Open XH5F file.
      ! data methods
      procedure, pass(self) :: close_block             !< Close block.
      procedure, pass(self) :: close_grid              !< Close grid.
      procedure, pass(self) :: open_grid               !< Open grid.
      procedure, pass(self) :: open_block              !< Open block.
      generic               :: save_block_field =>      &
                               save_block_field_3D_R8P, &
                               save_block_field_0D_R8P !< Save field in block.
      ! private methods
      procedure, pass(self), private :: save_block_field_3D_R8P !< Save field in block, kind R8P, rank 3D.
      procedure, pass(self), private :: save_block_field_0D_R8P !< Save field in block, kind R8P, rank 0D.
endtype xh5f_file_object

contains
   ! files methods
   subroutine close_file(self)
   !< Close XH5F file.
   class(xh5f_file_object), intent(inout) :: self !< File handler.

   call self%hdf5%close_file
   call self%xdmf%close_domain_tag
   call self%xdmf%close_file
   endsubroutine close_file

   subroutine open_file(self, filename_hdf5, filename_xdmf)
   !< Open XH5F file.
   !< @NOTE MPI init must be invoked before this routine.
   class(xh5f_file_object), intent(inout) :: self               !< File handler.
   character(*),            intent(in)    :: filename_hdf5      !< File name of HDF5 file.
   character(*),            intent(in)    :: filename_xdmf      !< File name of XDMF file.
   character(:), allocatable              :: blocks_group_name_ !< Name of blocks group, local var.
   logical                                :: is_mpi_initialized !< MPI env status.

   ! reset file handler
   select type(self)
   type is(xh5f_file_object)
      self = xh5f_file_object()
   endselect
   call MPI_INITIALIZED(is_mpi_initialized, self%error)
   if (.not.is_mpi_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call self%hdf5%open_file(filename=filename_hdf5)
   call self%xdmf%open_file(filename=filename_xdmf)
   call self%xdmf%open_domain_tag
   endsubroutine open_file

   ! data methods
   subroutine close_block(self)
   !< Close block.
   class(xh5f_file_object), intent(inout) :: self !< File handler.

   call self%hdf5%close_dspace
   call self%xdmf%close_grid_tag
   endsubroutine close_block

   subroutine close_grid(self, grid_type)
   !< Close grid.
   class(xh5f_file_object), intent(inout)        :: self      !< File handler.
   character(*),            intent(in), optional :: grid_type !< Grid type.

   call self%xdmf%close_grid_tag(grid_type=grid_type)
   endsubroutine close_grid

   subroutine open_grid(self, grid_name, grid_type, grid_collection_type, grid_section)
   !< Open grid.
   class(xh5f_file_object), intent(inout)        :: self                 !< File handler.
   character(*),            intent(in), optional :: grid_name            !< Grid name.
   character(*),            intent(in), optional :: grid_type            !< Grid type.
   character(*),            intent(in), optional :: grid_collection_type !< Grid collection type.
   character(*),            intent(in), optional :: grid_section         !< Grid section.

   call self%xdmf%open_grid_tag(grid_name=grid_name, grid_type=grid_type, &
                                grid_collection_type=grid_collection_type, grid_section=grid_section)
   endsubroutine open_grid

   subroutine open_block(self, block_type, block_name, nijk, emin, dxyz, time)
   !< Open block.
   class(xh5f_file_object), intent(inout)        :: self        !< File handler.
   character(*),            intent(in)           :: block_type  !< Block type.
   character(*),            intent(in), optional :: block_name  !< Block name.
   integer(I8P),            intent(in), optional :: nijk(3)     !< Dataspace datasets dimensions.
   real(R8P),               intent(in), optional :: emin(3)     !< Block minimum extents.
   real(R8P),               intent(in), optional :: dxyz(3)     !< Block space steps.
   real(R8P),               intent(in), optional :: time        !< Current time.
   type(string)                                  :: block_name_ !< Block name, local var.

   block_name_ = block_name_%tempname(prefix='block-') ; if (present(block_name)) block_name_ = trim(adjustl(block_name))
   select case(trim(adjustl(block_type)))
   case(XH5F_BLOCK_CARTESIAN)
   case(XH5F_BLOCK_CARTESIAN_UNIFORM)
      if ((.not.present(nijk)).or.(.not.present(emin)).or.(.not.present(dxyz))) then
         write(stderr, '(A)') 'error: opening XH5F block of type "'//XH5F_BLOCK_CARTESIAN_UNIFORM//&
                              '" needs "nijk", "emin", "dxyz" dummy arguments'
         stop
      endif
      call self%hdf5%open_dspace(dataspace_type=HDF5_DATASPACE_TYPE_SIMPLE, nijk=nijk)
      call self%xdmf%open_grid_tag(grid_name=block_name_%chars())
      call self%xdmf%open_geometry_tag(geometry_type=XDMF_GEOMETRY_TYPE_ODXYZ)
      if (present(time)) call self%xdmf%write_time_tag(time_value=trim(str(time)))
      call self%xdmf%write_dataitem_tag(content          = trim(str([emin(3),emin(2),emin(1)],separator=' ')), &
                                        item_dimensions  = '3',                                                &
                                        number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                    &
                                        number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                   &
                                        number_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call self%xdmf%write_dataitem_tag(content          = trim(str([dxyz(3),dxyz(2),dxyz(1)],separator=' ')), &
                                        item_dimensions  = '3',                                                &
                                        number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                    &
                                        number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                   &
                                        number_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call self%xdmf%close_geometry_tag
      call self%xdmf%write_topology_tag(topology_type=XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                                        topology_dimensions=trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')))
   case(XH5F_BLOCK_CURVILINEAR)
      write(stderr, '(A)') 'error: block of type "'//trim(adjustl(block_type))//'" is not yet supported'
      stop
   case default
      write(stderr, '(A)') 'error: block of type "'//trim(adjustl(block_type))//'" is unknown'
      stop
   endselect
   endsubroutine open_block

   ! private methods
   subroutine save_block_field_3D_R8P(self, xdmf_field_name, nijk, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R8P, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(I8P),            intent(in)           :: nijk(:)          !< Dataspace datasets dimensions.
   real(R8P),               intent(in)           :: field(:,:,:)     !< Field to be saved.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: field_center_    !< Field center, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_center_   =XDMF_ATTR_CENTER_CELL          ; if (present(field_center   )) field_center_   =trim(adjustl(field_center   ))
   field_format_   =XDMF_DATAITEM_NUMBER_FORMAT_HDF; if (present(field_format   )) field_format_   =trim(adjustl(field_format   ))
   hdf5_field_name_=trim(adjustl(xdmf_field_name)) ; if (present(hdf5_field_name)) hdf5_field_name_=trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, nijk=nijk, dset=field)
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(reshape(field,[nijk(1)*nijk(2)*nijk(3)]),separator=' '))
   case(XDMF_DATAITEM_NUMBER_FORMAT_BINARY)
      write(stderr, '(A)') 'error: field format "'//trim(adjustl(field_format_))//'" is not yet supported'
      stop
   case default
      write(stderr, '(A)') 'error: field format "'//trim(adjustl(field_format_))//'" is unknown'
      stop
   endselect
   call self%xdmf%open_attribute_tag(attribute_name   = trim(adjustl(xdmf_field_name)), &
                                     attribute_center = field_center_,                  &
                                     attribute_type   = XDMF_ATTR_TYPE_SCALAR)
   call self%xdmf%write_dataitem_tag(content          = dataitem_content,                                   &
                                     item_dimensions  = trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')), &
                                     number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                    &
                                     number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                   &
                                     number_format    = field_format_)
   call self%xdmf%close_attribute_tag
   endsubroutine save_block_field_3D_R8P

   subroutine save_block_field_0D_R8P(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind R8P, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   real(R8P),               intent(in)           :: field            !< Field to be saved.
   character(*),            intent(in), optional :: field_center     !< Field center (Cell, Node, Grid...).
   character(*),            intent(in), optional :: field_format     !< Field format, HDF, XML, Binary.
   character(*),            intent(in), optional :: hdf5_field_name  !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format_    !< Field format, local var.
   character(:), allocatable                     :: field_center_    !< Field center, local var.
   character(:), allocatable                     :: hdf5_field_name_ !< Field name in HDF5 file, local var.
   character(:), allocatable                     :: dataitem_content !< Field content.

   field_center_   =XDMF_ATTR_CENTER_CELL          ; if (present(field_center   )) field_center_   =trim(adjustl(field_center   ))
   field_format_   =XDMF_DATAITEM_NUMBER_FORMAT_HDF; if (present(field_format   )) field_format_   =trim(adjustl(field_format   ))
   hdf5_field_name_=trim(adjustl(xdmf_field_name)) ; if (present(hdf5_field_name)) hdf5_field_name_=trim(adjustl(hdf5_field_name))
   select case(field_format_)
   case(XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%save_dataset(dset_name=hdf5_field_name_, dset=field)
      dataitem_content = self%hdf5%filename//':'//hdf5_field_name_
   case(XDMF_DATAITEM_NUMBER_FORMAT_XML)
      dataitem_content = trim(str(field))
   case(XDMF_DATAITEM_NUMBER_FORMAT_BINARY)
      write(stderr, '(A)') 'error: field format "'//trim(adjustl(field_format_))//'" is not yet supported'
      stop
   case default
      write(stderr, '(A)') 'error: field format "'//trim(adjustl(field_format_))//'" is unknown'
      stop
   endselect
   call self%xdmf%open_attribute_tag(attribute_name   = trim(adjustl(xdmf_field_name)), &
                                     attribute_center = field_center_,                  &
                                     attribute_type   = XDMF_ATTR_TYPE_SCALAR)
   call self%xdmf%write_dataitem_tag(content          = dataitem_content,                &
                                     item_dimensions  = '1',                             &
                                     number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT, &
                                     number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,&
                                     number_format    = field_format_)
   call self%xdmf%close_attribute_tag
   endsubroutine save_block_field_0D_R8P
endmodule motion_xh5f_file_object
