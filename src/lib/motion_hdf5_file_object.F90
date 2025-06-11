!< MOTIOn, HDF5 file object class.
module motion_hdf5_file_object
!< MOTIOn, HDF5 file object class.

use penf
use stringifor
use hdf5
use mpi

implicit none
private
public :: hdf5_file_object
public :: HDF5_PARAMETERS

type :: hdf5_parameters_object
   !< Global named constants (paramters) class (container) of HDF5 syntax.
   character(6) :: HDF5_DATASPACE_TYPE_SIMPLE = 'simple'
endtype hdf5_parameters_object
type(hdf5_parameters_object), parameter :: HDF5_PARAMETERS=hdf5_parameters_object() !< List of HDF5 named constants.

type :: hdf5_file_object
   !< HDF5 file object class.
   type(string)   :: filename           !< File name.
   integer(HID_T) :: hdf5=0_HID_T       !< HDF5 file identifier (logical unit).
   integer(HID_T) :: dspace_id=0_HID_T  !< HDF5 dataspace identifier.
   integer(I4P)   :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)   :: myrank=0_I4P       !< MPI ID process.
   integer(I4P)   :: error=0_I4P        !< IO Error status.
   contains
      ! file methods
      procedure, pass(self) :: close_file !< Close HDF5 file.
      procedure, pass(self) :: open_file  !< Open HDF5 file.
      ! data methods
      procedure, pass(self) :: close_dspace        !< Close HDF5 dataspace.
      procedure, pass(self) :: open_dspace         !< Open HDF5 dataspace.
      generic               :: save_dataset =>      &
                               save_dataset_4D_R8P, &
                               save_dataset_3D_R8P, &
                               save_dataset_0D_R8P !< Save dataset in dataspace.
      ! private methods
      procedure, pass(self), private :: save_dataset_4D_R8P !< Save dataset in dataspace, kind R8P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_R8P !< Save dataset in dataspace, kind R8P, rank 3D.
      procedure, pass(self), private :: save_dataset_0D_R8P !< Save dataset in dataspace, kind R8P, rank 1D.
endtype hdf5_file_object

contains
   ! files methods
   subroutine close_file(self)
   !< Close HDF5 file.
   class(hdf5_file_object), intent(inout) :: self !< File handler.

   ! close the file
   call h5fclose_f(self%hdf5, self%error)
   ! close fortran interface
   call h5close_f(self%error)
   endsubroutine close_file

   subroutine open_file(self, filename)
   !< Open HDF5 file.
   class(hdf5_file_object), intent(inout) :: self               !< File handler.
   character(*),            intent(in)    :: filename           !< File name.
   logical                                :: is_mpi_initialized !< MPI env status.

   ! reset file handler
   select type(self)
   type is(hdf5_file_object)
      self = hdf5_file_object()
   endselect
   call MPI_INITIALIZED(is_mpi_initialized, self%error)
   if (.not.is_mpi_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%filename = trim(adjustl(filename))
   ! open fortran interface
   call h5open_f(self%error)
   ! create a new file using default properties
   call h5fcreate_f(self%filename%chars(), H5F_ACC_TRUNC_F, self%hdf5, self%error)
   endsubroutine open_file

   ! data methods
   subroutine close_dspace(self)
   !< Close HDF5 dataspace.
   class(hdf5_file_object), intent(inout) :: self !< File handler.

   ! terminate access to the data space
   call h5sclose_f(self%dspace_id, self%error)
   endsubroutine close_dspace

   subroutine open_dspace(self, dataspace_type, nd)
   !< Open HDF5 dataspace.
   class(hdf5_file_object), intent(inout)        :: self           !< File handler.
   character(*),            intent(in)           :: dataspace_type !< Dataspace type.
   integer(HSIZE_T),        intent(in), optional :: nd(:)          !< Dataspace datasets dimensions.

   ! create the dataspace for fields
   select case(trim(adjustl(dataspace_type)))
   case(HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE)
      if (present(nd)) call h5screate_simple_f(size(nd), nd, self%dspace_id, self%error)
   endselect
   endsubroutine open_dspace

   subroutine save_dataset_4D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:,:) !< Dataset to be saved.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_NATIVE_DOUBLE, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_R8P

   subroutine save_dataset_3D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:) !< Dataset to be saved.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_NATIVE_DOUBLE, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_R8P

   subroutine save_dataset_0D_R8P(self, dset_name, dset)
   !< Save dataset in dataspace, kind R8P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   real(R8P),               intent(in)    :: dset      !< Dataset to be saved.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_NATIVE_DOUBLE, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_R8P
endmodule motion_hdf5_file_object
