!< MOTIOn, HDF5 file object class.
module motion_hdf5_file_object
!< MOTIOn, HDF5 file object class.

use motion_file_abst_object
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
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
   character(6) :: HDF5_DATASPACE_TYPE_SIMPLE = 'simple' !< Create simple type dataspace.
endtype hdf5_parameters_object
type(hdf5_parameters_object), parameter :: HDF5_PARAMETERS=hdf5_parameters_object() !< List of HDF5 named constants.

! HDF5 fortran interface parameters
integer(HID_T) :: H5T_R8P=0_HID_T            !< HDF5 kind mapping R8P kind.
integer(HID_T) :: H5T_R4P=0_HID_T            !< HDF5 kind mapping R4P kind.
integer(HID_T) :: H5T_I8P=0_HID_T            !< HDF5 kind mapping I8P kind.
integer(HID_T) :: H5T_I4P=0_HID_T            !< HDF5 kind mapping I4P kind.
integer(HID_T) :: H5T_I2P=0_HID_T            !< HDF5 kind mapping I2P kind.
integer(HID_T) :: H5T_I1P=0_HID_T            !< HDF5 kind mapping I1P kind.
logical        :: IS_H5F_INITIALIZED=.false. !< Status of HDF5 fortran interface inizialization.

type, extends(file_base_object) :: hdf5_file_object
   !< HDF5 file object class.
   integer(HID_T) :: hdf5=0_HID_T      !< HDF5 file identifier (logical unit).
   integer(HID_T) :: dspace_id=0_HID_T !< HDF5 dataspace identifier.
   contains
      ! public methods
      procedure, pass(self) :: h5f_initialize !< Initialize HDF5 fortran interface.
      procedure, pass(self) :: h5f_finalize   !< Finalize   HDF5 fortran interface.
      ! file methods
      procedure, pass(self) :: close_file !< Close HDF5 file.
      procedure, pass(self) :: open_file  !< Open HDF5 file.
      ! data methods
      procedure, pass(self) :: close_dspace !< Close HDF5 dataspace.
      procedure, pass(self) :: open_dspace  !< Open HDF5 dataspace.
      ! inquire dataset
      procedure, pass(self) :: does_dataset_exist !< Return true if dataset exist in file.
      procedure, pass(self) :: get_dataset_dims   !< Get dataset dimensions.
      ! load/save dataset
      generic               :: load_dataset =>      &
                               load_dataset_7D_R8P, &
                               load_dataset_6D_R8P, &
                               load_dataset_5D_R8P, &
                               load_dataset_4D_R8P, &
                               load_dataset_3D_R8P, &
                               load_dataset_2D_R8P, &
                               load_dataset_1D_R8P, &
                               load_dataset_0D_R8P, &
                               load_dataset_7D_R4P, &
                               load_dataset_6D_R4P, &
                               load_dataset_5D_R4P, &
                               load_dataset_4D_R4P, &
                               load_dataset_3D_R4P, &
                               load_dataset_2D_R4P, &
                               load_dataset_1D_R4P, &
                               load_dataset_0D_R4P, &
                               load_dataset_7D_I8P, &
                               load_dataset_6D_I8P, &
                               load_dataset_5D_I8P, &
                               load_dataset_4D_I8P, &
                               load_dataset_3D_I8P, &
                               load_dataset_2D_I8P, &
                               load_dataset_1D_I8P, &
                               load_dataset_0D_I8P, &
                               load_dataset_7D_I4P, &
                               load_dataset_6D_I4P, &
                               load_dataset_5D_I4P, &
                               load_dataset_4D_I4P, &
                               load_dataset_3D_I4P, &
                               load_dataset_2D_I4P, &
                               load_dataset_1D_I4P, &
                               load_dataset_0D_I4P, &
                               load_dataset_7D_I2P, &
                               load_dataset_6D_I2P, &
                               load_dataset_5D_I2P, &
                               load_dataset_4D_I2P, &
                               load_dataset_3D_I2P, &
                               load_dataset_2D_I2P, &
                               load_dataset_1D_I2P, &
                               load_dataset_0D_I2P, &
                               load_dataset_7D_I1P, &
                               load_dataset_6D_I1P, &
                               load_dataset_5D_I1P, &
                               load_dataset_4D_I1P, &
                               load_dataset_3D_I1P, &
                               load_dataset_2D_I1P, &
                               load_dataset_1D_I1P, &
                               load_dataset_0D_I1P !< Load dataset in dataspace.
      generic               :: save_dataset =>      &
                               save_dataset_7D_R8P, &
                               save_dataset_6D_R8P, &
                               save_dataset_5D_R8P, &
                               save_dataset_4D_R8P, &
                               save_dataset_3D_R8P, &
                               save_dataset_2D_R8P, &
                               save_dataset_1D_R8P, &
                               save_dataset_0D_R8P, &
                               save_dataset_7D_R4P, &
                               save_dataset_6D_R4P, &
                               save_dataset_5D_R4P, &
                               save_dataset_4D_R4P, &
                               save_dataset_3D_R4P, &
                               save_dataset_2D_R4P, &
                               save_dataset_1D_R4P, &
                               save_dataset_0D_R4P, &
                               save_dataset_7D_I8P, &
                               save_dataset_6D_I8P, &
                               save_dataset_5D_I8P, &
                               save_dataset_4D_I8P, &
                               save_dataset_3D_I8P, &
                               save_dataset_2D_I8P, &
                               save_dataset_1D_I8P, &
                               save_dataset_0D_I8P, &
                               save_dataset_7D_I4P, &
                               save_dataset_6D_I4P, &
                               save_dataset_5D_I4P, &
                               save_dataset_4D_I4P, &
                               save_dataset_3D_I4P, &
                               save_dataset_2D_I4P, &
                               save_dataset_1D_I4P, &
                               save_dataset_0D_I4P, &
                               save_dataset_7D_I2P, &
                               save_dataset_6D_I2P, &
                               save_dataset_5D_I2P, &
                               save_dataset_4D_I2P, &
                               save_dataset_3D_I2P, &
                               save_dataset_2D_I2P, &
                               save_dataset_1D_I2P, &
                               save_dataset_0D_I2P, &
                               save_dataset_7D_I1P, &
                               save_dataset_6D_I1P, &
                               save_dataset_5D_I1P, &
                               save_dataset_4D_I1P, &
                               save_dataset_3D_I1P, &
                               save_dataset_2D_I1P, &
                               save_dataset_1D_I1P, &
                               save_dataset_0D_I1P !< Save dataset in dataspace.
      ! private methods
      procedure, pass(self), private :: load_dataset_7D_R8P !< Load dataset in dataspace, kind R8P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_R8P !< Load dataset in dataspace, kind R8P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_R8P !< Load dataset in dataspace, kind R8P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_R8P !< Load dataset in dataspace, kind R8P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_R8P !< Load dataset in dataspace, kind R8P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_R8P !< Load dataset in dataspace, kind R8P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_R8P !< Load dataset in dataspace, kind R8P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_R8P !< Load dataset in dataspace, kind R8P, rank 0D.
      procedure, pass(self), private :: load_dataset_7D_R4P !< Load dataset in dataspace, kind R4P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_R4P !< Load dataset in dataspace, kind R4P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_R4P !< Load dataset in dataspace, kind R4P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_R4P !< Load dataset in dataspace, kind R4P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_R4P !< Load dataset in dataspace, kind R4P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_R4P !< Load dataset in dataspace, kind R4P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_R4P !< Load dataset in dataspace, kind R4P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_R4P !< Load dataset in dataspace, kind R4P, rank 0D.
      procedure, pass(self), private :: load_dataset_7D_I8P !< Load dataset in dataspace, kind I8P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_I8P !< Load dataset in dataspace, kind I8P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_I8P !< Load dataset in dataspace, kind I8P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_I8P !< Load dataset in dataspace, kind I8P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_I8P !< Load dataset in dataspace, kind I8P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_I8P !< Load dataset in dataspace, kind I8P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_I8P !< Load dataset in dataspace, kind I8P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_I8P !< Load dataset in dataspace, kind I8P, rank 0D.
      procedure, pass(self), private :: load_dataset_7D_I4P !< Load dataset in dataspace, kind I4P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_I4P !< Load dataset in dataspace, kind I4P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_I4P !< Load dataset in dataspace, kind I4P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_I4P !< Load dataset in dataspace, kind I4P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_I4P !< Load dataset in dataspace, kind I4P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_I4P !< Load dataset in dataspace, kind I4P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_I4P !< Load dataset in dataspace, kind I4P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_I4P !< Load dataset in dataspace, kind I4P, rank 0D.
      procedure, pass(self), private :: load_dataset_7D_I2P !< Load dataset in dataspace, kind I2P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_I2P !< Load dataset in dataspace, kind I2P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_I2P !< Load dataset in dataspace, kind I2P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_I2P !< Load dataset in dataspace, kind I2P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_I2P !< Load dataset in dataspace, kind I2P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_I2P !< Load dataset in dataspace, kind I2P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_I2P !< Load dataset in dataspace, kind I2P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_I2P !< Load dataset in dataspace, kind I2P, rank 0D.
      procedure, pass(self), private :: load_dataset_7D_I1P !< Load dataset in dataspace, kind I1P, rank 7D.
      procedure, pass(self), private :: load_dataset_6D_I1P !< Load dataset in dataspace, kind I1P, rank 6D.
      procedure, pass(self), private :: load_dataset_5D_I1P !< Load dataset in dataspace, kind I1P, rank 5D.
      procedure, pass(self), private :: load_dataset_4D_I1P !< Load dataset in dataspace, kind I1P, rank 4D.
      procedure, pass(self), private :: load_dataset_3D_I1P !< Load dataset in dataspace, kind I1P, rank 3D.
      procedure, pass(self), private :: load_dataset_2D_I1P !< Load dataset in dataspace, kind I1P, rank 2D.
      procedure, pass(self), private :: load_dataset_1D_I1P !< Load dataset in dataspace, kind I1P, rank 1D.
      procedure, pass(self), private :: load_dataset_0D_I1P !< Load dataset in dataspace, kind I1P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_R8P !< Save dataset in dataspace, kind R8P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_R8P !< Save dataset in dataspace, kind R8P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_R8P !< Save dataset in dataspace, kind R8P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_R8P !< Save dataset in dataspace, kind R8P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_R8P !< Save dataset in dataspace, kind R8P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_R8P !< Save dataset in dataspace, kind R8P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_R8P !< Save dataset in dataspace, kind R8P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_R8P !< Save dataset in dataspace, kind R8P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_R4P !< Save dataset in dataspace, kind R4P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_R4P !< Save dataset in dataspace, kind R4P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_R4P !< Save dataset in dataspace, kind R4P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_R4P !< Save dataset in dataspace, kind R4P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_R4P !< Save dataset in dataspace, kind R4P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_R4P !< Save dataset in dataspace, kind R4P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_R4P !< Save dataset in dataspace, kind R4P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_R4P !< Save dataset in dataspace, kind R4P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_I8P !< Save dataset in dataspace, kind I8P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_I8P !< Save dataset in dataspace, kind I8P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_I8P !< Save dataset in dataspace, kind I8P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_I8P !< Save dataset in dataspace, kind I8P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_I8P !< Save dataset in dataspace, kind I8P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_I8P !< Save dataset in dataspace, kind I8P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_I8P !< Save dataset in dataspace, kind I8P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_I8P !< Save dataset in dataspace, kind I8P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_I4P !< Save dataset in dataspace, kind I4P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_I4P !< Save dataset in dataspace, kind I4P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_I4P !< Save dataset in dataspace, kind I4P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_I4P !< Save dataset in dataspace, kind I4P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_I4P !< Save dataset in dataspace, kind I4P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_I4P !< Save dataset in dataspace, kind I4P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_I4P !< Save dataset in dataspace, kind I4P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_I4P !< Save dataset in dataspace, kind I4P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_I2P !< Save dataset in dataspace, kind I2P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_I2P !< Save dataset in dataspace, kind I2P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_I2P !< Save dataset in dataspace, kind I2P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_I2P !< Save dataset in dataspace, kind I2P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_I2P !< Save dataset in dataspace, kind I2P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_I2P !< Save dataset in dataspace, kind I2P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_I2P !< Save dataset in dataspace, kind I2P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_I2P !< Save dataset in dataspace, kind I2P, rank 0D.
      procedure, pass(self), private :: save_dataset_7D_I1P !< Save dataset in dataspace, kind I1P, rank 7D.
      procedure, pass(self), private :: save_dataset_6D_I1P !< Save dataset in dataspace, kind I1P, rank 6D.
      procedure, pass(self), private :: save_dataset_5D_I1P !< Save dataset in dataspace, kind I1P, rank 5D.
      procedure, pass(self), private :: save_dataset_4D_I1P !< Save dataset in dataspace, kind I1P, rank 4D.
      procedure, pass(self), private :: save_dataset_3D_I1P !< Save dataset in dataspace, kind I1P, rank 3D.
      procedure, pass(self), private :: save_dataset_2D_I1P !< Save dataset in dataspace, kind I1P, rank 2D.
      procedure, pass(self), private :: save_dataset_1D_I1P !< Save dataset in dataspace, kind I1P, rank 1D.
      procedure, pass(self), private :: save_dataset_0D_I1P !< Save dataset in dataspace, kind I1P, rank 0D.
endtype hdf5_file_object

contains
   ! public methods
   subroutine h5f_initialize(self, verbose, prefix)
   !< Initialize HDF5 fortran interface.
   class(hdf5_file_object), intent(inout)        :: self    !< File handler.
   logical,                 intent(in), optional :: verbose !< Verbose output.
   character(*),            intent(in), optional :: prefix  !< Output prefix.
   character(:), allocatable                     :: prefix_ !< Output prefix, local var.

   if (.not.IS_H5F_INITIALIZED) then
      ! open fortran interface
      call h5open_f(self%error)
      H5T_R8P = h5kind_to_type(R8P,H5_REAL_KIND)
      H5T_R4P = h5kind_to_type(R4P,H5_REAL_KIND)
      H5T_I8P = h5kind_to_type(I8P,H5_INTEGER_KIND)
      H5T_I4P = h5kind_to_type(I4P,H5_INTEGER_KIND)
      H5T_I2P = h5kind_to_type(I2P,H5_INTEGER_KIND)
      H5T_I1P = h5kind_to_type(I1P,H5_INTEGER_KIND)
      IS_H5F_INITIALIZED=.true.
      if (present(verbose)) then
         if (verbose) then
            prefix_ = '' ; if (present(prefix)) prefix_ = trim(adjustl(prefix))
            print '(A)', prefix_//'H5T_R8P kind mapping: "'//trim(str(H5T_R8P))//'"'
            print '(A)', prefix_//'H5T_R4P kind mapping: "'//trim(str(H5T_R4P))//'"'
            print '(A)', prefix_//'H5T_I8P kind mapping: "'//trim(str(H5T_I8P))//'"'
            print '(A)', prefix_//'H5T_I4P kind mapping: "'//trim(str(H5T_I4P))//'"'
            print '(A)', prefix_//'H5T_I2P kind mapping: "'//trim(str(H5T_I2P))//'"'
            print '(A)', prefix_//'H5T_I1P kind mapping: "'//trim(str(H5T_I1P))//'"'
         endif
      endif
   endif
   endsubroutine h5f_initialize

   subroutine h5f_finalize(self)
   !< Finalize HDF5 fortran interface.
   class(hdf5_file_object), intent(inout) :: self !< File handler.

   if (.not.IS_H5F_INITIALIZED) then
      ! close fortran interface
      call h5close_f(self%error)
      IS_H5F_INITIALIZED=.false.
   endif
   endsubroutine h5f_finalize

   ! files methods
   subroutine close_file(self)
   !< Close HDF5 file.
   class(hdf5_file_object), intent(inout) :: self !< File handler.

   ! close the file
   call h5fclose_f(self%hdf5, self%error)
   endsubroutine close_file

   subroutine open_file(self, filename, act)
   !< Open HDF5 file.
   class(hdf5_file_object), intent(inout)        :: self     !< File handler.
   character(*),            intent(in)           :: filename !< File name.
   character(*),            intent(in), optional :: act      !< File action ['readonly, overwrite'...].
   character(:), allocatable                     :: act_     !< File action, local var.

   act_ = FILE_PARAMETERS%FILE_ACTION_READONLY ; if (present(act)) act_ = trim(adjustl(act))
   ! reset file handler
   select type(self)
   type is(hdf5_file_object)
      self = hdf5_file_object()
   endselect
   call self%file_base_object%initialize
   self%filename = trim(adjustl(filename))
   ! open fortran interface (if necessary)
   call self%h5f_initialize
   ! open file
   select case(act_)
   case(FILE_PARAMETERS%FILE_ACTION_READONLY)
      ! open file in readonly, exit fail if it does not exist
      call h5fopen_f(self%filename%chars(), H5F_ACC_RDONLY_F, self%hdf5, self%error)
   case(FILE_PARAMETERS%FILE_ACTION_OVERWRITE)
      ! create new file, overwrite if already exist
      call h5fcreate_f(self%filename%chars(), H5F_ACC_TRUNC_F, self%hdf5, self%error)
   case(FILE_PARAMETERS%FILE_ACTION_NEWFILE)
      ! create new file, exit fail if already exist
      call h5fcreate_f(self%filename%chars(), H5F_ACC_EXCL_F, self%hdf5, self%error)
   case default
      write(stderr, '(A)') 'error: file action "'//act_//'" is unknown'
      call MPI_FINALIZE(self%error)
      stop
   endselect
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

   ! inquire dataset
   function does_dataset_exist(self, dset_name) result(does_exist)
   !< Return true if dataset exist in file.
   class(hdf5_file_object), intent(inout) :: self       !< File handler.
   character(*),            intent(in)    :: dset_name  !< Dataset name.
   logical                                :: does_exist !< Inquire result.

   call h5lexists_f(self%hdf5, trim(adjustl(dset_name)), does_exist, self%error)
   endfunction does_dataset_exist

   subroutine get_dataset_dims(self, dset_name, nd, nd_max)
   !< Get dataset dimensions.
   !< @NOTE Currently only simple dataspace is supported.
   class(hdf5_file_object),       intent(inout)         :: self       !< File handler.
   character(*),                  intent(in)            :: dset_name  !< Dataset name.
   integer(HSIZE_T), allocatable, intent(out)           :: nd(:)      !< Dataspace datasets dimensions.
   integer(HSIZE_T), allocatable, intent(out), optional :: nd_max(:)  !< Dataspace datasets maximum dimensions.
   integer(HID_T)                                       :: dset_id    !< Dataset identifier.
   integer(I4P)                                         :: rank       !< Dataset rank.
   integer(HSIZE_T), allocatable                        :: nd_max_(:) !< Dataspace datasets maximum dimensions, local var.

   if (self%does_dataset_exist(dset_name)) then
      call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
      call h5dget_space_f(dset_id, self%dspace_id, self%error)
      call h5sget_simple_extent_ndims_f(self%dspace_id, rank, self%error)
      allocate(nd(     rank))
      allocate(nd_max_(rank))
      call h5sget_simple_extent_dims_f(self%dspace_id, nd, nd_max_, self%error)
      call h5dclose_f(dset_id, self%error)
      call h5sclose_f(self%dspace_id, self%error)
      if (present(nd_max)) nd_max = nd_max_
   endif
   endsubroutine get_dataset_dims

   ! load dataset
   ! R8P
   subroutine load_dataset_7D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_R8P

   subroutine load_dataset_6D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_R8P

   subroutine load_dataset_5D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_R8P

   subroutine load_dataset_4D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_R8P

   subroutine load_dataset_3D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_R8P

   subroutine load_dataset_2D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_R8P

   subroutine load_dataset_1D_R8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R8P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R8P),               intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_R8P

   subroutine load_dataset_0D_R8P(self, dset_name, dset)
   !< Load dataset in dataspace, kind R8P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   real(R8P),               intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   real(R8P)                              :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R8P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_R8P

   ! R4P
   subroutine load_dataset_7D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_R4P

   subroutine load_dataset_6D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_R4P

   subroutine load_dataset_5D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_R4P

   subroutine load_dataset_4D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_R4P

   subroutine load_dataset_3D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_R4P

   subroutine load_dataset_2D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_R4P

   subroutine load_dataset_1D_R4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind R4P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R4P),               intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_R4P

   subroutine load_dataset_0D_R4P(self, dset_name, dset)
   !< Load dataset in dataspace, kind R4P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   real(R4P),               intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   real(R4P)                              :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_R4P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_R4P

   ! I8P
   subroutine load_dataset_7D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_I8P

   subroutine load_dataset_6D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_I8P

   subroutine load_dataset_5D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_I8P

   subroutine load_dataset_4D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_I8P

   subroutine load_dataset_3D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_I8P

   subroutine load_dataset_2D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_I8P

   subroutine load_dataset_1D_I8P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I8P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I8P),            intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_I8P

   subroutine load_dataset_0D_I8P(self, dset_name, dset)
   !< Load dataset in dataspace, kind I8P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I8P),            intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   integer(I8P)                           :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I8P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_I8P

   ! I4P
   subroutine load_dataset_7D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_I4P

   subroutine load_dataset_6D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_I4P

   subroutine load_dataset_5D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_I4P

   subroutine load_dataset_4D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_I4P

   subroutine load_dataset_3D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_I4P

   subroutine load_dataset_2D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_I4P

   subroutine load_dataset_1D_I4P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I4P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I4P),            intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_I4P

   subroutine load_dataset_0D_I4P(self, dset_name, dset)
   !< Load dataset in dataspace, kind I4P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I4P),            intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   integer(I4P)                           :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I4P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_I4P

   ! I2P
   subroutine load_dataset_7D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_I2P

   subroutine load_dataset_6D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_I2P

   subroutine load_dataset_5D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_I2P

   subroutine load_dataset_4D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_I2P

   subroutine load_dataset_3D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_I2P

   subroutine load_dataset_2D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_I2P

   subroutine load_dataset_1D_I2P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I2P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I2P),            intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_I2P

   subroutine load_dataset_0D_I2P(self, dset_name, dset)
   !< Load dataset in dataspace, kind I2P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I2P),            intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   integer(I2P)                           :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I2P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_I2P

   ! I1P
   subroutine load_dataset_7D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_7D_I1P

   subroutine load_dataset_6D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_6D_I1P

   subroutine load_dataset_5D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_5D_I1P

   subroutine load_dataset_4D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_4D_I1P

   subroutine load_dataset_3D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_3D_I1P

   subroutine load_dataset_2D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_2D_I1P

   subroutine load_dataset_1D_I1P(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind I1P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I1P),            intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_1D_I1P

   subroutine load_dataset_0D_I1P(self, dset_name, dset)
   !< Load dataset in dataspace, kind I1P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I1P),            intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   integer(I1P)                           :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_I1P, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine load_dataset_0D_I1P

   ! save dataset
   ! R8P
   subroutine save_dataset_7D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_R8P

   subroutine save_dataset_6D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_R8P

   subroutine save_dataset_5D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_R8P

   subroutine save_dataset_4D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_R8P

   subroutine save_dataset_3D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_R8P

   subroutine save_dataset_2D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_R8P

   subroutine save_dataset_1D_R8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R8P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   real(R8P),               intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_R8P

   subroutine save_dataset_0D_R8P(self, dset_name, dset)
   !< Save dataset in dataspace, kind R8P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   real(R8P),               intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R8P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_R8P

   ! R4P
   subroutine save_dataset_7D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_R4P

   subroutine save_dataset_6D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_R4P

   subroutine save_dataset_5D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_R4P

   subroutine save_dataset_4D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_R4P

   subroutine save_dataset_3D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_R4P

   subroutine save_dataset_2D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_R4P

   subroutine save_dataset_1D_R4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind R4P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   real(R4P),               intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_R4P

   subroutine save_dataset_0D_R4P(self, dset_name, dset)
   !< Save dataset in dataspace, kind R4P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   real(R4P),               intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_R4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_R4P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_R4P

   ! I8P
   subroutine save_dataset_7D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_I8P

   subroutine save_dataset_6D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_I8P

   subroutine save_dataset_5D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_I8P

   subroutine save_dataset_4D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_I8P

   subroutine save_dataset_3D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_I8P

   subroutine save_dataset_2D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_I8P

   subroutine save_dataset_1D_I8P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I8P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   integer(I8P),            intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_I8P

   subroutine save_dataset_0D_I8P(self, dset_name, dset)
   !< Save dataset in dataspace, kind I8P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I8P),            intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I8P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I8P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_I8P

   ! I4P
   subroutine save_dataset_7D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_I4P

   subroutine save_dataset_6D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_I4P

   subroutine save_dataset_5D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_I4P

   subroutine save_dataset_4D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_I4P

   subroutine save_dataset_3D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_I4P

   subroutine save_dataset_2D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_I4P

   subroutine save_dataset_1D_I4P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I4P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   integer(I4P),            intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_I4P

   subroutine save_dataset_0D_I4P(self, dset_name, dset)
   !< Save dataset in dataspace, kind I4P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I4P),            intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I4P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I4P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_I4P

   ! I2P
   subroutine save_dataset_7D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_I2P

   subroutine save_dataset_6D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_I2P

   subroutine save_dataset_5D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_I2P

   subroutine save_dataset_4D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_I2P

   subroutine save_dataset_3D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_I2P

   subroutine save_dataset_2D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_I2P

   subroutine save_dataset_1D_I2P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I2P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   integer(I2P),            intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_I2P

   subroutine save_dataset_0D_I2P(self, dset_name, dset)
   !< Save dataset in dataspace, kind I2P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I2P),            intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I2P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I2P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_I2P

   ! I1P
   subroutine save_dataset_7D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_7D_I1P

   subroutine save_dataset_6D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_6D_I1P

   subroutine save_dataset_5D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_5D_I1P

   subroutine save_dataset_4D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_4D_I1P

   subroutine save_dataset_3D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_3D_I1P

   subroutine save_dataset_2D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_2D_I1P

   subroutine save_dataset_1D_I1P(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind I1P, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   integer(I1P),            intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_1D_I1P

   subroutine save_dataset_0D_I1P(self, dset_name, dset)
   !< Save dataset in dataspace, kind I1P, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(I1P),            intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_I1P, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_I1P, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine save_dataset_0D_I1P

endmodule motion_hdf5_file_object
