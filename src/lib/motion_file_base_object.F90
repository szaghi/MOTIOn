!< MOTIOn, abstract file object class.
module motion_file_abst_object
!< MOTIOn, abstract file object class.

use penf
use stringifor
use mpi

implicit none
private
public :: file_base_object
public :: FILE_PARAMETERS

type :: file_parameters_object
   !< Global named constants (paramters) class (container) for files handling.
   character(8) :: FILE_ACTION_READONLY  = 'readonly'  !< Open file in readonly, exit fail if it does not exist.
   character(9) :: FILE_ACTION_OVERWRITE = 'overwrite' !< Open new file, overwrite is already exist.
   character(7) :: FILE_ACTION_NEWFILE   = 'newfile'   !< Open new file, exit fail if already exist.
endtype file_parameters_object
type(file_parameters_object), parameter :: FILE_PARAMETERS=file_parameters_object() !< List of file named constants.

type :: file_base_object
   !< Abstract file object class.
   type(string) :: filename=string()  !< File name.
   integer(I4P) :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P) :: myrank=0_I4P       !< MPI ID process.
   integer(I4P) :: error=0_I4P        !< IO Error status.
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize file class.
endtype file_base_object

contains
   ! public methods
   subroutine initialize(self)
   !< Initialize file class.
   class(file_base_object), intent(inout) :: self               !< File handler.
   logical                                :: is_mpi_initialized !< MPI env status.

   ! reset file handler
   select type(self)
   type is(file_base_object)
      self = file_base_object()
   endselect
   ! initialize MPI env
   call MPI_INITIALIZED(is_mpi_initialized, self%error)
   if (.not.is_mpi_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   endsubroutine initialize
endmodule motion_file_abst_object
