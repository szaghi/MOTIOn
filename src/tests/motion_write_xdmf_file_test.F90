!< MOTIOn test: write XDMF file test.

module motion_test_domain_object
!< MOTIOn test: prototype of domain to be saved.
use penf
use stringifor
use mpi

implicit none
private
public :: domain_object

type :: domain_object
   !< Prototype of domain to be saved.
   !< For the sake of simplicity all domain blocks have the same cells number and spacing.
   integer(I4P)              :: procs_number=1_I4P             !< Number of MPI processes.
   integer(I4P)              :: myrank=0_I4P                   !< MPI ID process.
   integer(I4P)              :: nb_proc=0_I4P                  !< Number of blocks for each MPI process.
   integer(I4P)              :: nb=0_I4P                       !< Total number of blocks (all processes).
   integer(I4P)              :: mynb(2)=[0_I4P,0_I4P]          !< Number of blocks of current process [start-b,end-b].
   integer(I4P)              :: nv=0_I4P                       !< Number of fields variables.
   integer(I4P)              :: gc=0_I4P                       !< Ghost cells number.
   integer(I8P)              :: nijk(3)=[0_I4P,0_I4P,0_I4P]    !< Blocks dimensions.
   real(R8P)                 :: dxyz(3)=[0._R8P,0._R8P,0._R8P] !< Blocks space steps.
   real(R8P)                 :: time=0._R8P                    !< Time of domain solution.
   real(R8P),    allocatable :: emin(:,:)                      !< Blocks minimum extents [3,nb].
   real(R8P),    allocatable :: field(:,:,:,:,:)               !< Blocks fields [1:nv,1-gc:ni+gc,1-gc:nj+gc,1-gc:nk+gc,1:nb].
   type(string), allocatable :: field_name(:)                  !< Fields names [1:nv].
   integer(I4P)              :: error                          !< Error status.
   contains
      procedure, pass(self) :: initialize !< Initialize domain.
endtype domain_object

contains
   subroutine initialize(self, nb, nv, gc, nijk, dxyz, field_name, time)
   !< Initialize domain.
   class(domain_object), intent(inout) :: self           !< Domain.
   integer(I4P),         intent(in)    :: nb             !< Number of blocks.
   integer(I4P),         intent(in)    :: nv             !< Number of fields variables.
   integer(I4P),         intent(in)    :: gc             !< Ghost cells number.
   integer(I8P),         intent(in)    :: nijk(3)        !< Blocks dimensions.
   real(R8P),            intent(in)    :: dxyz(3)        !< Blocks space steps.
   character(*),         intent(in)    :: field_name(nv) !< Fields names.
   real(R8P),            intent(in)    :: time           !< Time of domain solution.
   integer(I4P)                        :: i,j,k,b,v      !< Counter.

   call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%nb_proc = nb / self%procs_number
   self%mynb = [self%nb_proc * self%myrank + 1, self%nb_proc * (self%myrank + 1)]
   self%nb   = nb
   self%nv   = nv
   self%gc   = gc
   self%nijk = nijk
   self%dxyz = dxyz
   self%time = time
   if (allocated(self%emin      )) deallocate(self%emin      ) ; allocate(self%emin(1:3,1:self%nb_proc))
   if (allocated(self%field     )) deallocate(self%field     ) ; allocate(self%field(1:nv,            &
                                                                                     1-gc:nijk(1)+gc, &
                                                                                     1-gc:nijk(2)+gc, &
                                                                                     1-gc:nijk(3)+gc, &
                                                                                     1:self%nb_proc))
   if (allocated(self%field_name)) deallocate(self%field_name) ; allocate(self%field_name(1:nv))
   self%emin = 0._R8P
   ! shift in x
   self%emin(1,1) = (nijk(1)-1) * dxyz(1) * self%myrank * self%nb_proc
   do b=2, self%nb_proc
      self%emin(1,b) = self%emin(1,b-1) + (nijk(1)-1)*dxyz(1)
   enddo
   do b=1, self%nb_proc
      do k=1, nijk(3)
         do j=1, nijk(2)
            do i=1, nijk(1)
               do v=1, nv
                  self%field(v,i,j,k,b) = (v + i + j + k + b) * 1._R8P
               enddo
            enddo
         enddo
      enddo
   enddo
   do v=1, nv
      self%field_name(v) = trim(adjustl(field_name(v)))
   enddo
   endsubroutine initialize
endmodule motion_test_domain_object

program motion_write_xdmf_file_test
!< MOTIOn test: write XDMF file test.
use motion
use motion_test_domain_object
use penf
use mpi

implicit none

type(domain_object) :: a_domain       !< A domain to play with.
logical             :: test_passed(1) !< List of passed tests.

call a_domain%initialize(nb=6_I4P,                                               &
                         nv=3_I4P,                                               &
                         gc=2_I4P,                                               &
                         nijk=[16_I8P,16_I8P,16_I8P],                            &
                         dxyz=[1._R8P,1._R8P,1._R8P],                            &
                         field_name=['pressure   ','temperature','density    '], &
                         time=0.56_R8P)

call write_hdf5_xdmf(basename='motion_write_xdmf_file_test',      domain=a_domain)
call write_xh5f(     basename='motion_write_xdmf_file_test-xh5f', domain=a_domain)
call MPI_FINALIZE(a_domain%error)

! test_passed = xdmf%error == 0
! print "(A,L1)", new_line('a')//'Are all tests passed? ', all(test_passed)
contains
   subroutine write_hdf5_xdmf(basename, domain)
   !< Write HDF5/XDMF files separately without XH5F file handler.
   character(*),        intent(in) :: basename      !< Basename of HDF5/XDMF file names.
   type(domain_object), intent(in) :: domain        !< Domain to be saved.
   character(:), allocatable       :: filename_hdf5 !< File name of HDF5 file.
   character(:), allocatable       :: filename_xdmf !< File name of XDMF file.
   type(hdf5_file_object)          :: hdf5          !< HDF5 file handler.
   type(xdmf_file_object)          :: xdmf          !< XDMF file handler.
   integer(I4P)                    :: b, v          !< Counter.

   associate(myrank=>domain%myrank, mynb=>domain%mynb, nijk=>domain%nijk, field=>domain%field, field_name=>domain%field_name, &
             emin=>domain%emin, dxyz=>domain%dxyz)
   filename_hdf5 = trim(adjustl(basename))//'-mpi_'//trim(strz(myrank,2))//'.h5'
   filename_xdmf = trim(adjustl(basename))//'-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf'
   call hdf5%open_file(filename=filename_hdf5)
   call xdmf%open_file(filename=filename_xdmf)
   call xdmf%open_domain_tag
   call xdmf%open_grid_tag(grid_name='blocks', grid_type=XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xdmf%open_grid_tag(grid_name='mpi_'//trim(strz(myrank,2)), grid_type=XDMF_GRID_TYPE_COLLECTION)
   do b=1, domain%nb_proc
      call hdf5%open_dspace(dataspace_type=HDF5_DATASPACE_TYPE_SIMPLE, nijk=nijk)
      call xdmf%open_grid_tag(grid_name='block_'//trim(strz(mynb(1)-1+b,2)))
      call xdmf%open_geometry_tag(geometry_type=XDMF_GEOMETRY_TYPE_ODXYZ)
      call xdmf%write_time_tag(time_value=trim(str(domain%time)))
      call xdmf%write_dataitem_tag(content          = trim(str([emin(3,b),emin(2,b),emin(1,b)],separator=' ')), &
                                   item_dimensions  = '3',                                                      &
                                   number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                          &
                                   number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                         &
                                   number_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%write_dataitem_tag(content          = trim(str([dxyz(3),dxyz(2),dxyz(1)],separator=' ')), &
                                   item_dimensions  = '3',                                                &
                                   number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                    &
                                   number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                   &
                                   number_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%close_geometry_tag
      call xdmf%write_topology_tag(topology_type=XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                                   topology_dimensions=trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')))
      call xdmf%open_attribute_tag(attribute_name   = 'Time',                &
                                   attribute_center = XDMF_ATTR_CENTER_GRID, &
                                   attribute_type   = XDMF_ATTR_TYPE_SCALAR)
      call xdmf%write_dataitem_tag(content          = trim(str(domain%time)),          &
                                   item_dimensions  = '1',                             &
                                   number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT, &
                                   number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,&
                                   number_format    = XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%close_attribute_tag
      do v=1, domain%nv
         call hdf5%save_dataset(dset_name='block_'//trim(strz(mynb(1)-1+b,2))//'-'//field_name(v)%chars(), &
                                nijk=nijk, dset=field(v,1:nijk(1),1:nijk(2),1:nijk(3),b))
         call xdmf%open_attribute_tag(attribute_name   = trim(adjustl(field_name(v)%chars())), &
                                      attribute_center = XDMF_ATTR_CENTER_CELL,                &
                                      attribute_type   = XDMF_ATTR_TYPE_SCALAR)
         call xdmf%write_dataitem_tag(content          = hdf5%filename//':'//'block_'//trim(strz(mynb(1)-1+b,2))// &
                                                         '-'//field_name(v)%chars(),                               &
                                      item_dimensions  = trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')),       &
                                      number_type      = XDMF_DATAITEM_NUMBER_TYPE_FLOAT,                          &
                                      number_precision = XDMF_DATAITEM_NUMBER_PRECISION_8,                         &
                                      number_format    = XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call xdmf%close_attribute_tag
      enddo
      call hdf5%close_dspace
      call xdmf%close_grid_tag
   enddo
   call xdmf%close_grid_tag
   call xdmf%close_grid_tag(grid_type=XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call hdf5%close_file
   call xdmf%close_domain_tag
   call xdmf%close_file
   endassociate
   endsubroutine write_hdf5_xdmf

   subroutine write_xh5f(basename, domain)
   !< Write HDF5/XDMF files by XH5F handler.
   character(*),        intent(in) :: basename      !< Basename of HDF5/XDMF file names.
   type(domain_object), intent(in) :: domain        !< Domain to be saved.
   character(:), allocatable       :: filename_hdf5 !< File name of HDF5 file.
   character(:), allocatable       :: filename_xdmf !< File name of XDMF file.
   type(xh5f_file_object)          :: xh5f          !< XH5F file handler.
   integer(I4P)                    :: b, v          !< Counter.

   associate(myrank=>domain%myrank, mynb=>domain%mynb, nijk=>domain%nijk, field=>domain%field, field_name=>domain%field_name)
   filename_hdf5 = trim(adjustl(basename))//'-mpi_'//trim(strz(myrank,2))//'.h5'
   filename_xdmf = trim(adjustl(basename))//'-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf'
   call xh5f%open_file(filename_hdf5=filename_hdf5, filename_xdmf=filename_xdmf)
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
   endassociate
   endsubroutine write_xh5f
endprogram motion_write_xdmf_file_test
