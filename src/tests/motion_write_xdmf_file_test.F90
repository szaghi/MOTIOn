!< MOTIOn test: write XDMF file test.

module motion_test_domain_object
!< MOTIOn test: prototype of domain to be saved.
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use motion
use penf
use stringifor
use hdf5
use mpi

implicit none
private
public :: domain_object

type :: domain_object
   !< Prototype of domain to be saved.
   !< For the sake of simplicity all domain blocks have the same cells number and spacing.
   character(:), allocatable :: topology                                !< Topology.
   integer(I4P)              :: procs_number=1_I4P                      !< Number of MPI processes.
   integer(I4P)              :: myrank=0_I4P                            !< MPI ID process.
   integer(I4P)              :: nb_proc=0_I4P                           !< Number of blocks for each MPI process.
   integer(I4P)              :: nb=0_I4P                                !< Total number of blocks (all processes).
   integer(I4P)              :: mynb(2)=[0_I4P,0_I4P]                   !< Number of blocks of current process [start-b,end-b].
   integer(I4P)              :: nv=0_I4P                                !< Number of fields variables (scalar, vector...).
   integer(I4P)              :: nvscalar=0_I4P                          !< Last index of scalar fields variables.
   integer(I4P)              :: gc=0_I4P                                !< Ghost cells number.
   real(R8P)                 :: time=0._R8P                             !< Time of domain solution.
   integer(HSIZE_T)          :: nijk(3)=[0_HSIZE_T,0_HSIZE_T,0_HSIZE_T] !< Blocks dimensions.
   real(R8P)                 :: dxyz(3)=[0._R8P,0._R8P,0._R8P]          !< Blocks space steps for cartesian uniform grid.
   real(R8P),    allocatable :: dx(:),dy(:),dz(:)                       !< Blocks space steps for cartesian         grid.
   real(R8P),    allocatable :: emin(:,:)                               !< Blocks minimum extents [3,nb].
   real(R8P),    allocatable :: x(:,:),y(:,:),z(:,:)                    !< Nodes coordinates for cartesian grid [0:nijk(ijk),nb].
   real(R8P),    allocatable :: nodes(:,:,:,:,:)                        !< Nodes coordinates for curvilinear grid [3,0:nijk,nb].
   real(R8P),    allocatable :: field(:,:,:,:,:)                        !< Fields [1:nv,1-gc:ni+gc,1-gc:nj+gc,1-gc:nk+gc,1:nb].
   type(string), allocatable :: field_name(:)                           !< Fields names [1:nv].
   integer(I4P)              :: error=0_I4P                             !< Error status.
   contains
      procedure, pass(self) :: initialize !< Initialize domain.
endtype domain_object

contains
   subroutine initialize(self, topology, nb, nv, nvscalar, gc, nijk, field_name, time, dxyz, dx, dy, dz)
   !< Initialize domain.
   class(domain_object), intent(inout)        :: self               !< Domain.
   character(*),         intent(in)           :: topology           !< Topology.
   integer(I4P),         intent(in)           :: nb                 !< Number of blocks.
   integer(I4P),         intent(in)           :: nv                 !< Number of fields variables.
   integer(I4P),         intent(in)           :: nvscalar           !< Last index of scalar fields variables.
   integer(I4P),         intent(in)           :: gc                 !< Ghost cells number.
   integer(HSIZE_T),     intent(in)           :: nijk(3)            !< Blocks dimensions.
   character(*),         intent(in)           :: field_name(nv)     !< Fields names.
   real(R8P),            intent(in)           :: time               !< Time of domain solution.
   real(R8P),            intent(in), optional :: dxyz(3)            !< Blocks space steps for cartesian uniform grids.
   real(R8P),            intent(in), optional :: dx(:),dy(:),dz(:)  !< Blocks space steps for cartesian         grids.
   integer(I4P)                               :: i,j,k,b,v          !< Counter.
   real(R8P)                                  :: lx                 !< Single block x length for cartesian grid.
   logical                                    :: is_mpi_initialized !< MPI env status.

   ! reset domain
   select type(self)
   type is(domain_object)
      self = domain_object()
   endselect
   call MPI_INITIALIZED(is_mpi_initialized, self%error)
   if (.not.is_mpi_initialized) call MPI_INIT(self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   self%nb_proc  = nb / self%procs_number
   self%mynb     = [self%nb_proc * self%myrank + 1, self%nb_proc * (self%myrank + 1)]
   self%topology = trim(adjustl(topology))
   self%nb       = nb
   self%nv       = nv
   self%nvscalar = nvscalar
   self%gc       = gc
   self%nijk     = nijk
   self%time     = time
   ! inizialize grid
   if (allocated(self%emin)) deallocate(self%emin) ; allocate(self%emin(1:3,1:self%nb_proc))
   self%emin = 0._R8P
   select case(self%topology)
   case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN)
      if (present(dx).and.present(dy).and.present(dz)) then
         self%dx = dx
         self%dy = dy
         self%dz = dz
         if (allocated(self%x)) deallocate(self%x) ; allocate(self%x(0:nijk(1),1:self%nb_proc))
         if (allocated(self%y)) deallocate(self%y) ; allocate(self%y(0:nijk(2),1:self%nb_proc))
         if (allocated(self%z)) deallocate(self%z) ; allocate(self%z(0:nijk(3),1:self%nb_proc))
         lx = 0._R8P
         do i=1, nijk(1)
            lx = lx + dx(i)
         enddo
         ! shift in x
         self%emin(1,1) = lx * self%myrank * self%nb_proc
         do b=2, self%nb_proc
            self%emin(1,b) = self%emin(1,b-1) + lx
         enddo
         ! compute nodes coordinates
         do b=1, self%nb_proc
            self%x(0,b) = self%emin(1,b)
            do i=1, nijk(1)
               self%x(i,b) = self%x(i-1,b) + dx(i)
            enddo
            self%y(0,b) = self%emin(2,b)
            do j=1, nijk(2)
               self%y(j,b) = self%y(j-1,b) + dy(j)
            enddo
            self%z(0,b) = self%emin(3,b)
            do k=1, nijk(3)
               self%z(k,b) = self%z(k-1,b) + dz(k)
            enddo
         enddo
      else
         write(stderr, '(A)') 'error: block of type "'//self%topology//'" need dx,dy,dz to be passed'
         call MPI_FINALIZE(self%error)
         stop
      endif
   case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM)
      if (present(dxyz)) then
         self%dxyz = dxyz
         ! shift in x
         self%emin(1,1) = (nijk(1)-1) * dxyz(1) * self%myrank * self%nb_proc
         do b=2, self%nb_proc
            self%emin(1,b) = self%emin(1,b-1) + (nijk(1)-1)*dxyz(1)
         enddo
      else
         write(stderr, '(A)') 'error: block of type "'//self%topology//'" need dxyz to be passed'
         call MPI_FINALIZE(self%error)
         stop
      endif
   case(XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR)
      if (present(dx).and.present(dy).and.present(dz)) then
         self%dx = dx
         self%dy = dy
         self%dz = dz
         if (allocated(self%nodes)) deallocate(self%nodes)
         allocate(self%nodes(1:3,0:nijk(1),0:nijk(2),0:nijk(3),1:self%nb_proc))
         lx = 0._R8P
         do i=1, nijk(1)
            lx = lx + dx(i)
         enddo
         ! shift in x
         self%emin(1,1) = lx * self%myrank * self%nb_proc
         do b=2, self%nb_proc
            self%emin(1,b) = self%emin(1,b-1) + lx
         enddo
         ! compute nodes coordinates
         do b=1, self%nb_proc
            self%nodes(1,:,:,:,b) = self%emin(1,b)
            self%nodes(2,:,:,:,b) = self%emin(2,b)
            self%nodes(3,:,:,:,b) = self%emin(3,b)
            do k=0, nijk(3)
            do j=0, nijk(2)
            do i=0, nijk(1)
               if (i>0) self%nodes(1,i,j,k,b) = self%nodes(1,i-1,0,0,b) + dx(i)
               if (j>0) self%nodes(2,i,j,k,b) = self%nodes(2,0,j-1,0,b) + dy(j)
               if (k>0) self%nodes(3,i,j,k,b) = self%nodes(3,0,0,k-1,b) + dz(k)
            enddo
            enddo
            enddo
            ! do j=1, nijk(2)
            !    self%nodes(2,:,j,:,b) = self%nodes(2,0,j-1,0,b) + dy(j)
            ! enddo
            ! do k=1, nijk(3)
            !    self%nodes(3,:,:,k,b) = self%nodes(3,0,0,k-1,b) + dz(k)
            ! enddo
         enddo
      else
         write(stderr, '(A)') 'error: block of type "'//self%topology//'" need dx,dy,dz to be passed'
         call MPI_FINALIZE(self%error)
         stop
      endif
   case default
      write(stderr, '(A)') 'error: block of type "'//self%topology//'" is unknown'
      call MPI_FINALIZE(self%error)
      stop
   endselect
   ! initialize fields
   if (allocated(self%field     )) deallocate(self%field     ) ; allocate(self%field(1:nv,            &
                                                                                     1-gc:nijk(1)+gc, &
                                                                                     1-gc:nijk(2)+gc, &
                                                                                     1-gc:nijk(3)+gc, &
                                                                                     1:self%nb_proc))
   if (allocated(self%field_name)) deallocate(self%field_name) ; allocate(self%field_name(1:nv))
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
use hdf5
use mpi

implicit none

type(domain_object) :: domain_cartesian_uniform !< Cartesian uniform domain.
type(domain_object) :: domain_cartesian         !< Cartesian domain.
type(domain_object) :: domain_curvilinear       !< Curvilinear domain.
integer(I4P)        :: error                    !< Error status.
logical             :: test_passed(1)           !< List of passed tests.

call domain_cartesian_uniform%initialize(topology=XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM,  &
                                         nb=6_I4P,                                               &
                                         nv=33_I4P,                                              &
                                         nvscalar=3_I4P,                                         &
                                         gc=2_I4P,                                               &
                                         nijk=[16_HSIZE_T,16_HSIZE_T,16_HSIZE_T],                &
                                         field_name=['pressure   ','temperature','density    ',  &
                                                     'velocity_x ','velocity_y ','velocity_z ',  &
                                                     'stress6_1  ','stress6_2  ','stress6_3  ',  &
                                                     'stress6_4  ','stress6_5  ','stress6_6  ',  &
                                                     'stress_1   ','stress_2   ','stress_3   ',  &
                                                     'stress_4   ','stress_5   ','stress_6   ',  &
                                                     'stress_7   ','stress_8   ','stress_9   ',  &
                                                     'matrix_1   ','matrix_2   ','matrix_3   ',  &
                                                     'matrix_4   ','matrix_5   ','matrix_6   ',  &
                                                     'matrix_7   ','matrix_8   ','matrix_9   ',  &
                                                     'matrix_10  ','matrix_11  ','matrix_12  '], &
                                         time=0.56_R8P,                                          &
                                         dxyz=[1._R8P,1._R8P,1._R8P])

call write_hdf5_xdmf(basename='motion_write_xdmf_file_test-cartesian_unifrom',      domain=domain_cartesian_uniform)
call write_xh5f(     basename='motion_write_xdmf_file_test-cartesian_unifrom-xh5f', domain=domain_cartesian_uniform)

call domain_cartesian%initialize(topology=XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN,          &
                                 nb=6_I4P,                                               &
                                 nv=33_I4P,                                              &
                                 nvscalar=3_I4P,                                         &
                                 gc=2_I4P,                                               &
                                 nijk=[3_HSIZE_T,4_HSIZE_T,5_HSIZE_T],                   &
                                 field_name=['pressure   ','temperature','density    ',  &
                                             'velocity_x ','velocity_y ','velocity_z ',  &
                                             'stress6_1  ','stress6_2  ','stress6_3  ',  &
                                             'stress6_4  ','stress6_5  ','stress6_6  ',  &
                                             'stress_1   ','stress_2   ','stress_3   ',  &
                                             'stress_4   ','stress_5   ','stress_6   ',  &
                                             'stress_7   ','stress_8   ','stress_9   ',  &
                                             'matrix_1   ','matrix_2   ','matrix_3   ',  &
                                             'matrix_4   ','matrix_5   ','matrix_6   ',  &
                                             'matrix_7   ','matrix_8   ','matrix_9   ',  &
                                             'matrix_10  ','matrix_11  ','matrix_12  '], &
                                 time=0.56_R8P,                                          &
                                 dx=[0.1_R8P, 0.2_R8P, 0.3_R8P],                         &
                                 dy=[0.1_R8P, 0.2_R8P, 0.3_R8P, 0.4_R8P],                &
                                 dz=[0.1_R8P, 0.2_R8P, 0.3_R8P, 0.4_R8P, 0.5_R8P])
call write_xh5f(basename='motion_write_xdmf_file_test-cartesian-xh5f', domain=domain_cartesian)

call domain_curvilinear%initialize(topology=XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR,        &
                                   nb=6_I4P,                                               &
                                   nv=33_I4P,                                              &
                                   nvscalar=3_I4P,                                         &
                                   gc=2_I4P,                                               &
                                   nijk=[3_HSIZE_T,4_HSIZE_T,5_HSIZE_T],                   &
                                   field_name=['pressure   ','temperature','density    ',  &
                                               'velocity_x ','velocity_y ','velocity_z ',  &
                                               'stress6_1  ','stress6_2  ','stress6_3  ',  &
                                               'stress6_4  ','stress6_5  ','stress6_6  ',  &
                                               'stress_1   ','stress_2   ','stress_3   ',  &
                                               'stress_4   ','stress_5   ','stress_6   ',  &
                                               'stress_7   ','stress_8   ','stress_9   ',  &
                                               'matrix_1   ','matrix_2   ','matrix_3   ',  &
                                               'matrix_4   ','matrix_5   ','matrix_6   ',  &
                                               'matrix_7   ','matrix_8   ','matrix_9   ',  &
                                               'matrix_10  ','matrix_11  ','matrix_12  '], &
                                   time=0.56_R8P,                                          &
                                   dx=[0.1_R8P, 0.2_R8P, 0.3_R8P],                         &
                                   dy=[0.1_R8P, 0.2_R8P, 0.3_R8P, 0.4_R8P],                &
                                   dz=[0.1_R8P, 0.2_R8P, 0.3_R8P, 0.4_R8P, 0.5_R8P])
call write_xh5f(basename='motion_write_xdmf_file_test-curvilinear-xh5f', domain=domain_curvilinear)

call MPI_FINALIZE(error)

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
   call xdmf%open_grid_tag(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xdmf%open_grid_tag(grid_name='mpi_'//trim(strz(myrank,2)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
   do b=1, domain%nb_proc
      call hdf5%open_dspace(dataspace_type=HDF5_PARAMETERS%HDF5_DATASPACE_TYPE_SIMPLE, nd=nijk)
      call xdmf%open_grid_tag(grid_name='block_'//trim(strz(mynb(1)-1+b,2)))
      call xdmf%open_geometry_tag(geometry_type=XDMF_PARAMETERS%XDMF_GEOMETRY_TYPE_ODXYZ)
      call xdmf%write_time_tag(time_value=trim(str(domain%time)))
      call xdmf%write_dataitem_tag(content          = trim(str([emin(3,b),emin(2,b),emin(1,b)],separator=' ')), &
                                   item_dimensions  = '3',                                                      &
                                   number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,          &
                                   number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,         &
                                   number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%write_dataitem_tag(content          = trim(str([dxyz(3),dxyz(2),dxyz(1)],separator=' ')), &
                                   item_dimensions  = '3',                                                &
                                   number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,    &
                                   number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,   &
                                   number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%close_geometry_tag
      call xdmf%write_topology_tag(topology_type=XDMF_PARAMETERS%XDMF_TOPOLOGY_TYPE_3DCORECTMESH, &
                                   topology_dimensions=trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')))
      call xdmf%open_attribute_tag(attribute_name   = 'Time',                                &
                                   attribute_center = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID, &
                                   attribute_type   = XDMF_PARAMETERS%XDMF_ATTR_TYPE_SCALAR)
      call xdmf%write_dataitem_tag(content          = trim(str(domain%time)),                          &
                                   item_dimensions  = '1',                                             &
                                   number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT, &
                                   number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,&
                                   number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      call xdmf%close_attribute_tag
      do v=1, domain%nv
         call hdf5%save_dataset(dset_name='block_'//trim(strz(mynb(1)-1+b,2))//'-'//field_name(v)%chars(), &
                                nd=nijk, dset=field(v,1:nijk(1),1:nijk(2),1:nijk(3),b))
         call xdmf%open_attribute_tag(attribute_name   = trim(adjustl(field_name(v)%chars())),  &
                                      attribute_center = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL, &
                                      attribute_type   = XDMF_PARAMETERS%XDMF_ATTR_TYPE_SCALAR)
         call xdmf%write_dataitem_tag(content          = hdf5%filename//':'//'block_'//trim(strz(mynb(1)-1+b,2))// &
                                                         '-'//field_name(v)%chars(),                               &
                                      item_dimensions  = trim(str([nijk(3),nijk(2),nijk(1)],separator=' ')),       &
                                      number_type      = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_TYPE_FLOAT,          &
                                      number_precision = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_PRECISION_8,         &
                                      number_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
         call xdmf%close_attribute_tag
      enddo
      call hdf5%close_dspace
      call xdmf%close_grid_tag
   enddo
   call xdmf%close_grid_tag
   call xdmf%close_grid_tag(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
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

   associate(nvs=>domain%nvscalar, mynb=>domain%mynb, nijk=>domain%nijk, field=>domain%field, field_name=>domain%field_name)
   filename_hdf5 = trim(adjustl(basename))//'-mpi_'//trim(strz(domain%myrank,2))//'.h5'
   filename_xdmf = trim(adjustl(basename))//'-mpi_procs_'//trim(strz(domain%procs_number,2))//'.xdmf'
   call xh5f%open_file(filename_hdf5=filename_hdf5, filename_xdmf=filename_xdmf)
   call xh5f%open_grid(grid_name='blocks', grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%open_grid(grid_name='mpi_'//trim(strz(domain%myrank,2)), grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION)
   do b=1, domain%nb_proc
      select case(domain%topology)
      case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN_UNIFORM)
         call xh5f%open_block(block_type = domain%topology,                     &
                              block_name = 'block_'//trim(strz(mynb(1)-1+b,2)), & ! global block numeration
                              nijk       = nijk,                                &
                              emin       = domain%emin(:,b),                    &
                              dxyz       = domain%dxyz,                         &
                              time       = domain%time)
      case(XH5F_PARAMETERS%XH5F_BLOCK_CARTESIAN)
         call xh5f%open_block(block_type = domain%topology,                     &
                              block_name = 'block_'//trim(strz(mynb(1)-1+b,2)), & ! global block numeration
                              nijk       = nijk,                                &
                              x          = domain%x(0:,b),                      &
                              y          = domain%y(0:,b),                      &
                              z          = domain%z(0:,b),                      &
                              time       = domain%time)
      case(XH5F_PARAMETERS%XH5F_BLOCK_CURVILINEAR)
         call xh5f%open_block(block_type = domain%topology,                     &
                              block_name = 'block_'//trim(strz(mynb(1)-1+b,2)), & ! global block numeration
                              nijk       = nijk,                                &
                              nodes      = domain%nodes(1:,0:,0:,0:,b),         &
                              time       = domain%time)
      endselect
      call xh5f%save_block_field(xdmf_field_name = 'Time',                                &
                                 field           = domain%time,                           &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_GRID, &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
      ! scalar fields
      do v=1, nvs
         call xh5f%save_block_field(xdmf_field_name = field_name(v)%chars(),                           &
                                    nd              = nijk,                                            &
                                    field           = field(v,1:nijk(1),1:nijk(2),1:nijk(3),b),        &
                                    field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,           &
                                    field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                    hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//field_name(v)%chars())
      enddo
      ! vector field
      call xh5f%save_block_field(xdmf_field_name = 'velocity',                                         &
                                 nd              = [3_HSIZE_T,nijk(1),nijk(2),nijk(3)],                &
                                 field           = field(nvs+1:nvs+3,1:nijk(1),1:nijk(2),1:nijk(3),b), &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,              &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,    &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//'velocity')
      ! tensor6 field
      call xh5f%save_block_field(xdmf_field_name = 'stress6',                                          &
                                 nd              = [6_HSIZE_T,nijk(1),nijk(2),nijk(3)],                &
                                 field           = field(nvs+4:nvs+9,1:nijk(1),1:nijk(2),1:nijk(3),b), &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,              &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,    &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//'stress6')
      ! tensor field
      call xh5f%save_block_field(xdmf_field_name = 'stress',                                             &
                                 nd              = [9_HSIZE_T,nijk(1),nijk(2),nijk(3)],                  &
                                 field           = field(nvs+10:nvs+18,1:nijk(1),1:nijk(2),1:nijk(3),b), &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,                &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,      &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//'stress')
      ! matrix field
      call xh5f%save_block_field(xdmf_field_name = 'matrix',                                        &
                                 nd              = [12_HSIZE_T,nijk(1),nijk(2),nijk(3)],            &
                                 field           = field(nvs+19:,1:nijk(1),1:nijk(2),1:nijk(3),b),  &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,           &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF, &
                                 hdf5_field_name = 'block_'//trim(strz(mynb(1)-1+b,2))//'-'//'matrix')
      call xh5f%close_block
   enddo
   call xh5f%close_grid
   call xh5f%close_grid(grid_type=XDMF_PARAMETERS%XDMF_GRID_TYPE_COLLECTION_ASYNC)
   call xh5f%close_file
   endassociate
   endsubroutine write_xh5f
endprogram motion_write_xdmf_file_test
