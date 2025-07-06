   !< MOTIOn, XH5F load/save block field agnostic impementation
   subroutine LOAD_BLOCK_FIELD_KKP_0D(self, field, xdmf_field_name, hdf5_field_name)
   !< Load field in block, kind KKP, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self            !< File handler.
   VARTYPE(KKP),            intent(inout)        :: field           !< Field.
   character(*),            intent(in), optional :: xdmf_field_name !< Field name in XDMF file.
   character(*),            intent(in), optional :: hdf5_field_name !< Field name in HDF5 file.
   character(:), allocatable                     :: field_format    !< Field format, local var.
   character(:), allocatable                     :: field_name      !< Field name.

   select case(field_format)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF)
      call self%hdf5%load_dataset(dset_name=field_name, dset=field)
   case(XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_XML)
   case default
      write(stderr, '(A)') 'error: load block field error, field_format "'//field_format//'" unknown'
      call MPI_FINALIZE(self%error)
      stop
   endselect
   endsubroutine LOAD_BLOCK_FIELD_KKP_0D 
 
   subroutine LOAD_BLOCK_FIELD_KKP_3D(self, nd, field, xdmf_field_name, hdf5_field_name)
   !< Load field in block, kind KKP, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self            !< File handler.
   integer(HSIZE_T),        intent(in)           :: nd(:)           !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout)        :: field(:,:,:)    !< Field.
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
   endsubroutine LOAD_BLOCK_FIELD_KKP_3D 

   subroutine LOAD_BLOCK_FIELD_KKP_4D(self, nd, field, xdmf_field_name, hdf5_field_name)
   !< Load field in block, kind KKP, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self            !< File handler.
   integer(HSIZE_T),        intent(in)           :: nd(:)           !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout)        :: field(:,:,:,:)  !< Field.
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
   endsubroutine LOAD_BLOCK_FIELD_KKP_4D

   subroutine SAVE_BLOCK_FIELD_KKP_0D(self, xdmf_field_name, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind KKP, rank 0D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   VARTYPE(KKP),            intent(in)           :: field            !< Field.
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
   endsubroutine SAVE_BLOCK_FIELD_KKP_0D
 
   subroutine SAVE_BLOCK_FIELD_KKP_3D(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind KKP, rank 3D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)           :: field(:,:,:)     !< Field.
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
   endsubroutine SAVE_BLOCK_FIELD_KKP_3D
 
   subroutine SAVE_BLOCK_FIELD_KKP_4D(self, xdmf_field_name, nd, field, field_center, field_format, hdf5_field_name)
   !< Save field in block, kind KKP, rank 4D.
   class(xh5f_file_object), intent(inout)        :: self             !< File handler.
   character(*),            intent(in)           :: xdmf_field_name  !< Field name in XDMF file.
   integer(HSIZE_T),        intent(in)           :: nd(:)            !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)           :: field(:,:,:,:)   !< Field.
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
   endsubroutine SAVE_BLOCK_FIELD_KKP_4D
