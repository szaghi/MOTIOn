   !< MOTIOn, HDF5 load/save dataset agnostic impementation
   subroutine LOAD_DATASET_KKP_0D(self, dset_name, dset)
   !< Load dataset in dataspace, kind KKP, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   VARTYPE(KKP),            intent(inout) :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.
   VARTYPE(KKP)                           :: dset_(1)  !< Dataset, local var.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset_, [1_HSIZE_T], self%error)
   dset = dset_(1)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_0D
 
   subroutine LOAD_DATASET_KKP_1D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_1D
 
   subroutine LOAD_DATASET_KKP_2D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_2D
 
   subroutine LOAD_DATASET_KKP_3D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_3D
 
   subroutine LOAD_DATASET_KKP_4D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_4D
 
   subroutine LOAD_DATASET_KKP_5D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_5D
 
   subroutine LOAD_DATASET_KKP_6D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_6D
  
    subroutine LOAD_DATASET_KKP_7D(self, dset_name, nd, dset)
   !< Load dataset in dataspace, kind KKP, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(inout) :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dopen_f(self%hdf5, trim(adjustl(dset_name)), dset_id, self%error)
   call h5dread_f( dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine LOAD_DATASET_KKP_7D

   subroutine SAVE_DATASET_KKP_0D(self, dset_name, dset)
   !< Save dataset in dataspace, kind KKP, rank 0D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   VARTYPE(KKP),            intent(in)    :: dset      !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, [dset], [1_HSIZE_T], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_0D
 
   subroutine SAVE_DATASET_KKP_1D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 1D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd        !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:)   !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, [nd], self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_1D
 
   subroutine SAVE_DATASET_KKP_2D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 2D.
   class(hdf5_file_object), intent(inout) :: self      !< File handler.
   character(*),            intent(in)    :: dset_name !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)     !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:) !< Dataset.
   integer(HID_T)                         :: dset_id   !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_2D
 
   subroutine SAVE_DATASET_KKP_3D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 3D.
   class(hdf5_file_object), intent(inout) :: self        !< File handler.
   character(*),            intent(in)    :: dset_name   !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)       !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id     !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_3D
 
   subroutine SAVE_DATASET_KKP_4D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 4D.
   class(hdf5_file_object), intent(inout) :: self          !< File handler.
   character(*),            intent(in)    :: dset_name     !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)         !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id       !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_4D
 
   subroutine SAVE_DATASET_KKP_5D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 5D.
   class(hdf5_file_object), intent(inout) :: self            !< File handler.
   character(*),            intent(in)    :: dset_name       !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)           !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id         !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_5D
 
   subroutine SAVE_DATASET_KKP_6D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 6D.
   class(hdf5_file_object), intent(inout) :: self              !< File handler.
   character(*),            intent(in)    :: dset_name         !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)             !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id           !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_6D
 
   subroutine SAVE_DATASET_KKP_7D(self, dset_name, nd, dset)
   !< Save dataset in dataspace, kind KKP, rank 7D.
   class(hdf5_file_object), intent(inout) :: self                !< File handler.
   character(*),            intent(in)    :: dset_name           !< Dataset name.
   integer(HSIZE_T),        intent(in)    :: nd(:)               !< Dataspace datasets dimensions.
   VARTYPE(KKP),            intent(in)    :: dset(:,:,:,:,:,:,:) !< Dataset.
   integer(HID_T)                         :: dset_id             !< Dataset identifier.

   call h5dcreate_f(self%hdf5, trim(adjustl(dset_name)), H5T_KKP, self%dspace_id, dset_id, self%error)
   call h5dwrite_f(dset_id, H5T_KKP, dset, nd, self%error)
   call h5dclose_f(dset_id, self%error)
   endsubroutine SAVE_DATASET_KKP_7D
