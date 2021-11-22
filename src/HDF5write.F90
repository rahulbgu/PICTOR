module HDF5write
	use parameters
	use vars
	use HDF5
contains 
!---------------------------------------------------------------------------------------------------------------------------------
!    Some high level subroutines to write data in HDF5 files 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine HDF5writeINT(fid,dspace_id,varname,var)
     INTEGER(HID_T)                 :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim1
     character(len=*)               :: varname
     integer                        :: var,err
     data_dim1(1)=1
     call h5dcreate_f(fid,varname,H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,data_dim1,err)
     call h5dclose_f(dset_id,err) 
end subroutine HDF5writeINT

subroutine HDF5writeArrINT(fid,dspace_id,sizex,varname,var)
     INTEGER(HID_T)                 :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim1
     character(len=*)               :: varname
     integer                        :: err,sizex
     integer, dimension(sizex)      :: var
     data_dim1(1)=sizex
    call h5dcreate_f(fid,varname,H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,data_dim1,err)
     call h5dclose_f(dset_id,err) 
end subroutine HDF5writeArrINT

subroutine HDF5writeINT3(fid,dspace_id,sizex,sizey,sizez,varname,var)
     INTEGER(HID_T)                      :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(3)      :: data_dim
     character(len=*)                    :: varname
     integer                             :: sizex,sizey,sizez,err
     integer,dimension(sizex,sizey,sizez):: var
     
     data_dim(1)=sizex
     data_dim(2)=sizey
     data_dim(3)=sizez
    call h5dcreate_f(fid,varname,H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,data_dim,err)
     call h5dclose_f(dset_id,err) 
end subroutine HDF5writeINT3

subroutine WriteArrReal1(FileName,FileID,sizex,var,varname,extn)
     character(len=*), optional, intent(in) :: FileName
     integer(HID_T), optional,   intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex
     real     , dimension(sizex),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     integer, optional,          intent(in) :: extn
     INTEGER(HID_T)                         :: dspace_id,dset_id,fid_this 
     INTEGER(HSIZE_T),dimension(1)          :: data_dim
     integer                                :: rank,err
     
     
     data_dim(1)=sizex
     if(present(FileID)) then ! FileName must be present otherwise otherwise   
          fid_this=FileID
     else 
          if(present(extn)) then
                call GetFileID(FileName,extn,fid_this)
          else 
               call GetFileID(FileName,0,fid_this) ! by default filename ends with _TimeStep
          end if
     end if
     rank=1
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(fid_this,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_REAL,real(var),data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)     
     
end subroutine WriteArrReal1

subroutine WriteArrReal1_DP(FileName,FileID,sizex,var,varname,extn) ! same as above, but in double precison 
     character(len=*), optional, intent(in) :: FileName
     integer(HID_T), optional,   intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex
     real(dbpsn),dimension(sizex),intent(in):: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     integer, optional,          intent(in) :: extn
     INTEGER(HID_T)                         :: dspace_id,dset_id,fid_this 
     INTEGER(HSIZE_T),dimension(1)          :: data_dim
     integer                                :: rank,err
     
     
     data_dim(1)=sizex
     if(present(FileID)) then ! FileName must be present otherwise 
          fid_this=FileID
     else 
          if(present(extn)) then
                call GetFileID(FileName,extn,fid_this)
          else 
               call GetFileID(FileName,0,fid_this) ! by default filename ends with _TimeStep
          end if
     end if
     rank=1
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(fid_this,varname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)     
     
end subroutine WriteArrReal1_DP

subroutine WriteArrReal2(FileName,FileID,sizex,sizey,var,varname,extn)
     character(len=*), optional, intent(in) :: FileName
     integer(HID_T), optional,   intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex,sizey
     real,dimension(sizex,sizey),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     integer, optional,          intent(in) :: extn
     INTEGER(HID_T)                         :: dspace_id,dset_id,fid_this 
     INTEGER(HSIZE_T),dimension(2)          :: data_dim
     integer                                :: rank,err
     
     
     data_dim(1)=sizex
     data_dim(2)=sizey
     if(present(FileID)) then ! FileName must be present otherwise otherwise   
          fid_this=FileID
     else 
          if(present(extn)) then
                call GetFileID(FileName,extn,fid_this)
          else 
               call GetFileID(FileName,0,fid_this) ! by default filename ends with _TimeStep
          end if
     end if
     rank=2
    call h5screate_simple_f(rank,data_dim,dspace_id,err)
    call h5dcreate_f(fid_this,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
    call h5dwrite_f(dset_id,H5T_NATIVE_REAL,real(var),data_dim,err)
    call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)          
end subroutine WriteArrReal2

subroutine WriteArrReal2_DP(FileName,FileID,sizex,sizey,var,varname,extn)
     character(len=*), optional, intent(in) :: FileName
     integer(HID_T), optional,   intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex,sizey
     real(dbpsn),dimension(sizex,sizey),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     integer, optional,          intent(in) :: extn
     INTEGER(HID_T)                         :: dspace_id,dset_id,fid_this 
     INTEGER(HSIZE_T),dimension(2)          :: data_dim
     integer                                :: rank,err
     
     
     data_dim(1)=sizex
     data_dim(2)=sizey
     if(present(FileID)) then ! FileName must be present otherwise otherwise   
          fid_this=FileID
     else 
          if(present(extn)) then
                call GetFileID(FileName,extn,fid_this)
          else 
               call GetFileID(FileName,0,fid_this) ! by default filename ends with _TimeStep
          end if
     end if
     rank=2
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(fid_this,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)          
end subroutine WriteArrReal2_DP

subroutine GetFileID(FileName,extn,fid)
     character(len=*), intent(in) :: FileName
     integer, intent(in)          :: extn
     integer(HID_T), intent(inout):: fid
     character (len=1024)         :: FileNameThis
     logical                      :: file_exists
     integer                      :: h5err
     
     
     if(extn.eq.1) then 
         write(FileNameThis,'(I0)') t
        FileNameThis=trim(data_folder)//"/"//trim(FileName)//"_"//FileNameThis
     else 
         FileNameThis=trim(data_folder)//"/"//FileName
     end if     
          
     inquire(File=FileNameThis,exist=file_exists)
     if(file_exists) then 
          call h5fopen_f(FileNameThis, H5F_ACC_RDWR_F,fid, h5err)
     else 
          call h5fcreate_f(FileNameThis, H5F_ACC_TRUNC_F,fid, h5err)
     end if
     
end subroutine GetFileID 
end module HDF5write