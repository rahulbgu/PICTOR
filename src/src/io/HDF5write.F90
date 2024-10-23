module HDF5write
	use parameters
	use vars
	use HDF5
	
	implicit none
	
	INTERFACE WriteReal
		module procedure WriteReal_SP, WriteReal_DP
	END INTERFACE 
	
	INTERFACE WriteArrReal1
		module procedure WriteArrReal1_SP, WriteArrReal1_DP 
	END INTERFACE
	
	INTERFACE WriteArrReal2
		module procedure WriteArrReal2_SP, WriteArrReal2_DP 
	END INTERFACE 
	
contains 
!---------------------------------------------------------------------------------------------------------------------------------
!    Some high level subroutines to write data in HDF5 files 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine WriteINT(fid,dspace_id,varname,var)
     INTEGER(HID_T)                 :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim1
     character(len=*)               :: varname
     integer                        :: var,err
     data_dim1(1)=1
     call h5dcreate_f(fid,varname,H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,data_dim1,err)
     call h5dclose_f(dset_id,err) 
end subroutine WriteINT

subroutine WriteINT1(fid,sizex,varname,var)
	 INTEGER(HID_T), intent(in)     :: fid 
	 character(len=*), intent(in)   :: varname
	 integer, intent(in)            :: sizex
     integer,dimension(:),intent(in):: var
     INTEGER(HID_T)                 :: dspace_id, dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim
     integer                        :: err, rank
	 
	 rank = 1
     data_dim(1)=sizex
	 call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(fid,varname,H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var(1:sizex),data_dim1,err)
     call h5dclose_f(dset_id,err) 
	 call h5sclose_f(dspace_id,err)  
end subroutine WriteINT1

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

!---------------------------------------------------------------------------------------------
! Write a single real varaible (SP = single precision, DP=double precision)
! use the interface WriteArrReal 
!---------------------------------------------------------------------------------------------
subroutine WriteReal_SP(fid,dspace_id,varname,var)
     INTEGER(HID_T)                 :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim1
     character(len=*)               :: varname
	 real                           :: var
     integer                        :: err
     data_dim1(1)=1
     call h5dcreate_f(fid,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_REAL,var,data_dim1,err)
     call h5dclose_f(dset_id,err) 
end subroutine WriteReal_SP

subroutine WriteReal_DP(fid,dspace_id,varname,var)
     INTEGER(HID_T)                 :: dspace_id,fid,dset_id 
     INTEGER(HSIZE_T), dimension(1) :: data_dim1
     character(len=*)               :: varname
	 real(dbpsn)                    :: var
     integer                        :: err
     data_dim1(1)=1
     call h5dcreate_f(fid,varname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,data_dim1,err)
     call h5dclose_f(dset_id,err) 
end subroutine WriteReal_DP

!---------------------------------------------------------------------------------------------
! One dimensional array of reals (SP = single precision, DP=double precision)
! use the interface WriteArrReal1 
!---------------------------------------------------------------------------------------------
subroutine WriteArrReal1_SP(FileID,sizex,var,varname)
     integer(HID_T),             intent(in) :: FileID  
     integer,                    intent(in) :: sizex
     real     , dimension(sizex),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     INTEGER(HID_T)                         :: dspace_id,dset_id
     INTEGER(HSIZE_T),dimension(1)          :: data_dim
     integer                                :: rank,err
      
     data_dim(1)=sizex
     rank=1
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(FileID,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_REAL,var(1:sizex),data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)     
     
end subroutine WriteArrReal1_SP

subroutine WriteArrReal1_DP(FileID,sizex,var,varname) ! same as above, but in double precison 
     integer(HID_T),             intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex
     real(dbpsn),dimension(sizex),intent(in):: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     INTEGER(HID_T)                         :: dspace_id,dset_id
     INTEGER(HSIZE_T),dimension(1)          :: data_dim
     integer                                :: rank,err
     
     data_dim(1)=sizex
     rank=1
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(FileID,varname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var(1:sizex),data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)     
     
end subroutine WriteArrReal1_DP

!---------------------------------------------------------------------------------------------
! Two dimensional array of reals (SP = single precision, DP=double precision)
! use the interface WriteArrReal2 
!---------------------------------------------------------------------------------------------

subroutine WriteArrReal2_SP(FileID,sizex,sizey,var,varname)
     integer(HID_T),             intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex,sizey
     real,dimension(sizex,sizey),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     INTEGER(HID_T)                         :: dspace_id,dset_id
     INTEGER(HSIZE_T),dimension(2)          :: data_dim
     integer                                :: rank,err
     
     rank=2
     data_dim(1)=sizex
     data_dim(2)=sizey
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(FileID,varname,H5T_NATIVE_REAL,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_REAL,var,data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)          
end subroutine WriteArrReal2_SP

subroutine WriteArrReal2_DP(FileID,sizex,sizey,var,varname)
     integer(HID_T),             intent(in) :: FileID  !provide either FileName or FileID   
     integer,                    intent(in) :: sizex,sizey
     real(dbpsn),dimension(sizex,sizey),intent(in) :: var  
     character(len=*),           intent(in) :: varname ! Name of saved variable in the file
     INTEGER(HID_T)                         :: dspace_id,dset_id
     INTEGER(HSIZE_T),dimension(2)          :: data_dim
     integer                                :: rank,err
     
     rank=2
     data_dim(1)=sizex
     data_dim(2)=sizey
     call h5screate_simple_f(rank,data_dim,dspace_id,err)
     call h5dcreate_f(FileID,varname,H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,data_dim,err)
     call h5dclose_f(dset_id,err)     
     call h5sclose_f(dspace_id,err)          
end subroutine WriteArrReal2_DP


!---------------------------------------------------------------------------------------------------------------------------------
!
!  Collectively (parallel) write data in HDF5 files 
!
!----------------------------------------------------------------------------------------------------------------------------------

subroutine WriteReal_Collective_SP(fid,varname,var,dspace_id,memspace,plist,offset)
	 INTEGER(HID_T)                 :: fid, dspace_id, memspace, plist
	 INTEGER(HID_T)                 :: dset_id, dspace_this
	 INTEGER(HSIZE_T), dimension(1) :: offset, data_dim1
	 real                           :: var
	 character(len=*)               :: varname
	 integer                        :: err
 
	 data_dim1(1)=1
     call h5dcreate_f(fid,varname, H5T_NATIVE_REAL, dspace_id, dset_id,err)               
     call h5dget_space_f(dset_id, dspace_this, err)
     call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset,data_dim1,err)       
     call h5dwrite_f(dset_id, H5T_NATIVE_REAL,var, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
     call h5sclose_f(dspace_this,err)
     call h5dclose_f(dset_id,err)
 
end subroutine WriteReal_Collective_SP

subroutine WriteINT_Collective(fid,varname,var,dspace_id,memspace,plist,offset)
	 INTEGER(HID_T)                 :: fid, dspace_id, memspace, plist
	 INTEGER(HID_T)                 :: dset_id, dspace_this
	 INTEGER(HSIZE_T), dimension(1) :: offset, data_dim1
	 integer                        :: var
	 character(len=*)               :: varname
	 integer                        :: err
 
	 data_dim1(1)=1
	 call h5dcreate_f(fid,varname, H5T_NATIVE_INTEGER, dspace_id, dset_id,err)               
	 call h5dget_space_f(dset_id, dspace_this, err)
	 call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset,data_dim1,err)       
	 call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,var, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
	 call h5sclose_f(dspace_this,err)
	 call h5dclose_f(dset_id,err)
 
end subroutine WriteINT_Collective


end module HDF5write