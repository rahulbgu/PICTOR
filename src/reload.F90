module reload
     use parameters
     use vars
     use hdf5
	 use prtl_tag
	 INTEGER(HID_T) :: h5real_psn, h5real_dbpsn
contains 
	
     subroutine CheckRestart
		 logical :: file_exists
		 fname=trim(data_folder)//"/restart"//"/TimeStep"
		 INQUIRE(FILE=fname, EXIST=file_exists)
		 if(file_exists) restart=.true.
	 end subroutine CheckRestart	
	 
	 subroutine InitRestart
         fname=trim(data_folder)//"/restart"//"/TimeStep"
         data_dim1(1)=1
		 call h5open_f(err)
         call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
         call h5dopen_f(fid,'restart_time', dset_id, err)
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,restart_time,data_dim1,err)
		 call h5dclose_f(dset_id,err) 
		 call h5fclose_f(fid,err)
         call h5close_f(err)
		 
		 if(proc.eq.0) then
			 print*,'-------  Restarting the simulation!  --------'
			 print*,'Loading the data saved at Time Step =',restart_time
			 print*,'======================================================='
		 end if
		 
         write(str1,'(I0)') restart_time
		 write(str2,'(I0)') proc
		 
		 call h5open_f(err)
		 
		 call h5tcopy_f(H5T_NATIVE_DOUBLE,h5real_dbpsn,err)
         select case(psn)
         case(kind(1.0d0))
            call h5tcopy_f(H5T_NATIVE_DOUBLE,h5real_psn,err)
         case(kind(1.0e0))
            call h5tcopy_f(H5T_NATIVE_REAL,h5real_psn,err)
         end select
		 call h5close_f(err)
		 
		 tstart=restart_time+1
		 
	 end subroutine InitRestart
	 subroutine RestartMPI

         fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
		 call h5open_f(err)
         call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
		 call HDF5readINT(fid,'lproc',lproc)
		 call HDF5readINT(fid,'rproc',rproc)
		 call HDF5readINT(fid,'tproc',tproc)
		 call HDF5readINT(fid,'bproc',bproc)
		 call HDF5readINT(fid,'uproc',uproc)
		 call HDF5readINT(fid,'dproc',dproc)
		 
		 call h5fclose_f(fid,err)
         call h5close_f(err)
	 end subroutine RestartMPI 
     subroutine RestartInitParam

          fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
          call h5open_f(err)
		  
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
		  call HDF5readINT(fid,'mx',mx)
		  call HDF5readINT(fid,'my',my)
		  call HDF5readINT(fid,'mz',mz)
		  call HDF5readINT(fid,'prtl_arr_size',prtl_arr_size)
		  call HDF5readINT(fid,'used_prtl_arr_size',used_prtl_arr_size)
		  call HDF5readINT(fid,'test_prtl_arr_size',test_prtl_arr_size)
		  call HDF5readINT(fid,'used_test_prtl_arr_size',used_test_prtl_arr_size)
		  call HDF5readINT(fid,'np',np)
          call HDF5readINT(fid,'ntp',ntp)
		 
          call h5fclose_f(fid,err)
          call h5close_f(err)
     end subroutine RestartInitParam
     subroutine RestartInitFlds

          fname=trim(data_folder)//"/restart"//"/fld_"//trim(str2)//"_"//trim(str1)
          data_dim3(1)=mx
          data_dim3(2)=my
          data_dim3(3)=mz
		  
          call h5open_f(err)		  
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)		  
          call h5dopen_f(fid,'Ex', dset_id, err)		  
          call h5dread_f(dset_id,h5real_psn,Ex,data_dim3,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'Ey', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,Ey,data_dim3,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'Ez', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,Ez,data_dim3,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'Bx', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,Bx,data_dim3,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'By', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,By,data_dim3,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'Bz', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,Bz,data_dim3,err)
          call h5dclose_f(dset_id,err)
 		  call h5fclose_f(fid,err)
          call h5close_f(err)		  
     end subroutine RestartInitFlds
	 
     subroutine RestartInitPrtl

          call h5open_f(err)
		  
          fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
          data_dim1(1)=1
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)		  	  
          call h5dopen_f(fid,'Nflvr', dset_id, err)		  
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,Nflvr,data_dim1,err)
          call h5dclose_f(dset_id,err) 
		  
          allocate(flvrqm(Nflvr),FlvrCharge(Nflvr),FlvrSaveFldData(Nflvr),FlvrType(Nflvr),FlvrSaveRatio(Nflvr))
          allocate(CurrentTagID(Nflvr))
          data_dim1(1)=Nflvr
          call h5dopen_f(fid,'flvrqm', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,flvrqm,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrCharge', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,FlvrCharge,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrSaveFldData', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,FlvrSaveFldData,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrType', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,FlvrType,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrSaveRatio', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,FlvrSaveRatio,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'CurrentTagID', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagID,data_dim1,err)
          call h5dclose_f(dset_id,err)
		  
		  call InitPrtlTag 
          call h5dopen_f(fid,'TagBlock', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,TagBlock,data_dim1,err)
          call h5dclose_f(dset_id,err) 
		  
          data_dim1(1)=6
          call h5dopen_f(fid,'inflowBC_speed', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,inflowBC_speed,data_dim1,err)
          call h5dclose_f(dset_id,err) 
		  
		  
		  !Now Load Prtl Data 
          fname=trim(data_folder)//"/restart"//"/prtl_"//trim(str2)//"_"//trim(str1)
          data_dim1(1)=used_prtl_arr_size
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
          call h5dopen_f(fid,'qp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,qp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'flvp', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,flvp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'tagp', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,tagp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'xp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,xp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'yp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,yp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'zp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,zp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'up', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,up(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'vp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,vp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'wp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,wp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
          call h5dopen_f(fid,'var1p', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,var1p(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
 		  		  		  
		  call h5fclose_f(fid,err)
          call h5close_f(err)
     end subroutine RestartInitPrtl
	 
	 subroutine RestartAllVars
		 integer :: nSubDomainsZthis

         fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
         call h5open_f(err)
         call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)

		 nSubDomainsZthis=nSubDomainsZ	 
#ifdef twoD
		 nSubDomainsZthis=1
#endif	
   
         call h5dopen_f(fid,'xborders', dset_id, err)
		 data_dim1(1)=nSubDomainsX+1
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,xborders(0:nSubDomainsX),data_dim1,err)
         call h5dclose_f(dset_id,err)
         call h5dopen_f(fid,'yborders', dset_id, err)
		 data_dim1(1)=nSubDomainsY+1
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,yborders(0:nSubDomainsY),data_dim1,err)
         call h5dclose_f(dset_id,err)
         call h5dopen_f(fid,'zborders', dset_id, err)
		 data_dim1(1)=nSubDomainsZthis+1
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,zborders(0:nSubDomainsZthis),data_dim1,err)
         call h5dclose_f(dset_id,err)
	 
		 data_dim1(1)=nSubDomainsX*nSubDomainsY*nSubDomainsZthis
         call h5dopen_f(fid,'procxind', dset_id, err)
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,procxind(0:nSubDomainsX*nSubDomainsY*nSubDomainsZthis-1),data_dim1,err)
         call h5dclose_f(dset_id,err)
         call h5dopen_f(fid,'procyind', dset_id, err)
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,procyind(0:nSubDomainsX*nSubDomainsY*nSubDomainsZthis-1),data_dim1,err)
         call h5dclose_f(dset_id,err)
         call h5dopen_f(fid,'proczind', dset_id, err)
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,proczind(0:nSubDomainsX*nSubDomainsY*nSubDomainsZthis-1),data_dim1,err)
         call h5dclose_f(dset_id,err)
		    
		 call HDF5readINT(fid,'fdataxi_box',fdataxi_box)
		 call HDF5readINT(fid,'fdataxf_box',fdataxf_box)
		 call HDF5readINT(fid,'fdatayi_box',fdatayi_box)
		 call HDF5readINT(fid,'fdatayf_box',fdatayf_box)
		 call HDF5readINT(fid,'fdatazi_box',fdatazi_box)
		 call HDF5readINT(fid,'fdatazf_box',fdatazf_box)
		 
		 call HDF5readINT(fid,'load_balancing_type',load_balancing_type)
		 
		 call HDF5readRealDP(fid,'BC_Xmin_Prtl',BC_Xmin_Prtl)
		 call HDF5readRealDP(fid,'BC_Xmax_Prtl',BC_Xmax_Prtl)
		 call HDF5readRealDP(fid,'BC_Ymin_Prtl',BC_Ymin_Prtl)
		 call HDF5readRealDP(fid,'BC_Ymax_Prtl',BC_Ymax_Prtl)
		 call HDF5readRealDP(fid,'BC_Zmin_Prtl',BC_Zmin_Prtl)
		 call HDF5readRealDP(fid,'BC_Zmax_Prtl',BC_Zmax_Prtl)
		 call HDF5readRealDP(fid,'BC_Xmin_Fld',BC_Xmin_Fld)
		 call HDF5readRealDP(fid,'BC_Xmax_Fld',BC_Xmax_Fld)
		 call HDF5readRealDP(fid,'BC_Ymin_Fld',BC_Ymin_Fld)
		 call HDF5readRealDP(fid,'BC_Ymax_Fld',BC_Ymax_Fld)
		 call HDF5readRealDP(fid,'BC_Zmin_Fld',BC_Zmin_Fld)
		 call HDF5readRealDP(fid,'BC_Zmax_Fld',BC_Zmax_Fld)
		 
	     data_dim3(1)=nSubDomainsX
	     data_dim3(2)=nSubDomainsY
		 data_dim3(3)=nSubDomainsZthis
	 		 
         call h5dopen_f(fid,'proc_grid', dset_id, err)
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,proc_grid(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZthis-1),data_dim3,err)
         call h5dclose_f(dset_id,err)
		 call h5fclose_f(fid,err)
		 call h5close_f(err)
		 
	 end subroutine RestartAllVars
	 
	 subroutine RestartAllocatePrtlMem

         allocate(qp(prtl_arr_size))
         allocate(xp(prtl_arr_size)) 
         allocate(yp(prtl_arr_size)) 
         allocate(zp(prtl_arr_size)) 
         allocate(up(prtl_arr_size)) 
         allocate(vp(prtl_arr_size)) 
         allocate(wp(prtl_arr_size)) 
         allocate(tagp(prtl_arr_size)) 
         allocate(flvp(prtl_arr_size))
         allocate(var1p(prtl_arr_size))
		 qp=0
		 flvp=0
		 
         allocate(qtp(test_prtl_arr_size))
         allocate(xtp(test_prtl_arr_size)) 
         allocate(ytp(test_prtl_arr_size)) 
         allocate(ztp(test_prtl_arr_size)) 
         allocate(utp(test_prtl_arr_size)) 
         allocate(vtp(test_prtl_arr_size)) 
         allocate(wtp(test_prtl_arr_size)) 
         allocate(tagtp(test_prtl_arr_size)) 
         allocate(flvtp(test_prtl_arr_size))
         allocate(var1tp(test_prtl_arr_size))
		 qtp=0 
		 flvtp=0
		 
	 end subroutine RestartAllocatePrtlMem
     
	 subroutine HDF5readINT(fid,varname,var)
	      INTEGER(HID_T)                 :: fid,dset_id 
	      INTEGER(HSIZE_T), dimension(1) :: data_dim1
	      character(len=*)               :: varname
	      integer                        :: var,err
	      data_dim1(1)=1
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,H5T_NATIVE_INTEGER,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readINT

	 subroutine HDF5readReal(fid,varname,var)
	      INTEGER(HID_T)                 :: fid,dset_id 
	      INTEGER(HSIZE_T), dimension(1) :: data_dim1
	      character(len=*)               :: varname
	      integer                        :: var,err
	      data_dim1(1)=1
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,h5real_psn,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readReal
	 
	 subroutine HDF5readRealDP(fid,varname,var)
	      INTEGER(HID_T)                 :: fid,dset_id 
	      INTEGER(HSIZE_T), dimension(1) :: data_dim1
	      character(len=*)               :: varname
		  real(dbpsn)                    :: var
	      integer                        :: err
	      data_dim1(1)=1
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,h5real_dbpsn,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readRealDP
end module reload