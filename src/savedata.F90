module savedata
     use parameters
     use vars
     use hdf5
	 use prtl_stats
     use savedata_routines
	 use HDF5write
     !use setup
     use comm_savedata                                      
     use fields, only : UpdateCurrentsAllEdges,SyncCurrentEdges
#ifdef gpu
     use comm_fld_gpu
	 use savedata_gpu
#endif 
     implicit none    
	 INTEGER(HID_T) :: h5real_psn
     character (len=4), dimension(19)::pvnames=(/'x   ','y   ','z   ','u   ','v   ','w   ','q   ','a   ','flv ','var1','pEx ','pEy ','pEz ',&
                                                 'pBx ','pBy ','pBz ','pJx ','pJy ','pJz '/)
     character (len=5), dimension(26)::fvnames=(/'Ex   ','Ey   ','Ez   ','Bx   ','By   ','Bz   ','Jx   ','Jy   ','Jz   ','Di   ',&
                                                 'De   ','Jxi  ','Jyi  ','Jzi  ','Jxe  ','Jye  ','Jze  ','divE ','ExVxi','EyVyi',&
                                                            'EzVzi','ExVxe','EyVye','EzVze','Eni  ','Ene  '/) 
															                                
contains
     subroutine SaveOutputCollective !all data are written in parallel by all proc
           !call SaveOutput_custom
           if((prtl_save_period.gt.0).and.(modulo(t,prtl_save_period).eq.0)) call save_particles_collective !!Note: the collective writing doesn't work !!!
           if((fld_save_period.gt.0).and.(modulo(t,fld_save_period).eq.0)) then !default: saves fld data over entire box
                call save_field_collective(fdataxi_box,fdataxf_box,fdatayi_box,fdatayf_box,fdatazi_box,fdatazf_box,fsave_ratio,0,'')
           end if
           if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call save_spec_master_all_prtl 

           if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call save_total_energy ! to save energy data
           if((modulo(t,1).eq.0).or.(t.eq.1)) call CrossCheck
           if((prtl_mean_save_period.gt.0).and.(modulo(t,prtl_mean_save_period).eq.0)) call save_prtl_mean ! mean of prtl related qunatities
           if((prtl_mean_save_period.eq.0).and.(modulo(t,spec_save_period).eq.0)) call save_prtl_mean
         
    end subroutine SaveOutputCollective
	
	subroutine InitSaveOutput
		integer :: err
		if(.not.restart) call CreateDataFolder
		call h5open_f(err)
        select case(psn)
        case(kind(1.0d0))
           call h5tcopy_f(H5T_NATIVE_DOUBLE,h5real_psn,err)
        case(kind(1.0e0))
           call h5tcopy_f(H5T_NATIVE_REAL,h5real_psn,err)
        end select
		call h5close_f(err)
	end subroutine InitSaveOutput
	
! !-----------The following subroutine are to crosscheck several things ------
    subroutine CrossCheck
          !call GetTotalNP
          !call CheckPrtlLimit
          !call CheckNegativeA
          !call CheckIonflv
          !call GetTotalIon
          !call GetNegativeQ
          !call GetPositiveQ
          !call GetTotalPositiveA
     end subroutine CrossCheck
!      subroutine CheckPrtlLimit
!           integer :: i
!           real(psn):: pxmin,pxmax,pymin,pymax
!           pxmin=1000
!           pymin=1000
!           pxmax=xmax-100
!           pymax=ymax-100
!           do i=1,prtl_arr_size
!                if(p(i)%x.lt.pxmin) pxmin=p(i)%x
!                if(p(i)%x.gt.pxmax) pxmax=p(i)%x
!                if(p(i)%y.lt.pymin) pymin=p(i)%y
!                if(p(i)%y.gt.pymax) pymax=p(i)%y
!           end do
!           !print *,'at proc',proc,'pxmin, pxmax , pymin, and pymax are : ', pxmin,pxmax,pymin,pymax
!           if(pxmin.lt.2.5) STOP 'xmin limit crossed at'
!           if(pymin.lt.2.5) STOP 'ymin limit crossed at'
!           if(pxmax.gt.xmax+0.5) STOP 'xmax limit corossed at'
!           if(pymax.gt.ymax+0.5) STOP 'ymax limit corossed at'
!      end subroutine CheckPrtlLimit
!      subroutine CheckNegativeA
!           integer :: i
!           do i=1,prtl_arr_size
!                if(p(i)%a.lt.0) then
!                     print*,'Found negative a at proc',proc,'p is',p(i)
!                     STOP
!                end if
!           end do
!      end subroutine CheckNegativeA
!      subroutine CheckIonflv
!           integer :: n
!           do n=1,prtl_arr_size
!                if((flvp(n).eq.1).and.(qp(n).ne.1)) then
!                     print*,'Found something which is flv 1 but charge not 1',qp(n),flvp(n)
!                     STOP
!                end if
!                if((qp(n).eq.0).and.(flvp(n).ne.0)) then
!                     print*,'Found something which has 0 charge but non-zero flv',qp(n),flvp(n)
!                     STOP
!                end if
!           end do
!      end subroutine CheckIonflv


!--------------------------------------------     
     
     subroutine SavePerformanceData
             if((performance_save_period.gt.0).and.(modulo(t,performance_save_period).eq.0)) call save_performance_collective
    end subroutine SavePerformanceData 
     subroutine CreateDataFolder
          if(proc.ne.0) return 
          cmd="mkdir "//trim(data_folder)
          call system(cmd)
          cmd="mkdir "//trim(data_folder)//"/restart"
          call system(cmd)
          if(performance_save_period.gt.0) then 
              cmd="mkdir "//trim(data_folder)//"/prfm"
              call system(cmd)
          end if
     end subroutine CreateDataFolder
     subroutine SaveRestartData
		 integer :: stat
          if((restart_save_period.gt.0).and.(modulo(t,restart_save_period).eq.0)) then
               call h5open_f(err)
               if(proc.eq.0) print*,'Saving the current state of the simulation ...'
               call StartTimer(41)
               call save_fields_restart
               call save_particles_restart
               call save_param_restart
			   call save_restart_time
               call h5close_f(err)
			   
			   call MPI_Barrier(MPI_COMM_WORLD, ierr) 
               
               !now delete older files, if any
        
                write(str1,'(I0)') t-restart_save_period
                write(str2,'(I0)') proc

				fname=trim(data_folder)//"/restart"//"/fld_"//trim(str2)//"_"//trim(str1)
				open(unit=1234, iostat=stat, file=fname, status='old')
				if(stat.eq.0) close(1234, status='delete')

			    fname=trim(data_folder)//"/restart"//"/prtl_"//trim(str2)//"_"//trim(str1)
				open(unit=1234, iostat=stat, file=fname, status='old')
				if(stat.eq.0) close(1234, status='delete')
				
			    fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
				open(unit=1234, iostat=stat, file=fname, status='old')
				if(stat.eq.0) close(1234, status='delete')
				
               call StopTimer(41)
               if(proc.eq.0) print*, 'Restart data saved at Time Step = ',t,'. Data write time: ',real(exec_time(41)),' sec'
          end if
     end subroutine SaveRestartData

!*************************************************************************************************************
!   Following subroutines use HDF5 to write data in serial of parallel 
!   savedata_subroutines contains several helping subroutines that are used to prepare data to be saved
!*************************************************************************************************************


!===================================================================================================
! Following subroutines write data in parallel 
!=================================================================================================== 
     
     subroutine save_particles_collective
          implicit none 
          integer:: i
          INTEGER(HSIZE_T), dimension(1) :: dmemspace_this
          
          
          call h5open_f(err) !should be deleted
          write(fname,'(I0)') t
          fname=trim(data_folder)//"/prtl_"//trim(fname)
          
          call h5open_f(err)
          prtl_arr_size_all=0
#ifdef CPU	  
          call GetSizeofCollectPrtl
#endif		  
#ifdef gpu
          call CalcPrtlGPU
#endif 		  
          
          prtl_arr_size_all(proc)=tosave_prtl_arr_size
          call ReduceToSavePrtlSize
          rank=1
          offset1(1)=0
          data_dim1(1)=0
          do i=0,nproc-1
              if(i.lt.proc) offset1(1)=offset1(1)+prtl_arr_size_all(i)
               data_dim1(1)=data_dim1(1)+prtl_arr_size_all(i)
          end do
		  if(data_dim1(1).eq.0) return
		  
          allocate(pdata_real(tosave_prtl_arr_size),pdata_int(tosave_prtl_arr_size))
          if(save_prtl_local_fld) allocate(pdata_local_field(tosave_prtl_arr_size,6))
		  
          dmemspace_this(1)=prtl_arr_size_all(proc)
          call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
          call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL , err)
          call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, err, access_prp = plist)
          call h5pclose_f(plist, err)

          call h5screate_simple_f(rank, dmemspace_this, memspace, err)
          call h5screate_simple_f(rank, data_dim1, dspace_id, err)

          call save_particles_real_arr_collective(1)          
          call save_particles_real_arr_collective(2)
          call save_particles_real_arr_collective(3)
          call save_particles_real_arr_collective(4)
          call save_particles_real_arr_collective(5)
          call save_particles_real_arr_collective(6)
          call save_particles_real_arr_collective(7)
          call save_particles_int_arr_collective(8) !integer type
          call save_particles_int_arr_collective(9) 
          call save_particles_real_arr_collective(10) 
          if(save_prtl_local_fld) then

#ifdef CPU		  
               call CalcPrtlLocalEMField
#endif			   
#ifdef gpu
               call CalcPrtlLocalFieldGPU
#endif			   
               call save_particles_real_arr_collective(11)
               call save_particles_real_arr_collective(12)
               call save_particles_real_arr_collective(13)
               call save_particles_real_arr_collective(14)
               call save_particles_real_arr_collective(15)
               call save_particles_real_arr_collective(16)
          end if
          if(save_prtl_local_curr) then
#ifdef CPU			  
			   call SyncCurrentEdges
               call CalcPrtlLocalCurr
#endif 

#ifdef gpu
			   call SyncCurrentEdgesGPU_Exclusive
			   call CalcPrtlLocalCurrGPU
#endif			   
               call save_particles_real_arr_collective(17)
               call save_particles_real_arr_collective(18)
               call save_particles_real_arr_collective(19)
          end if
          
          call h5sclose_f(dspace_id, err)
          call h5sclose_f(memspace, err)
          call h5fclose_f(fid,err)
          call h5close_f(err)
          deallocate(pdata_real,pdata_int)     
          if(save_prtl_local_fld) deallocate(pdata_local_field)                              
     end subroutine save_particles_collective
     subroutine save_particles_real_arr_collective(vid)
          integer :: vid
          INTEGER(HSIZE_T), dimension(1) :: local_data_dim
          INTEGER(HID_T) :: dspace_this
#ifdef CPU          
          call CollectPrtl(vid)
#endif
#ifdef gpu
          call CollectPrtlGPU(vid)
#endif	  
          if(vid.eq.1) pdata_real=pdata_real+xborders(procxind(proc))-3
          if(vid.eq.2) pdata_real=pdata_real+yborders(procyind(proc))-3
#ifndef twoD
          if(vid.eq.3) pdata_real=pdata_real+zborders(proczind(proc))-3
#endif     
          local_data_dim(1)=prtl_arr_size_all(proc)
          call h5dcreate_f(fid,pvnames(vid), H5T_NATIVE_REAL, dspace_id, dset_id,err)          
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,pdata_real, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace,xfer_prp = plist)  
          call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
     end subroutine save_particles_real_arr_collective
     subroutine save_particles_int_arr_collective(vid)
          integer :: vid
          INTEGER(HSIZE_T), dimension(1) :: local_data_dim
          INTEGER(HID_T) :: dspace_this
#ifdef CPU          
        call CollectPrtl(vid)
#endif	
#ifdef gpu
        call CollectPrtlGPU(vid)
#endif
          local_data_dim(1)=prtl_arr_size_all(proc)
           call h5dcreate_f(fid,pvnames(vid), H5T_NATIVE_INTEGER, dspace_id, dset_id,err)
         call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err)
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,pdata_int, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace,xfer_prp = plist)
          call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
     end subroutine save_particles_int_arr_collective
     
     
     subroutine save_field_collective(xi,xf,yi,yf,zi,zf,res,fldid,ext)
          integer :: xi,xf,yi,yf,zi,zf,res
          integer :: fldid !0 means the call to save all flds, save according to user defined true or false
          character(len=*) :: ext

          !set the parameters for saving the fld data
          fdataxi=xi
          fdataxf=xi+fsave_ratio*int((xf-xi)/fsave_ratio)
          fdatayi=yi
          fdatayf=yi+fsave_ratio*int((yf-yi)/fsave_ratio)
          fdatazi=zi
          fdatazf=zi+fsave_ratio*int((zf-zi)/fsave_ratio)
          fsave_ratio=res
          binlen=real(fsave_ratio)/2.0_psn
#ifdef gpu
          call InitSaveDataGPU(fdataxi,fdatayi,fdatazi,xborders(procxind(proc)),yborders(procyind(proc)),zborders(proczind(proc))) 
#endif
                                                   
         write(fname,'(I0)') t
          fname=trim(data_folder)//"/fld"//ext//"_"//fname
          
          call h5open_f(err)
         call GetSizeofCollectFld
          
          rank=3
          offset3(1)=int((istart-3+xborders(procxind(proc))-fdataxi)/fsave_ratio)
          offset3(2)=int((jstart-3+yborders(procyind(proc))-fdatayi)/fsave_ratio)
#ifdef twoD
        offset3(3)=0
#else           
          offset3(3)=int((kstart-3+zborders(proczind(proc))-fdatazi)/fsave_ratio)
#endif
          !print*,'At proc',proc,'offset is',offset3,fdatax
          
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
        call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL , err)
        call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, err, access_prp = plist)
        call h5pclose_f(plist, err)
          

        data_dim3(1)=fdatax
        data_dim3(2)=fdatay
        data_dim3(3)=fdataz
        call h5screate_simple_f(rank, data_dim3, memspace, err)
        data_dim3(1)=(xf-xi)/fsave_ratio+1
        data_dim3(2)=(yf-yi)/fsave_ratio+1
#ifdef twoD
        data_dim3(3)=1
#else
        data_dim3(3)=(zf-zi)/fsave_ratio+1
#endif          
        call h5screate_simple_f(rank, data_dim3, dspace_id, err)
        
		!Set the chunk-size for HDF5 files, an optimial chunk size on a specific machine can differ
! #ifdef twoD
!         chunk_dim3(1)=(xf-xi)/fsave_ratio+1!min(512,nx/fsave_ratio+1)
!         chunk_dim3(2)=(yf-yi)/fsave_ratio+1!min(2048,ny/fsave_ratio+1)
!         chunk_dim3(3)=1
! #else
!         chunk_dim3(1)=min(128,nx/fsave_ratio+1)
!         chunk_dim3(2)=min(128,ny/fsave_ratio+1)
!         chunk_dim3(3)=min(128,nz/fsave_ratio+1)
! #endif
		

          if(fldid.eq.0.or.fldid.eq.1) then      
              call save_fields_arr_collective(Ex,'Ex',1)
              call save_fields_arr_collective(Ey,'Ey',3)
              call save_fields_arr_collective(Ez,'Ez',5)
              call save_fields_arr_collective(Bx,'Bx',7)
              call save_fields_arr_collective(By,'By',9)
              call save_fields_arr_collective(Bz,'Bz',11)
          end if

         if((fldid.eq.0.and.save_tot_curr).or.fldid.eq.2) then
          !make sure that current edges are updated before taking averages
          !currents are averaged in the same way as E field, so same avgid(third argument)
                call SyncCurrentEdges  !defined in fields
                call save_fields_arr_collective(Jx,'Jx',13)
                call save_fields_arr_collective(Jy,'Jy',13)
                call save_fields_arr_collective(Jz,'Jz',13)
         end if
		

          if((fldid.eq.0.and.save_density).or.fldid.eq.3) call SaveDensityFldAll_collective

          if((fldid.eq.0.and.save_ch_flux_fld).or.fldid.eq.4) call SaveChFluxFldAll_collective

          if((fldid.eq.0.and.save_divE).or.fldid.eq.5) then
                 call CalcDivE
                 call save_fields_arr_collective(Jx,'divE',13)
          end if

          if((fldid.eq.0.and.save_EdotV_fld).or.fldid.eq.6) call SaveEdotVFldAll_collective

          if((fldid.eq.0.and.save_eng_spat).or.fldid.eq.7) call SaveEnergySpatAll_collective

          if((fldid.eq.0.and.save_velsq_fld).or.fldid.eq.8) call SaveMeanVelSquareFldAll_collective
          

          call SaveFldDataLimit !to save the limit of the data and fsave_ratio
    
          call h5sclose_f(dspace_id, err)
          call h5sclose_f(memspace, err)
          call h5fclose_f(fid,err)
          call h5close_f(err) 
          deallocate(fdata) ! fdata is allocated in GetSizeofCollectFld, at the begining of the fields data writing operation               
     end subroutine save_field_collective
     
     subroutine save_fields_arr_collective(arr,vname,avgid)
          implicit none 
          real(psn), dimension(mx,my,mz) :: arr
          character (len=*)              :: vname
          integer :: avgid
          INTEGER(HSIZE_T), dimension(3) :: local_data_dim
          INTEGER(HID_T) :: dspace_this
          
		  call CollectFld(arr,avgid)     
          local_data_dim(1)=fdatax
          local_data_dim(2)=fdatay
          local_data_dim(3)=fdataz
		  	  
!! chunked data-set		  
! 		  call h5pcreate_f(H5P_DATASET_CREATE_F, plist, err)
! 	      call h5pset_chunk_f(plist, rank, chunk_dim3, err)
!           call h5dcreate_f(fid,vname, H5T_NATIVE_REAL, dspace_id, dset_id,err,plist)
! 		  call h5pclose_f(plist, err)    
		  
          call h5dcreate_f(fid, vname, H5T_NATIVE_REAL, dspace_id, dset_id,err) ! no chunk    
          
		  call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset3,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)          
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,fdata, data_dim3, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
		  call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
          
     end subroutine save_fields_arr_collective
     
     subroutine SaveChFluxFldAll_collective
          integer :: i
          character(len=6):: vname
          do i=1,Nflvr
            if(FlvrSaveFldData(i).eq.0) cycle   
#ifdef CPU 			  
            if(FlvrType(i).eq.0) call CalcPrtlChargeFlux(i)
			if(FlvrType(i).eq.-1) call CalcTestPrtlChargeFlux(i)
#endif 
#ifdef gpu
			if(FlvrType(i).eq.0) call CalcPrtlChargeFluxGPU(i)
			if(FlvrType(i).eq.-1) call CalcTestPrtlChargeFluxGPU(i)
#endif 			
            call UpdateCurrentsAllEdges
            write(vname,'(I0)') i
            vname="Jx"//trim(vname)
            call save_fields_arr_collective(Jx,vname,13)
            write(vname,'(I0)') i
            vname="Jy"//trim(vname)
            call save_fields_arr_collective(Jy,vname,13)
            write(vname,'(I0)') i
            vname="Jz"//trim(vname)
            call save_fields_arr_collective(Jz,vname,13)
          end do
     end subroutine SaveChFluxFldAll_collective
     subroutine SaveEdotVFldAll_collective
          integer :: i
          character(len=6):: vname
          do i=1,Nflvr
              if(FlvrSaveFldData(i).eq.0) cycle     
          if(FlvrType(i).eq.0) call CalcPrtlEdotV(i)
          if(FlvrType(i).eq.-1) call CalcTestPrtlEdotV(i)
          call UpdateCurrentsAllEdges
            write(vname,'(I0)') i
            vname="ExVx"//trim(vname)
          call save_fields_arr_collective(Jx,vname,13)
            write(vname,'(I0)') i
            vname="EyVy"//trim(vname)
          call save_fields_arr_collective(Jy,vname,13)
            write(vname,'(I0)') i
            vname="EzVz"//trim(vname)
          call save_fields_arr_collective(Jz,vname,13)
         end do 
     end subroutine SaveEdotVFldAll_collective
     subroutine SaveEnergySpatAll_collective
          integer :: i
          character(len=6) :: vname
        do i=1,Nflvr  
            if(FlvrSaveFldData(i).eq.0) cycle     
            if(FlvrType(i).eq.0) call CalcPrtlEnergySpatial(i)
            if(FlvrType(i).eq.-1) call CalcTestPrtlEnergySpatial(i)
            call UpdateCurrentsAllEdges
            write(vname,'(I0)') i
            vname="En"//trim(vname)
          call save_fields_arr_collective(Jx,vname,13)
          !call save_fields_arr_collective(Jy,26,13) ! think, may be the loop can be optimised to process 3 flvr at once 
         end do 
    end subroutine SaveEnergySpatAll_collective
     
     subroutine SaveDensityFldAll_collective
          integer :: i
          character(len=6) :: vname
		  Jx=0
		  Jy=0 !To avod floating point error beuse only Jx is being used but Jy Jz are also communicated, optimization possible 
		  Jz=0
        do i=1,Nflvr  
            if(FlvrSaveFldData(i).eq.0) cycle     
#ifdef CPU			
            if(FlvrType(i).eq.0) call CalcPrtlDensity(i) ! calculate "charge" density
            if(FlvrType(i).eq.-1) call CalcTestPrtlDensity(i) 
#endif 
#ifdef gpu
            if(FlvrType(i).eq.0) call CalcPrtlDensityGPU(i) ! calculate "charge" density
            if(FlvrType(i).eq.-1) call CalcTestPrtlDensityGPU(i) 
#endif 			
			
            call UpdateCurrentsAllEdges
            write(vname,'(I0)') i
            vname="D"//trim(vname)
            call save_fields_arr_collective(Jx,vname,13) 
	  end do 
     end subroutine SaveDensityfldAll_collective
        
     subroutine SaveMeanVelSquareFldAll_collective
          integer :: i
          character(len=6):: vname
          do i=1,Nflvr
          if(FlvrSaveFldData(i).eq.0) cycle          
            if(FlvrType(i).eq.0) call CalcPrtlMeanVelSquare(i)
            if(FlvrType(i).eq.-1) call CalcTestPrtlMeanVelSquare(i)
            call UpdateCurrentsAllEdges
            write(vname,'(I0)') i
            vname="V2x"//trim(vname)
            call save_fields_arr_collective(Jx,vname,13)
            write(vname,'(I0)') i
            vname="V2y"//trim(vname)
            call save_fields_arr_collective(Jy,vname,13)
            write(vname,'(I0)') i
            vname="V2z"//trim(vname)
            call save_fields_arr_collective(Jz,vname,13)
          end do
     end subroutine SaveMeanVelSquareFldAll_collective
     
     subroutine SaveFldDataLimit
     end subroutine SaveFldDataLimit
     
 
!########### End of subroutines to save field and particle data in parallel #######################!


!################################## End of subroutines that write particle and field data in serial ##############################     
!----------------------------------------------------------------------------------------------------------------------------
!Subroutines to save average of several quantities related to prtl 
!----------------------------------------------------------------------------------------------------------------------------
subroutine save_prtl_mean
     INTEGER(HID_T) :: fid
     character (len=1024) :: FileName
     integer              :: err
     real, dimension(Nflvr):: meanKE, meanExVx,meanEyVy,meanEzVz
     real, dimension(Nflvr)::meanPx2,meanPy2,meanPz2
     real, dimension(Nflvr,Speed_spec_binlen)::SpeedBin_meanExVx,SpeedBin_meanEyVy,SpeedBin_meanEzVz
     !real, dimension(Nflvr)::meanLR

     call SumPrtlQntyMain
     call ReduceSumPrtlQntyAll
     if(proc.eq.0) then !now save the data on master proc.
          meanKE=real(sumQxKE/sumQ)
          meanExVx=real(sumQxExVx/sumQ)
          meanEyVy=real(sumQxEyVy/sumQ)
          meanEzVz=real(sumQxEzVz/sumQ)
          meanPx2=real(sumQxPx2/sumQ)
          meanPy2=real(sumQxPy2/sumQ)
          meanPz2=real(sumQxPz2/sumQ)
          
          call CalcSpeedBinMeanEV(Speed_SumQxExVx,SpeedBin_meanExVx)
          call CalcSpeedBinMeanEV(Speed_SumQxEyVy,SpeedBin_meanEyVy)
          call CalcSpeedBinMeanEV(Speed_SumQxEzVz,SpeedBin_meanEzVz)     
          
          !meanLR=real(sumQxLR/sumQ)     
               
     
          write(FileName,'(I0)') t
         FileName=trim(data_folder)//"/PrtlMean_"//trim(FileName)
         call h5open_f(err)     
          call h5fcreate_f(FileName,H5F_ACC_TRUNC_F,fid, err)
         
         call WriteArrReal1_DP(FileID=fid,sizex=Nflvr,var=sumNprtl,varname='sumNprtl')
         call WriteArrReal1_DP(FileID=fid,sizex=Nflvr,var=sumQ,varname='sumQ')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanKE,varname='meanKE')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanExVx,varname='meanExVx')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanEyVy,varname='meanEyVy')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanEzVz,varname='meanEzVz')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanPx2,varname='meanPx2')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanPy2,varname='meanPy2')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanPz2,varname='meanPz2')
         call WriteArrReal2(FileID=fid,sizex=Nflvr,sizey=Speed_spec_binlen,var=SpeedBin_meanExVx,varname='SpeedBin_meanExVx')
         call WriteArrReal2(FileID=fid,sizex=Nflvr,sizey=Speed_spec_binlen,var=SpeedBin_meanEyVy,varname='SpeedBin_meanEyVy')
         call WriteArrReal2(FileID=fid,sizex=Nflvr,sizey=Speed_spec_binlen,var=SpeedBin_meanEzVz,varname='SpeedBin_meanEzVz')
         !call WriteArrReal1(FileID=fid,sizex=Nflvr,var=meanLR,varname='meanLR') !mean Larmor radius in units of electron skin depth  
          
         call h5fclose_f(fid,err)
         call h5close_f(err)     
    
     end if
end subroutine save_prtl_mean


!---------------------------------------------------------------------------------------------------------------------------------
!     Subroutines to save spectrum of particles 
!---------------------------------------------------------------------------------------------------------------------------------

subroutine save_spec_master_all_prtl  
      if(proc.eq.0) then 
           call h5open_f(err)
           write(fname,'(I0)') t
           fname=trim(data_folder)//"/spec_"//trim(fname)
           call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
      end if
      
      if((prtl_spec_type.eq.1).or.(prtl_spec_type.eq.3)) call save_Gamma_spec_master(fid)
      if((prtl_spec_type.eq.2).or.(prtl_spec_type.eq.3)) call save_Speed_spec_master(fid)
 
      if(proc.eq.0) then 
          call h5fclose_f(fid,err)
          call h5close_f(err)     
      end if
     
end subroutine save_spec_master_all_prtl  
     
     subroutine save_Gamma_spec_master(fid)
          INTEGER(HID_T) :: fid
          call CalcGmaxGminLocal_all_prtl          
          call GetGmaxGminGlobal
          call CreateGammaSpecBin
          call CalcGammaSpectrum_all_prtl
          call ReduceGammaSpectrum
          
           if(proc.eq.0) then
              rank=3
              data_dim3(1)=Nflvr
              data_dim3(2)=Gamma_spec_binlen
              data_dim3(3)=1
              call h5screate_simple_f(rank,data_dim3,dspace_id,err)
              call h5dcreate_f(fid,'gspec',H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
              call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,spec_gamma,data_dim3,err)
              call h5dclose_f(dset_id,err)
              call h5sclose_f(dspace_id,err)
                 
              rank=1
              data_dim1(1)=Gamma_spec_binlen
              call h5screate_simple_f(rank,data_dim1,dspace_id,err)
              call h5dcreate_f(fid,'gbin',H5T_NATIVE_REAL,dspace_id,dset_id,err)
              call h5dwrite_f(dset_id,H5T_NATIVE_REAL,Gamma_spec_bin,data_dim1,err)
              call h5dclose_f(dset_id,err)
              call h5sclose_f(dspace_id,err)
         end if

     end subroutine save_Gamma_spec_master 
     
     
     subroutine save_Speed_spec_master(fid)
          INTEGER(HID_T) :: fid
          call CalcSpeedSpectrum_all_prtl
          call ReduceSpeedSpectrum
          
           if(proc.eq.0) then
              rank=3
              data_dim3(1)=Nflvr
              data_dim3(2)=Speed_spec_binlen
              data_dim3(3)=1
              call h5screate_simple_f(rank,data_dim3,dspace_id,err)
              call h5dcreate_f(fid,'vspec',H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
              call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,spec_speed,data_dim3,err)
              call h5dclose_f(dset_id,err)
              call h5sclose_f(dspace_id,err)
                 
              rank=1
              data_dim1(1)=Speed_spec_binlen
              call h5screate_simple_f(rank,data_dim1,dspace_id,err)
              call h5dcreate_f(fid,'vbin',H5T_NATIVE_REAL,dspace_id,dset_id,err)
              call h5dwrite_f(dset_id,H5T_NATIVE_REAL,Speed_spec_bin,data_dim1,err)
              call h5dclose_f(dset_id,err)
              call h5sclose_f(dspace_id,err)
         end if

     end subroutine save_Speed_spec_master 
	 
	 subroutine SaveGammaSpecInSubDomain(xi,xf,yi,yf,zi,zf,ext)
	     integer   :: xi,xf,yi,yf,zi,zf
	     character (len=1024) :: FileName
         character(len=*) :: ext
		 integer          :: err
		  
		 call CalcGmaxGminLocalInSubDomain(xi,xf,yi,yf,zi,zf)      
         call GetGmaxGminGlobal
         call CreateGammaSpecBin
         call CalcGammaSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
         call ReduceGammaSpectrum
		 if(proc.eq.0) then 
	          write(fname,'(I0)') t
	          FileName=trim(data_folder)//"/spec"//ext//"_"//fname
	          call h5open_f(err)     
	          call h5fcreate_f(FileName,H5F_ACC_TRUNC_F,fid, err)
		      call WriteArrReal1(FileID=fid,sizex=Gamma_spec_binlen,var=Gamma_spec_bin,varname='gbin')
		      call WriteArrReal2_DP(FileID=fid,sizex=Nflvr,sizey=Gamma_spec_binlen,var=spec_gamma,varname='gspec')
	          call h5fclose_f(fid,err)
	          call h5close_f(err)  
		 end if
	 end subroutine SaveGammaSpecInSubDomain

	 subroutine SaveSpeedSpecInSubDomain(xi,xf,yi,yf,zi,zf,ext)
	     integer   :: xi,xf,yi,yf,zi,zf
	     character (len=1024) :: FileName
         character(len=*) :: ext
		 integer          :: err
		  
         call CalcSpeedSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
         call ReduceSpeedSpectrum
		 if(proc.eq.0) then 
	          write(fname,'(I0)') t
	          FileName=trim(data_folder)//"/spec"//ext//"_"//fname
	          call h5open_f(err)     
	          call h5fcreate_f(FileName,H5F_ACC_TRUNC_F,fid, err)
		      call WriteArrReal1(FileID=fid,sizex=Speed_spec_binlen,var=Speed_spec_bin,varname='vbin')
		      call WriteArrReal2_DP(FileID=fid,sizex=Nflvr,sizey=Speed_spec_binlen,var=spec_speed,varname='vspec')
	          call h5fclose_f(fid,err)
	          call h5close_f(err)  
		 end if
	 end subroutine SaveSpeedSpecInSubDomain
!--------------------------------------End of spectrum subroutines------------------------------------------------------------------



     subroutine SaveParameters
       fname=trim(data_folder)//"/param"
       call h5open_f(err)
      rank=1
       data_dim1(1)=1
       call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
      call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
     
      call h5dcreate_f(fid,'c',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,c,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'compe',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,c_ompe,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'resgrid',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,fsave_ratio,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'nx',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,nx,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'ny',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,ny,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'nz',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,nz,data_dim1,err)
       call h5dclose_f(dset_id,err)

      call h5dcreate_f(fid,'qi',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,qi,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'qmi',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,qmi,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'mi',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,massi,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'me',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,masse,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'gamma',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,g0,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'epc',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,epc,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'psaveratio',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,psave_ratio,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'Nflvr',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,Nflvr,data_dim1,err)
       call h5dclose_f(dset_id,err)
        
      call h5dcreate_f(fid,'extBx',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,Bx_ext0,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'extBy',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,By_ext0,data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'extBz',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,Bz_ext0,data_dim1,err)
       call h5dclose_f(dset_id,err)
       

     call h5dcreate_f(fid,'dim',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
#ifdef twoD
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,2,data_dim1,err)
#else
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,3,data_dim1,err)
#endif
     call h5dclose_f(dset_id,err)
      
      call WriteArrReal1(FileID=fid,sizex=Nflvr,var=real(flvrqm),varname='flvrqm')

       call h5sclose_f(dspace_id,err)
       call h5fclose_f(fid,err)
       call h5close_f(err)
     end subroutine SaveParameters
	 
     subroutine SaveParam(ParamName, ParamValue)
         character(len=*) :: ParamName
		 real(psn) :: ParamValue
		 
		 fname=trim(data_folder)//"/param"
         call h5open_f(err)
         rank=1
         data_dim1(1)=1
         call h5fopen_f(fname,H5F_ACC_RDWR_F,fid, err)
         call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
     
         call h5dcreate_f(fid, ParamName, h5real_psn,dspace_id,dset_id,err)
         call h5dwrite_f(dset_id,h5real_psn, ParamValue, data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5sclose_f(dspace_id,err)
         call h5fclose_f(fid,err)
         call h5close_f(err)
	 end subroutine SaveParam
		 
	 
	 
	 
     subroutine save_total_energy
       call GatherEnergy
       if(proc.eq.0) then 
       write(fname,'(I0)') t
       fname=trim(data_folder)//"/energy_"//trim(fname)
         call h5open_f(err)
      rank=1
         data_dim1(1)=1
       call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
      call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
       
      call h5dcreate_f(fid,'KEe',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,energy(1),data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'KEi',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,energy(2),data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'Benergy',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,energy(3),data_dim1,err)
       call h5dclose_f(dset_id,err)
       
      call h5dcreate_f(fid,'Eenergy',h5real_psn,dspace_id,dset_id,err)
      call h5dwrite_f(dset_id,h5real_psn,energy(4),data_dim1,err)
       call h5dclose_f(dset_id,err)
       
       call h5sclose_f(dspace_id,err)
       call h5fclose_f(fid,err)
       call h5close_f(err)
              
       end if
     end subroutine save_total_energy
     
	 
	 
     subroutine save_performance_collective
          INTEGER(HID_T) :: dspace_this

          write(str1,'(I0)') t
          !write(str2,'(I0)') proc
          fname=trim(data_folder)//"/prfm/"//trim(str1)
          rank=1
         call h5open_f(err)
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
        call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL , err)
          call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, err, access_prp = plist)
        call h5pclose_f(plist, err)

          data_dim1(1)=nproc
        call h5screate_simple_f(rank, data_dim1, dspace_id, err)
          data_dim1(1)=1
        call h5screate_simple_f(rank, data_dim1, memspace, err)
          
          offset1(1)=proc
        
          call h5dcreate_f(fid,'tot', H5T_NATIVE_REAL, dspace_id, dset_id,err)               
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)          
          !!Warning::there seems to be some confusion about size of the data to be specified while writing the data
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,real(exec_time(31)), data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
          
          call h5dcreate_f(fid,'movdep', H5T_NATIVE_REAL, dspace_id, dset_id,err)               
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,real(exec_time(4)), data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
		
        call h5dcreate_f(fid,'movfld', H5T_NATIVE_REAL, dspace_id, dset_id,err)               
        call h5dget_space_f(dset_id, dspace_this, err)
        call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
        call h5dwrite_f(dset_id, H5T_NATIVE_REAL,real(exec_time(18)+exec_time(2)+exec_time(5)+exec_time(7)), data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
        call h5sclose_f(dspace_this,err)
        call h5dclose_f(dset_id,err)
          
          call h5dcreate_f(fid,'prtlexchange', H5T_NATIVE_REAL, dspace_id, dset_id,err)               
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,real(exec_time(12)), data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
		  
          call h5dcreate_f(fid,'fldexchange', H5T_NATIVE_REAL, dspace_id, dset_id,err)               
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,real(exec_time(18)), data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
          
          call h5dcreate_f(fid,'np', H5T_NATIVE_INTEGER, dspace_id, dset_id,err)               
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,data_dim1,err)
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,np, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
                    
          call h5sclose_f(dspace_id, err)     
          call h5pclose_f(plist, err)
          call h5sclose_f(memspace, err)
          call h5fclose_f(fid,err)
          call h5close_f(err)     
     end subroutine save_performance_collective
     
!=====================================================================================================================     
!                              #RESTART# 
!The following subroutines are used to save essential data to restart the simulation for an intermediate time step
!=====================================================================================================================     

subroutine save_fields_restart
     INTEGER(HID_T) :: h5psn
     write(str1,'(I0)') t
     write(str2,'(I0)') proc
     fname=trim(data_folder)//"/restart"//"/fld_"//trim(str2)//"_"//trim(str1)
     rank=3
     data_dim3(1)=mx
     data_dim3(2)=my
     data_dim3(3)=mz
     
     select case(psn)
     case(kind(1.0d0))
         call h5tcopy_f(H5T_NATIVE_DOUBLE,h5psn,err)
     case(kind(1.0e0))
        call h5tcopy_f(H5T_NATIVE_REAL,h5psn,err)
     end select
     
     call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
     call h5screate_simple_f(rank,data_dim3,dspace_id,err)               
     call h5dcreate_f(fid,'Ex',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,Ex,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'Ey',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,Ey,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'Ez',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,Ez,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'Bx',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,Bx,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'By',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,By,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'Bz',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,Bz,data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5sclose_f(dspace_id,err)     
     call h5fclose_f(fid,err)       
end subroutine save_fields_restart
subroutine save_particles_restart
     INTEGER(HID_T) :: h5psn
     write(str1,'(I0)') t
     write(str2,'(I0)') proc
     fname=trim(data_folder)//"/restart"//"/prtl_"//trim(str2)//"_"//trim(str1)
     rank=1
     data_dim1(1)=used_prtl_arr_size
     select case(psn)
     case(kind(1.0d0))
         call h5tcopy_f(H5T_NATIVE_DOUBLE,h5psn,err)
     case(kind(1.0e0))
        call h5tcopy_f(H5T_NATIVE_REAL,h5psn,err)
     end select

     call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)

     call h5screate_simple_f(rank,data_dim1,dspace_id,err)
     call h5dcreate_f(fid,'qp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,qp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'flvp',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,flvp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'tagp',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,tagp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
	 call h5dcreate_f(fid,'xp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,xp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'yp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,yp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'zp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,zp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'up',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,up(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'vp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,vp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'wp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,wp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'var1p',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,var1p(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     
	 call h5sclose_f(dspace_id,err)	 
     call h5fclose_f(fid,err)
end subroutine save_particles_restart

subroutine save_param_restart
     INTEGER(HID_T) :: h5psn
	 integer :: size1, nSubDomainsZThis
	 write(str1,'(I0)') t
     write(str2,'(I0)') proc
     fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
     rank=1
     data_dim1(1)=1
     call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)     

     call HDF5writeINT(fid,dspace_id,'mx',mx)
     call HDF5writeINT(fid,dspace_id,'my',my)
     call HDF5writeINT(fid,dspace_id,'mz',mz)
     call HDF5writeINT(fid,dspace_id,'prtl_arr_size',prtl_arr_size)
     call HDF5writeINT(fid,dspace_id,'test_prtl_arr_size',test_prtl_arr_size)
     call HDF5writeINT(fid,dspace_id,'used_prtl_arr_size',used_prtl_arr_size)
     call HDF5writeINT(fid,dspace_id,'used_test_prtl_arr_size',used_test_prtl_arr_size)
	
	  
     call HDF5writeINT(fid,dspace_id,'np',np)
     call HDF5writeINT(fid,dspace_id,'ntp',ntp)
     call HDF5writeINT(fid,dspace_id,'Nflvr',Nflvr)
     call HDF5writeINT(fid,dspace_id,'nSubDomainsX',nSubDomainsX)
     call HDF5writeINT(fid,dspace_id,'nSubDomainsY',nSubDomainsY)
     call HDF5writeINT(fid,dspace_id,'nSubDomainsZ',nSubDomainsZ)
     call HDF5writeINT(fid,dspace_id,'fdataxi_box',fdataxi_box)
     call HDF5writeINT(fid,dspace_id,'fdataxf_box',fdataxf_box)
     call HDF5writeINT(fid,dspace_id,'fdatayi_box',fdatayi_box)
     call HDF5writeINT(fid,dspace_id,'fdatayf_box',fdatayf_box)
     call HDF5writeINT(fid,dspace_id,'fdatazi_box',fdatazi_box)
     call HDF5writeINT(fid,dspace_id,'fdatazf_box',fdatazf_box)
	 call HDF5writeINT(fid,dspace_id,'load_balancing_type', load_balancing_type)
	 
     call HDF5writeINT(fid,dspace_id,'proc',proc)
     call HDF5writeINT(fid,dspace_id,'lproc',lproc)
     call HDF5writeINT(fid,dspace_id,'rproc',rproc)
     call HDF5writeINT(fid,dspace_id,'tproc',tproc)
     call HDF5writeINT(fid,dspace_id,'bproc',bproc)
     call HDF5writeINT(fid,dspace_id,'uproc',uproc)
     call HDF5writeINT(fid,dspace_id,'dproc',dproc)
 
	 !save BC vars
	 call HDF5writeRealDP(fid,dspace_id,'BC_Xmin_Prtl',BC_Xmin_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Xmax_Prtl',BC_Xmax_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Ymin_Prtl',BC_Ymin_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Ymax_Prtl',BC_Ymax_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Zmin_Prtl',BC_Zmin_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Zmax_Prtl',BC_Zmax_Prtl)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Xmin_Fld',BC_Xmin_Fld)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Xmax_Fld',BC_Xmax_Fld)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Ymin_Fld',BC_Ymin_Fld)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Ymax_Fld',BC_Ymax_Fld)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Zmin_Fld',BC_Zmin_Fld)
	 call HDF5writeRealDP(fid,dspace_id,'BC_Zmax_Fld',BC_Zmax_Fld)
	 
	 
	 call h5sclose_f(dspace_id,err)   
	 
	 !Flavour information for all species
     select case(psn)
     case(kind(1.0d0))
        call h5tcopy_f(H5T_NATIVE_DOUBLE,h5psn,err)
     case(kind(1.0e0))
        call h5tcopy_f(H5T_NATIVE_REAL,h5psn,err)
     end select
     data_dim1(1)=Nflvr	 
 	 call h5screate_simple_f(rank,data_dim1,dspace_id,err)
     call h5dcreate_f(fid,'flvrqm',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,flvrqm,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrCharge',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,FlvrCharge,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrSaveFldData',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,FlvrSaveFldData,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrType',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,FlvrType,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrSpare',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,FlvrSpare,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'TagCounter',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,TagCounter,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'CurrentTagID',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagID,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5sclose_f(dspace_id,err)

     rank=1
     data_dim1(1)=nSubDomainsX+1
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)
     call h5dcreate_f(fid,'xborders',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,xborders(0:nSubDomainsX),data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5sclose_f(dspace_id,err)

     rank=1
     data_dim1(1)=nSubDomainsY+1
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)
     call h5dcreate_f(fid,'yborders',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,yborders(0:nSubDomainsY),data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5sclose_f(dspace_id,err)

     rank=1
	 nSubDomainsZThis=nSubDomainsZ
#ifdef twoD
     nSubDomainsZThis=1
#endif		 
     data_dim1(1)=nSubDomainsZThis+1
 
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)
     call h5dcreate_f(fid,'zborders',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,zborders(0:nSubDomainsZThis),data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5sclose_f(dspace_id,err)

     rank=1
#ifdef twoD
     size1=nSubDomainsX*nSubDomainsY
#else
     size1=nSubDomainsX*nSubDomainsY*nSubDomainsZ
#endif
     data_dim1(1)=size1
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)
	 call h5dcreate_f(fid,'procxind',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
	 call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,procxind(0:size1-1),data_dim1,err)
     call h5dclose_f(dset_id,err) 
	 call h5dcreate_f(fid,'procyind',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
	 call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,procyind(0:size1-1),data_dim1,err)
     call h5dclose_f(dset_id,err)
	 call h5dcreate_f(fid,'proczind',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
	 call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,proczind(0:size1-1),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5sclose_f(dspace_id,err)

     rank=3
#ifdef twoD
     size1=1
#else
     size1=nSubDomainsZ
#endif
     data_dim3(1)=nSubDomainsX
     data_dim3(2)=nSubDomainsY
     data_dim3(3)=size1
     call h5screate_simple_f(rank,data_dim3,dspace_id,err)
	 call h5dcreate_f(fid,'proc_grid',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
	 call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,proc_grid(0:nSubDomainsX-1,0:nSubDomainsY-1,0:size1-1),data_dim3,err)
     call h5dclose_f(dset_id,err)
     call h5sclose_f(dspace_id,err)   
     
     call h5fclose_f(fid,err)
               
end subroutine save_param_restart

subroutine save_restart_time
	if(proc.ne.0) return
    fname=trim(data_folder)//"/restart"//"/TimeStep"
    rank=1
    data_dim1(1)=1
    call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
    call h5screate_simple_f(rank,data_dim1,dspace_id,err)     

    call HDF5writeINT(fid,dspace_id,'restart_time',t)
    call h5sclose_f(dspace_id,err)
	call h5fclose_f(fid,err)	
end subroutine save_restart_time

!============================================End of restart data save subroutines ====================================================

      
	 
!-------------------------------------------------------------------------------------------------------
! Some older subroutines 
!------------------------------------------------------------------------------------------------------	 
!      subroutine SaveOutputMaster ! data is first transferred to the root proc and then root writes the data
!           call SaveOutput_custom
!           if((prtl_save_period.gt.0).and.(modulo(t,prtl_save_period).eq.0)) call save_particles_master
!           if((fld_save_period.gt.0).and.(modulo(t,fld_save_period).eq.0))   call save_field_master
!           if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call save_spec_master_all_prtl 
!           !call save_total_energy ! to save energy data
!      end subroutine SaveOutputMaster	
!========================================================================================================
! The following subroutines are used to write field and particle data in serial: All processors send their data to the 
! master processor and then the master processor writes all the data 
! Warning : These subroutines may be borken (Parallel version is constantly updated and is the one genrally used)  
!======================================================================================================== 
!
!      subroutine save_particles_master
!           integer:: i
!           call h5open_f(err)
!           if(proc.eq.0) then
!           write(fname,'(I0)') t
!           fname=trim(data_folder)//"/prtl_"//trim(fname)
!           do i=1,nproc-1
!             call RecvPrtlSizeAtMaster(i)
!           end do
!            call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
!           call GetSizeofCollectPrtl
!           call save_particles_real_arr_master(1)
!           call save_particles_real_arr_master(2)
!           call save_particles_real_arr_master(3)
!           call save_particles_real_arr_master(4)
!           call save_particles_real_arr_master(5)
!           call save_particles_real_arr_master(6)
!           call save_particles_real_arr_master(7)
!           call save_particles_int_arr_master(8) !integer type
! #ifdef mulflvr
!           call save_particles_int_arr_master(9)
!           call save_particles_int_arr_master(10)
! #endif
!           if(save_prtl_local_fld) then
!                call CalcPrtlLocalEMField
!                call save_particles_real_arr_master(11)
!                call save_particles_real_arr_master(12)
!                call save_particles_real_arr_master(13)
!                call save_particles_real_arr_master(14)
!                call save_particles_real_arr_master(15)
!                call save_particles_real_arr_master(16)
!           end if
!
!
!
!
!            call h5fclose_f(fid,err)
!           deallocate(pdata_real,pdata_int)
!           if(save_prtl_local_fld) deallocate(pdata_local_field)
!
!           else
!                 call GetSizeofCollectPrtl
!                 call SendPrtlSizeToMaster
!                if(save_prtl_local_fld) call CalcPrtlLocalEMField
!                 call SendPrtlToMaster
!                deallocate(pdata_real,pdata_int)
!                if(save_prtl_local_fld) deallocate(pdata_local_field)
!           end if
!           call h5close_f(err)
!      end subroutine save_particles_master
!      subroutine save_particles_real_arr_master(vid)
!           integer :: vid,curr_dset_size
!           call CollectPrtl(vid)
!           if(vid.eq.1) pdata_real=pdata_real+xborders(procxind(proc))-3
!           if(vid.eq.2) pdata_real=pdata_real+yborders(procyind(proc))-3
! #ifndef twoD
!            if(vid.eq.3) pdata_real=pdata_real+zborders(proczind(proc))-3
! #endif
!
!            rank=1
!            data_dim1(1)=tosave_prtl_arr_size
!            max_dim1(1)=H5S_UNLIMITED_F
!            chunk_dim1(1)=max(tosave_prtl_arr_size,10) !to make sure that chunk_dim doesn't become 0(in case there is nothing to write at proc 0)
!           call h5screate_simple_f(rank,data_dim1,dspace_id,err,max_dim1)
!           call h5pcreate_f(H5P_DATASET_CREATE_F, plist, err)
!           call h5pset_chunk_f(plist, RANK, chunk_dim1, err)
!            call h5dcreate_f(fid,pvnames(vid), H5T_NATIVE_DOUBLE, dspace_id, dset_id,err,plist)
!           call h5dwrite_f(dset_id,H5T_NATIVE_REAL,pdata_real,data_dim1,err)
!           deallocate(pdata_real)
!             curr_dset_size=tosave_prtl_arr_size
!                  do i=1,nproc-1
!                      call RecvPrtlAtMaster(i,vid,1)
!                     if(vid.eq.1) pdata_real=pdata_real+xborders(procxind(i))-3
!                     if(vid.eq.2) pdata_real=pdata_real+yborders(procyind(i))-3
! #ifndef twoD
!                      if(vid.eq.3) pdata_real=pdata_real-3
! #endif
!                        offset1(1)=curr_dset_size
!                        dcount1(1)=prtl_arr_size_all(i)
!                        curr_dset_size=curr_dset_size+prtl_arr_size_all(i)
!                        new_size1(1)=curr_dset_size
!                        call h5dset_extent_f(dset_id, new_size1, err)
!                        data_dim1(1)=prtl_arr_size_all(i)
!                      call h5screate_simple_f(rank, data_dim1,memspace,err)
!                       call h5dget_space_f(dset_id, dspace_id, err)
!                        call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,offset1,data_dim1,err)
!                        call h5dwrite_f(dset_id, H5T_NATIVE_REAL,pdata_real,data_dim1,err,memspace,dspace_id)
!                     deallocate(pdata_real)
!                    call h5sclose_f(dspace_id,err)
!                     call h5sclose_f(memspace, err)
!                end do
!           call h5dclose_f(dset_id,err)
!            call h5pclose_f(plist,err)
!           allocate(pdata_real(tosave_prtl_arr_size))
!      end subroutine save_particles_real_arr_master
!
!      subroutine save_particles_int_arr_master(vid)
!           integer ::vid,curr_dset_size
!
!           !allocate(pdata_int(tosave_prtl_arr_size))
!           call CollectPrtl(vid)
!
!            rank=1
!            data_dim1(1)=tosave_prtl_arr_size
!            max_dim1(1)=H5S_UNLIMITED_F
!            chunk_dim1(1)=max(tosave_prtl_arr_size,10)
!           call h5screate_simple_f(rank,data_dim1,dspace_id,err,max_dim1)
!           call h5pcreate_f(H5P_DATASET_CREATE_F, plist, err)
!           call h5pset_chunk_f(plist, RANK, chunk_dim1, err)
!            call h5dcreate_f(fid,pvnames(vid), H5T_NATIVE_INTEGER, dspace_id, dset_id,err,plist)
!            call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,pdata_int,data_dim1,err)
!           deallocate(pdata_int)
!           curr_dset_size=tosave_prtl_arr_size
!                  do i=1,nproc-1
!                      call RecvPrtlAtMaster(i,vid,2)  ! 2 for integer type
!                        offset1(1)=curr_dset_size
!                        dcount1(1)=prtl_arr_size_all(i)
!                        curr_dset_size=curr_dset_size+prtl_arr_size_all(i)
!                        new_size1(1)=curr_dset_size
!                        call h5dset_extent_f(dset_id, new_size1, err)
!                        data_dim1(1)=prtl_arr_size_all(i)
!                      call h5screate_simple_f(rank, data_dim1,memspace,err)
!                       call h5dget_space_f(dset_id, dspace_id, err)
!                        call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,offset1,data_dim1,err)
!                        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,pdata_int,data_dim1,err,memspace,dspace_id)
!                     deallocate(pdata_int)
!                    call h5sclose_f(dspace_id,err)
!                     call h5sclose_f(memspace, err)
!                end do
!                call h5dclose_f(dset_id,err)
!                call h5pclose_f(plist,err)
!                allocate(pdata_int(tosave_prtl_arr_size))
!
!      end subroutine save_particles_int_arr_master
!
!
!      subroutine save_field_master
!           integer:: i
!           integer,dimension(proc,3) :: fld_size_all
!           call h5open_f(err)
!           if(proc.eq.0) then
!             write(fname,'(I0)') t
!             fname=trim(data_folder)//"/fld_"//trim(fname)
!
!             do i=1,nproc-1
!               call RecvFldSizeAtMaster(i)
!             end do
!              call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
!             call GetSizeofCollectFld
!             !fld_size_all(0,1)=fdatax !!!NOTE: fortran gives error if I define 0th element here (WHY??)
!             !fld_size_all(0,2)=fdatay
!             !fld_size_all(0,3)=fdataz
!             !fld_size_all(0,1)=fdatax
!
!             call save_fields_arr_master(Ex,1,1) ! third argument tells the averaging rule to calculate quantities at grid points
!             call save_fields_arr_master(Ey,2,3)
!             call save_fields_arr_master(Ez,3,5)
!             call save_fields_arr_master(Bx,4,7)
!             call save_fields_arr_master(By,5,9)
!             call save_fields_arr_master(Bz,6,11)
!
!             if(save_tot_curr) then
!             !make sure that current edges are updated before taking averages
!             ! currents are averaged in the same way as E field, so same avgid(third argument)
!                 call SyncCurrentEdges  !defined in fields
!                 call save_fields_arr_master(Jx,7,1)
!                 call save_fields_arr_master(Jy,8,3)
!                 call save_fields_arr_master(Jz,9,5)
!            end if
!
!               !call CalcDensity ! calculate charge density
!             !used the available Current matrices to save charge density instead.Now update Edges
!           !call UpdateCurrentsAllEdges
!
!                !!!!WARNING : for 3D XY edge update is also needed
!
!             !call save_fields_arr_master(Jx,10,13)
!             !call save_fields_arr_master(Jy,11,13)
!
!             if(save_ch_flux_fld) then
!                  call CalcChargeFlux(1)
!                  call UpdateCurrentsAllEdges
!                  call save_fields_arr_master(Jx,12,13)
!                  call save_fields_arr_master(Jy,13,13)
!                  call save_fields_arr_master(Jz,14,13)
!                  call CalcChargeFlux(-1)
!                  call UpdateCurrentsAllEdges
!                  call save_fields_arr_master(Jx,15,13)
!                  call save_fields_arr_master(Jy,16,13)
!                  call save_fields_arr_master(Jz,17,13)
!             end if
!
!             if(save_divE) then
!                  call CalcDivE
!                  call save_fields_arr_master(Jx,18,13)
!             end if
!
!
!             call h5fclose_f(fid,err)
!             deallocate(fdata)
!           else
!                call GetSizeofCollectFld
!                 call SendFldSizeToMaster
!                 call SendEMFldToMaster
!
!                if(save_tot_curr) then
!                     call SyncCurrentEdges
!                    call SendCurrToMaster
!               end if
!
!
!                !call CalcDensity
!               !call UpdateCurrentsAllEdges
!
!                !use the available Current matrices Jx,Jy,Jz to save charge density
!
!                !call SendDensityToMaster
!
!                 if(save_ch_flux_fld) then
!                    call CalcChargeFlux(1)
!                  call UpdateCurrentsAllEdges
!                  call SendVelFldToMaster_ion
!                  call CalcChargeFlux(-1)
!                  call UpdateCurrentsAllEdges
!                  call SendVelFldToMaster_elec
!                end if
!                if(save_divE) call SenddivEToMaster
!
!                deallocate(fdata)
!           end if
!           call h5close_f(err)
!      end subroutine save_field_master
!
!      subroutine save_fields_arr_master(arr,vid,avgid)
!           real(psn),dimension(mx,my,mz) :: arr
!           integer :: vid,avgid
!           rank=3
!            data_dim3(1)=fdatax
!           data_dim3(2)=fdatay
!           data_dim3(3)=fdataz
! #ifdef twoD
!         data_dim3(3)=1
! #endif
!            max_dim3(1)=H5S_UNLIMITED_F
!            max_dim3(2)=H5S_UNLIMITED_F
!            max_dim3(3)=H5S_UNLIMITED_F
!
!            chunk_dim3(1)=fdatax!mx-4
!           chunk_dim3(2)=fdatay!my
!           chunk_dim3(3)=fdataz!mz
!
!           call CollectFld(arr,avgid)
!           call h5screate_simple_f(rank,data_dim3,dspace_id,err,max_dim3)
!           call h5pcreate_f(H5P_DATASET_CREATE_F, plist, err)
!           call h5pset_chunk_f(plist, RANK, chunk_dim3, err)
!            call h5dcreate_f(fid,fvnames(vid), H5T_NATIVE_REAL, dspace_id, dset_id,err,plist)
!            call h5dwrite_f(dset_id,H5T_NATIVE_REAL,fdata,data_dim3,err)
!           deallocate(fdata)
!           offset3(1)=0
!           offset3(2)=0
!           offset3(3)=0
!           fld_size_all(0,1)=fdatax
!           fld_size_all(0,2)=fdatay
!           fld_size_all(0,2)=fdataz
!
!           new_size3=data_dim3
!                  do i=1,nproc-1
!                      call RecvFldAtMaster(i,vid)
!                        offset3(1)=offset3(1)+fld_size_all(i-1,1)
!                           if(i.lt.nSubDomainsX) new_size3(1)=new_size3(1)+fld_size_all(i,1)
!                        if(modulo(i,nSubDomainsX).eq.0) then
!                             new_size3(2)=new_size3(2)+fld_size_all(i,2)
!                             offset3(1)=0
!                             offset3(2)=offset3(2)+fld_size_all(i-1,2)
!                        end if
!                        call h5dset_extent_f(dset_id, new_size3, err)
!                           data_dim3=fld_size_all(i,:)
!                      call h5screate_simple_f(rank, data_dim3,memspace,err)
!                       call h5dget_space_f(dset_id, dspace_id, err)
!                        call h5sselect_hyperslab_f(dspace_id,H5S_SELECT_SET_F,offset3,data_dim3,err)
!                        call h5dwrite_f(dset_id, H5T_NATIVE_REAL,fdata,data_dim3,err,memspace,dspace_id)
!                     deallocate(fdata)
!                    call h5sclose_f(dspace_id,err)
!                     call h5sclose_f(memspace, err)
!                end do
!                call h5dclose_f(dset_id,err)
!                call h5pclose_f(plist,err)
!               allocate(fdata(fdatax,fdatay,fdataz)) !reallocate the array to write another variable
!      end subroutine save_fields_arr_master
!-------------------------------------------------------------------------------------------------------
!
! End of subroutines to write data in a serial way 
!
!-------------------------------------------------------------------------------------------------------
	
	 

end module savedata