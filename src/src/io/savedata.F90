module savedata
     use parameters
     use vars
     use hdf5
	 use prtl_stats
     use savedata_routines
	 use HDF5write
	 use comm_fld  
	 use movdep_routines                                  
#ifdef gpu
     use comm_fld_gpu
	 use savedata_gpu
#endif 
     implicit none    
	 INTEGER(HID_T) :: h5real_psn	 
	 real(dbpsn), dimension(6) :: spec_lims	, psave_lims
															
 	 INTERFACE SaveParam
 		 module procedure SaveParam_SP, SaveParam_DP, SaveParam_INT
 	 END INTERFACE 
													                                
contains
     subroutine SaveOutputCollective !all data are written in parallel by all proc

           if((prtl_save_period.gt.0).and.(modulo(t,prtl_save_period).eq.0)) then 
			call save_particles( psave_lims ) 
		 end if 
		   
           if((fld_save_period.gt.0).and.(modulo(t,fld_save_period).eq.0)) then !default: saves fld data over entire box
                call save_field_collective(box_bounds(1),box_bounds(2),box_bounds(3),box_bounds(4),box_bounds(5),box_bounds(6),fsave_ratio,'')
           end if
           
		 if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call save_spec (spec_lims,'')

           if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call save_total_energy ! to save energy data
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
		
		call SetSaveLimitsFullDomain(psave_lims)
		call SetSaveLimitsFullDomain(spec_lims)
	end subroutine InitSaveOutput
	     
     subroutine SavePerformanceData
             if((performance_save_period.gt.0).and.(modulo(t,performance_save_period).eq.0)) call save_performance_data
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
		 integer     :: stat
           integer     :: save_now, exit_now
           real(dbpsn) :: time_since_start 

           save_now = 0
           exit_now = 0
           
           if(saveAndExitAfterHours.gt.0) then 
               !to avoid frequent global communication, check should not occur in every time step 
               !currently "fld_save_period" is selected to check periodically the condition to save and exit
               if((fld_save_period.gt.0).and.(modulo(t,fld_save_period).eq.0)) then
                    if(proc.eq.0) then
                         !save the restart data only once right at the end
                         time_since_start = (get_current_time() - exec_time(1))/3600.0_dbpsn ! in hours
                         if(time_since_start.gt.saveAndExitAfterHours) then
                              save_now = 1
                              exit_now = 1
                         end if
                    end if
                    call MPI_Bcast( save_now,1, MPI_INTEGER, 0, MPI_COMM_WORLD)
                    call MPI_Bcast( exit_now,1, MPI_INTEGER, 0, MPI_COMM_WORLD)
               end if 
           else 
               !save the restart data periodically
               if((restart_save_period.gt.0).and.(modulo(t,restart_save_period).eq.0)) then
                    save_now =1
               end if
          end if
          
          if(save_now .eq. 1) then
               call h5open_f(err)
               if(proc.eq.0) print*,'Saving the current state of the simulation ...'
               call StartTimer(41)
               call save_fields_restart
               call save_particles_restart
               call save_param_restart
               call save_restart_time
               call h5close_f(err)
               
               call MPI_Barrier(MPI_COMM_WORLD, ierr) 
               
               !now delete older restart files, if present
     
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

          if( exit_now .eq. 1) then
               if(proc.eq.0) print*,'Exiting . Allocated execution time is exhausted.'
               call MPI_Barrier(MPI_COMM_WORLD, ierr)  
               call Abort(10)
          end if 
     end subroutine SaveRestartData

!===================================================================================================
! Set default limits of the domain for output
!=================================================================================================== 
	subroutine SetSaveLimitsFullDomain(lims)
		real(dbpsn), dimension(6) :: lims
		
		lims(1) = real(-huge(1),dbpsn)
		lims(2) = real(huge(1),dbpsn)
		lims(3) = real(-huge(1),dbpsn)
		lims(4) = real(huge(1),dbpsn)
		lims(5) = real(-huge(1),dbpsn)
		lims(6) = real(huge(1),dbpsn)
#ifdef twoD
          lims(5) = -1.0_dbpsn 
		lims(6) = 2.0_dbpsn
#endif	
	end subroutine SetSaveLimitsFullDomain 

!===================================================================================================
! write particle data in parallel 
!=================================================================================================== 
     
     subroutine save_particles(save_lims)
          real(dbpsn), dimension(6) :: save_lims
          INTEGER(HSIZE_T), dimension(1) :: dmemspace_this
          real, dimension(:), allocatable :: x,y,z,u,v,w,q,var1
		real, dimension(:), allocatable :: pFx,pFy,pFz 
		real(psn), dimension(:,:), allocatable  :: vecJx_this, vecJy_this, vecJz_this
		integer, dimension(:), allocatable :: flv, tag, proc_id
          real, dimension(:), allocatable :: qm, wt
          integer, dimension(:), allocatable :: save_flag ! flag particles (flag =1) that are to be saved 
		integer:: i
		  
          
		call SetSaveLimsLocal(psave_lims_local,save_lims)
		  
          call h5open_f(err) 

          write(fname,'(I0)') t
          fname=trim(data_folder)//"/prtl_"//trim(fname)

#ifdef gpu
          call CopyTaggedPrtlFromGPU
#endif
          
		  
		if( .not. allocated(prtl_arr_size_all)) allocate(prtl_arr_size_all(nproc))
          prtl_arr_size_all=0

          allocate(save_flag(used_prtl_arr_size))

#ifdef CPU	  
          call GetSizeofSavePrtl(save_flag)
#endif		  
! #ifdef gpu
!           call CalcPrtlGPU
! #endif	  
          
          prtl_arr_size_all(proc+1)=nprtl_save_this
          
		call MPI_ALLGATHER(nprtl_save_this,1,MPI_INTEGER,prtl_arr_size_all,1,MPI_INTEGER,MPI_COMM_WORLD)
          
		rank=1
          offset1(1)=0
          data_dim1(1)=0
          do i=1,nproc
               if(i.lt.proc+1) offset1(1)=offset1(1)+prtl_arr_size_all(i)
               data_dim1(1)=data_dim1(1)+prtl_arr_size_all(i)
          end do
		if(data_dim1(1).eq.0) return
		  
		allocate(q(nprtl_save_this),x(nprtl_save_this),y(nprtl_save_this),z(nprtl_save_this),u(nprtl_save_this),v(nprtl_save_this),w(nprtl_save_this),var1(nprtl_save_this),flv(nprtl_save_this))
          allocate(tag(nprtl_save_this),proc_id(nprtl_save_this))
          allocate(qm(nprtl_save_this),wt(nprtl_save_this))

            
		  
          dmemspace_this(1)=prtl_arr_size_all(proc+1)
          call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
          call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD%MPI_VAL, MPI_INFO_NULL%MPI_VAL , err)
          call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, err, access_prp = plist)
          call h5pclose_f(plist, err)

          call h5screate_simple_f(rank, dmemspace_this, memspace, err)
          call h5screate_simple_f(rank, data_dim1, dspace_id, err)

		  !extract complete attributes of the tagged particles in x,y,... arrays 
		  call GatherSavePrtl(nprtl_save_this,save_flag,q,x,y,z,u,v,w,var1,flv,tag,proc_id,qm,wt) 
          
		  call save_particles_real_arr_collective(nprtl_save_this,q,'q')
		  
		  !x, y, z are saved in the global cordinated, local q is resued 
		  q = x + xborders(procxind)-3          
          call save_particles_real_arr_collective(nprtl_save_this,q,'x')
		  q = y + yborders(procyind)-3
          call save_particles_real_arr_collective(nprtl_save_this,q,'y')
#ifndef twoD
		  q = z + zborders(proczind)-3
		  call save_particles_real_arr_collective(nprtl_save_this,q,'z')
          
#else 
          call save_particles_real_arr_collective(nprtl_save_this,z,'z')
#endif		  

          
		call save_particles_real_arr_collective(nprtl_save_this,u,'u')
          call save_particles_real_arr_collective(nprtl_save_this,v,'v')
          call save_particles_real_arr_collective(nprtl_save_this,w,'w')
          call save_particles_real_arr_collective(nprtl_save_this,var1,'var1') 
          call save_particles_int_arr_collective(nprtl_save_this,flv,'flv') 

          call save_particles_int_arr_collective(nprtl_save_this,tag,'tag') 
          call save_particles_int_arr_collective(nprtl_save_this,proc_id,'proc') 

          call save_particles_real_arr_collective(nprtl_save_this,qm,'qm') 
          call save_particles_real_arr_collective(nprtl_save_this,wt,'wt') 
		  
		deallocate(q,u,v,w,var1,flv)
          deallocate(tag,proc_id)
          deallocate(qm,wt)

		allocate(pFx(nprtl_save_this), pFy(nprtl_save_this), pFz(nprtl_save_this))
          
		  if(save_prtl_local_fld) then
			  
			   call GatherPrtlLocalField(nprtl_save_this,x,y,z,pFx,pFy,pFz,vecEx,vecEy,vecEz)
		   
               call save_particles_real_arr_collective(nprtl_save_this,pFx,'pEx')
               call save_particles_real_arr_collective(nprtl_save_this,pFy,'pEy')
               call save_particles_real_arr_collective(nprtl_save_this,pFz,'pEz')
			   
			   call GatherPrtlLocalField(nprtl_save_this,x,y,z,pFx,pFy,pFz,vecBx,vecBy,vecBz)
			   
               call save_particles_real_arr_collective(nprtl_save_this,pFx,'pBx')
               call save_particles_real_arr_collective(nprtl_save_this,pFy,'pBy')
               call save_particles_real_arr_collective(nprtl_save_this,pFz,'pBz')
          end if
          if(save_prtl_local_curr) then

			   ! make sure that the values in the ghost cells are updated and then compute VecJ_this 
#ifdef twoD
			   allocate(VecJx_this(4,mx*my), VecJy_this(4,mx*my), VecJz_this(4,mx*my))
#else			   
			   allocate(VecJx_this(8,mx*my*mz), VecJy_this(8,mx*my*mz), VecJz_this(8,mx*my*mz))
#endif			   
               VecJx_this = 0; VecJy_this = 0; VecJz_this = 0; 
			   
			   call SendRecvFlds(Jx,Jy,Jz)
			   call GatherVecEfld(VecJx_this,VecJy_this,VecJz_this,Jx,Jy,Jz)

               call GatherPrtlLocalField(nprtl_save_this,x,y,z,pFx,pFy,pFz,vecJx_this,vecJy_this,vecJz_this)
			   
			   call save_particles_real_arr_collective(nprtl_save_this,pFx,'pJx')
			   call save_particles_real_arr_collective(nprtl_save_this,pFy,'pJy')
			   call save_particles_real_arr_collective(nprtl_save_this,pFz,'pJz')
			   
			   deallocate(VecJx_this,VecJy_this,VecJz_this)
          end if
          
          call h5sclose_f(dspace_id, err)
          call h5sclose_f(memspace, err)
          call h5fclose_f(fid,err)
          call h5close_f(err)
		  
		deallocate(pFx,pFy,pFz)
		deallocate(x,y,z)
          deallocate(save_flag)
          deallocate(prtl_arr_size_all)
		  
		call MPI_Barrier(MPI_COMM_WORLD)
          if(proc.eq.0) call SaveLimits(fname,psave_lims)
                             
     end subroutine save_particles
	 
     subroutine save_particles_real_arr_collective(size,var,vname)
		integer                        :: size
		real, dimension(size)          :: var
		character (len=*)              :: vname
          integer :: vid
          INTEGER(HSIZE_T), dimension(1) :: local_data_dim
          INTEGER(HID_T) :: dspace_this
   
          local_data_dim(1)=prtl_arr_size_all(proc+1)
          call h5dcreate_f(fid,vname, H5T_NATIVE_REAL, dspace_id, dset_id,err)          
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,var, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace,xfer_prp = plist)  
          call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
     end subroutine save_particles_real_arr_collective
     subroutine save_particles_int_arr_collective(size,var,vname)
	      integer                           :: size
	      integer, dimension(size)          :: var
	      character (len=*)                 :: vname
          INTEGER(HSIZE_T), dimension(1)    :: local_data_dim
          INTEGER(HID_T) :: dspace_this

          local_data_dim(1)=prtl_arr_size_all(proc+1)
          call h5dcreate_f(fid,vname, H5T_NATIVE_INTEGER, dspace_id, dset_id,err)
          call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset1,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err)
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)
          call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER,var, data_dim1, err, file_space_id = dspace_this, mem_space_id = memspace,xfer_prp = plist)
          call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
     end subroutine save_particles_int_arr_collective
     
     
     subroutine save_field_collective(xi,xf,yi,yf,zi,zf,res,ext)
          integer :: xi,xf,yi,yf,zi,zf,res
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
          call InitSaveDataGPU(fdataxi,fdatayi,fdatazi,xborders(procxind),yborders(procyind),zborders(proczind)) 
#endif
                                                   
          write(fname,'(I0)') t
          fname=trim(data_folder)//"/fld"//ext//"_"//fname
          
          call h5open_f(err)
          call GetSizeofCollectFld
          
          rank=3
          offset3(1)=int((istart-3+xborders(procxind)-fdataxi)/fsave_ratio)
          offset3(2)=int((jstart-3+yborders(procyind)-fdatayi)/fsave_ratio)
#ifdef twoD
          offset3(3)=0
#else           
          offset3(3)=int((kstart-3+zborders(proczind)-fdatazi)/fsave_ratio)
#endif
          !print*,'At proc',proc,'offset is',offset3,fdatax
          
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
        call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD%MPI_VAL, MPI_INFO_NULL%MPI_VAL , err)
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
        
		!Set the chunk-size for the HDF5 dataset
#ifdef twoD
        chunk_dim3(1)=min(256,(xf-xi)/fsave_ratio+1)
        chunk_dim3(2)=min(256,(yf-yi)/fsave_ratio+1)
        chunk_dim3(3)=1
#else
        chunk_dim3(1)=min(64,(xf-xi)/fsave_ratio+1)
        chunk_dim3(2)=min(64,(yf-yi)/fsave_ratio+1)
        chunk_dim3(3)=min(64,(zf-zi)/fsave_ratio+1)
#endif
		

                
          call save_fields_arr_collective(Ex,'Ex',1)
          call save_fields_arr_collective(Ey,'Ey',3)
          call save_fields_arr_collective(Ez,'Ez',5)
          call save_fields_arr_collective(Bx,'Bx',7)
          call save_fields_arr_collective(By,'By',9)
          call save_fields_arr_collective(Bz,'Bz',11)
			  

          if(save_tot_curr) then
          !make sure that current edges are updated before taking averages
          !currents are averaged in the same way as E field, so same avgid(third argument)
                call SendRecvFlds(Jx,Jy,Jz)  !defined in fields
                call save_fields_arr_collective(Jx,'Jx',13)
                call save_fields_arr_collective(Jy,'Jy',13)
                call save_fields_arr_collective(Jz,'Jz',13)
          end if
		
          if(save_divE) then
                 call CalcDivE(Jx,Ex,Ey,Ez)
                 call save_fields_arr_collective(Jx,'divE',13)
          end if

          if(save_density) call SaveDenAllFlv(1,'D')

          if(save_ch_flux_fld) call SaveFluxAllFlv(1,'Vx','Vy','Vz')
		  
		if(save_velsq_fld) call SaveFluxAllFlv(2,'V2x','V2y','V2z')

          if(save_EdotV_fld) call SaveFluxAllFlv(3,'ExVx','EyVy','EzVz')

          if(save_eng_spat) call SaveDenAllFlv(2,'Eng')

         
          
          call h5sclose_f(dspace_id, err)
          call h5sclose_f(memspace, err)
		  
          call h5fclose_f(fid,err)
          call h5close_f(err) 
          deallocate(fdata) ! fdata is allocated in GetSizeofCollectFld, at the begining of the fields data writing operation 
		  
          call MPI_Barrier(MPI_COMM_WORLD)
          if(proc.eq.0) call SaveFldDataParam(fname,xi,xf,yi,yf,zi,zf,fsave_ratio)              
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
		  call h5pcreate_f(H5P_DATASET_CREATE_F, plist, err)
	      call h5pset_chunk_f(plist, rank, chunk_dim3, err)
          call h5dcreate_f(fid,vname, H5T_NATIVE_REAL, dspace_id, dset_id,err,plist)
		  call h5pclose_f(plist, err)
		  
          !call h5dcreate_f(fid, vname, H5T_NATIVE_REAL, dspace_id, dset_id,err) ! no chunking    
          
		  call h5dget_space_f(dset_id, dspace_this, err)
          call h5sselect_hyperslab_f(dspace_this,H5S_SELECT_SET_F,offset3,local_data_dim,err)
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)          
          call h5dwrite_f(dset_id, H5T_NATIVE_REAL,fdata, data_dim3, err, file_space_id = dspace_this, mem_space_id = memspace, xfer_prp = plist)     
		call h5pclose_f(plist, err)
          call h5sclose_f(dspace_this,err)
          call h5dclose_f(dset_id,err)
          
     end subroutine save_fields_arr_collective
	 
     subroutine SaveFluxAllFlv(qnty_id,vname_x,vname_y,vname_z)
		  integer :: qnty_id
          character(len=*):: vname_x, vname_y, vname_z
		  character(len=6):: vname
		  integer :: i
		  
          
		  do i=1,Nflvr
          	if(FlvrSaveFldData(i).eq.0) cycle   
            
			call CalcPrtlFlux(i,qnty_id,Jx,Jy,Jz)
			
			call UpdateCurrentAllFaces(Jx,Jy,Jz)
			call SendRecvFlds(Jx,Jy,Jz)
			
            write(vname,'(I0)') i
            vname=trim(vname_x)//trim(vname)
            call save_fields_arr_collective(Jx,vname,13)
            write(vname,'(I0)') i
            vname=trim(vname_y)//trim(vname)
            call save_fields_arr_collective(Jy,vname,13)
            write(vname,'(I0)') i
            vname=trim(vname_z)//trim(vname)
            call save_fields_arr_collective(Jz,vname,13)
          end do
     end subroutine SaveFluxAllFlv
	 
     subroutine SaveDenAllFlv(qnty_id,vname)
		  integer :: qnty_id
          character(len=*):: vname
		  character(len=6):: vname_this
		  integer :: i
		  
          
		  do i=1,Nflvr
          	if(FlvrSaveFldData(i).eq.0) cycle   
            
			call CalcPrtlDen(i,qnty_id,Jx)
		            
			call UpdateCurrentAllFaces(Jx,Jy,Jz)
			call SendRecvFlds(Jx,Jy,Jz)
			
            write(vname_this,'(I0)') i
            vname_this=trim(vname)//trim(vname_this)
            call save_fields_arr_collective(Jx,vname_this,13)

          end do
     end subroutine SaveDenAllFlv
	 
	 !save the x,y,z range of box and the downsampling ratio res = fsave_ratio
	 subroutine SaveFldDataParam(FileName,xi,xf,yi,yf,zi,zf,res)
		  character (len=*) :: FileName
		  integer :: xi,xf,yi,yf,zi,zf,res
		  
		  call h5open_f(err)
		  call h5fopen_f(FileName,H5F_ACC_RDWR_F,fid, err)
		  
          
          rank=1
          data_dim1(1)=1
		  call h5screate_simple_f(rank,data_dim1,dspace_id,err)
		  
		  call WriteINT(fid,dspace_id,'xi',xi)
		  call WriteINT(fid,dspace_id,'xf',xf)
		  call WriteINT(fid,dspace_id,'yi',yi)
		  call WriteINT(fid,dspace_id,'yf',yf)
		  call WriteINT(fid,dspace_id,'zi',zi)
		  call WriteINT(fid,dspace_id,'zf',zf)
		  call WriteINT(fid,dspace_id,'res',res)
		  
		  call h5sclose_f(dspace_id, err)
		  call h5fclose_f(fid,err)
	       call h5close_f(err)
	 end subroutine SaveFldDataParam
	 
	 !include limits of the spatial doamin in the outfile file
	 subroutine SaveLimits(FileName,lims)
		  character (len=*) :: FileName
		  real(dbpsn), dimension(6) :: lims		  
		  
		  call h5open_f(err)
		  call h5fopen_f(FileName,H5F_ACC_RDWR_F,fid, err)
		  
          rank=1
          data_dim1(1)=1
		  call h5screate_simple_f(rank,data_dim1,dspace_id,err)
		  
		  call WriteReal(fid,dspace_id,'xi',lims(1))
		  call WriteReal(fid,dspace_id,'xf',lims(2))
		  call WriteReal(fid,dspace_id,'yi',lims(3))
		  call WriteReal(fid,dspace_id,'yf',lims(4))
		  call WriteReal(fid,dspace_id,'zi',lims(5))
		  call WriteReal(fid,dspace_id,'zf',lims(6))
		  
		  call h5sclose_f(dspace_id, err)
		  call h5fclose_f(fid,err)
	      call h5close_f(err)
	 end subroutine SaveLimits
     

     
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
     !real, dimension(Nflvr) :: meanLR

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
          
                    
         write(FileName,'(I0)') t
         FileName=trim(data_folder)//"/PrtlMean_"//trim(FileName)
         call h5open_f(err)     
         call h5fcreate_f(FileName,H5F_ACC_TRUNC_F,fid, err)
         
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=sumNprtl,varname='sumNprtl')
         call WriteArrReal1(FileID=fid,sizex=Nflvr,var=sumQ,varname='sumQ')
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
          
         call h5fclose_f(fid,err)
         call h5close_f(err)     
    
     end if
end subroutine save_prtl_mean


!---------------------------------------------------------------------------------------------------------------------------------
!     Subroutines to save spectrum of particles 
!---------------------------------------------------------------------------------------------------------------------------------

     subroutine save_spec(lims,ext)
		  real(dbpsn), dimension(6) :: lims	
		  character(len=*) :: ext
 	      character (len=1024) :: FileName
		  
          write(fname,'(I0)') t
          FileName=trim(data_folder)//"/spec"//ext//"_"//fname
	 
	 	  call SaveGammaSpecInSubDomain(lims(1),lims(2),lims(3),lims(4),lims(5),lims(6),FileName)
	 	  call SaveSpeedSpecInSubDomain(lims(1),lims(2),lims(3),lims(4),lims(5),lims(6),FileName)
		  
		  if(proc.eq.0) call SaveLimits(FileName,lims)
     
	 end subroutine save_spec 
     
	 subroutine SaveGammaSpecInSubDomain(xi,xf,yi,yf,zi,zf,FileName)
	     real(dbpsn)   :: xi,xf,yi,yf,zi,zf
          character(len=*) :: FileName
		integer          :: err
		  
		 call CalcGmaxLocalInSubDomain(xi,xf,yi,yf,zi,zf)      
           call GetGmaxGlobal
           call InitGammaSpec
           call CalcGammaSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
           call AllReduceGammaSpectrum
		 if(proc.eq.0) then 
	          call h5open_f(err)     
	          call h5fcreate_f(FileName,H5F_ACC_TRUNC_F,fid, err)
		     call WriteArrReal1(FileID=fid,sizex=Gamma_spec_binlen,var=Gamma_spec_bin,varname='gbin')
		     call WriteArrReal2(FileID=fid,sizex=Nflvr,sizey=Gamma_spec_binlen,var=spec_gamma,varname='gspec')
	          call h5fclose_f(fid,err)
	          call h5close_f(err)  
		 end if
           deallocate( spec_gamma, Gamma_spec_bin )
	 end subroutine SaveGammaSpecInSubDomain

	 subroutine SaveSpeedSpecInSubDomain(xi,xf,yi,yf,zi,zf,FileName)
	     real(dbpsn)   :: xi,xf,yi,yf,zi,zf
          character(len=*) :: FileName
		integer          :: err
		  
          call CalcSpeedSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
          call AllReduceSpeedSpectrum
		if(proc.eq.0) then 
	          call h5open_f(err)     
	          call h5fopen_f(FileName,H5F_ACC_RDWR_F,fid, err)
		     call WriteArrReal1(FileID=fid,sizex=Speed_spec_binlen,var=Speed_spec_bin,varname='vbin')
		     call WriteArrReal2(FileID=fid,sizex=Nflvr,sizey=Speed_spec_binlen,var=spec_speed,varname='vspec')
	          call h5fclose_f(fid,err)
	          call h5close_f(err)  
		end if
	 end subroutine SaveSpeedSpecInSubDomain
!--------------------------------------End of spectrum subroutines------------------------------------------------------------------



     subroutine SaveParameters
		 integer :: sim_dim, cord_cyl, i
           integer, dimension(:), allocatable :: int_arr

		 fname=trim(data_folder)//"/param"
         call h5open_f(err)
         rank=1
         data_dim1(1)=1
         call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
         call h5screate_simple_f(rank,data_dim1,dspace_id,err)   
		 
		 call WriteReal(fid,dspace_id,'c',c)
		 call WriteReal(fid,dspace_id,'c_ompe',c_ompe)  
		 call WriteINT(fid,dspace_id,'resgrid',fsave_ratio)  
		 call WriteINT(fid,dspace_id,'nx',nx)  
		 call WriteINT(fid,dspace_id,'ny',ny)  
		 call WriteINT(fid,dspace_id,'nz',nz)  
		 call WriteReal(fid,dspace_id,'qi',qi)
		 call WriteReal(fid,dspace_id,'qmi',qmi)
		 call WriteReal(fid,dspace_id,'massi',massi)
		 call WriteReal(fid,dspace_id,'masse',masse)
           call WriteReal(fid,dspace_id,'mi_me',mi_me)
           call WriteReal(fid,dspace_id,'gamma0',g0)
           call WriteINT(fid,dspace_id,'epc',epc)
	      call WriteINT(fid,dspace_id,'psave_ratio',psave_ratio)
	      call WriteINT(fid,dspace_id,'Nflvr',Nflvr)
           call WriteReal(fid,dspace_id,'extBx',Bx_ext0)
           call WriteReal(fid,dspace_id,'extBy',By_ext0)
	      call WriteReal(fid,dspace_id,'extBz',Bz_ext0)
           call WriteReal(fid,dspace_id,'ne0',ne0)

           call WriteReal(fid,dspace_id,'grid_dx',grid_dx)
           call WriteReal(fid,dspace_id,'grid_dy',grid_dy)
           call WriteReal(fid,dspace_id,'grid_dz',grid_dz)
           call WriteReal(fid,dspace_id,'cell_volume',cell_volume)
           call WriteINT(fid,dspace_id,'epc_x',epc_x)
           call WriteINT(fid,dspace_id,'epc_y',epc_y)
           call WriteINT(fid,dspace_id,'epc_z',epc_z)

           call WriteReal(fid,dspace_id,'skin_depth_cm',skin_depth_cm)
           call WriteReal(fid,dspace_id,'elc_per_prtl',elc_per_prtl)
           call WriteReal(fid,dspace_id,'grid_spacing_cm',grid_spacing_cm)
           call WriteReal(fid,dspace_id,'E_SI',E_SI)




       
         sim_dim = 3
#ifdef twoD
         sim_dim = 2 
#endif		 
         call WriteINT(fid,dspace_id,'sim_dim',sim_dim)
		 
		 cord_cyl = 0
#ifdef cyl
         cord_cyl = 1 
#endif		 
         call WriteINT(fid,dspace_id,'cord_cyl',cord_cyl)
         
		 call h5sclose_f(dspace_id,err)
		 
         
		call WriteArrReal1(FileID=fid,sizex=Nflvr,var=real(flvrqm),varname='flvrqm')
          call WriteArrReal1(FileID=fid,sizex=Nflvr,var=real(FlvrCharge),varname='FlvrCharge')
 	     call WriteINT1(fid,Nflvr,'FlvrSaveFldData',FlvrSaveFldData)
 	     call WriteINT1(fid,Nflvr,'FlvrType',FlvrType)
 	     call WriteINT1(fid,Nflvr,'FlvrSaveRatio',FlvrSaveRatio)

          allocate(int_arr(Nflvr))
          int_arr = [(flvr_prpt(i)%Z_nucleus , i=1, Nflvr)]
          call WriteINT1(fid,Nflvr,'ionization_z_nucleus',int_arr)
          deallocate(int_arr)
	  
         call h5fclose_f(fid,err)
         call h5close_f(err)
     end subroutine SaveParameters
	 

!---------------------------------------------------------------------------------------------------------------------------------
!   Append the 'param' file and write setup specific paramters
!---------------------------------------------------------------------------------------------------------------------------------	 
	 
     subroutine SaveParam_SP(ParamName, ParamValue)
         character(len=*) :: ParamName
		 real :: ParamValue
		 
		 fname=trim(data_folder)//"/param"
         call h5open_f(err)
         rank=1
         data_dim1(1)=1
         call h5fopen_f(fname,H5F_ACC_RDWR_F,fid, err)
         call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
     
         call h5dcreate_f(fid, ParamName, H5T_NATIVE_REAL,dspace_id,dset_id,err)
         call h5dwrite_f(dset_id,H5T_NATIVE_REAL, ParamValue, data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5sclose_f(dspace_id,err)
         call h5fclose_f(fid,err)
         call h5close_f(err)
	 end subroutine SaveParam_SP
	 
     subroutine SaveParam_DP(ParamName, ParamValue)
         character(len=*) :: ParamName
		 real(dbpsn) :: ParamValue
		 
		 fname=trim(data_folder)//"/param"
         call h5open_f(err)
         rank=1
         data_dim1(1)=1
         call h5fopen_f(fname,H5F_ACC_RDWR_F,fid, err)
         call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
     
         call h5dcreate_f(fid, ParamName, H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
         call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE, ParamValue, data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5sclose_f(dspace_id,err)
         call h5fclose_f(fid,err)
         call h5close_f(err)
	 end subroutine SaveParam_DP
	 
     subroutine SaveParam_INT(ParamName, ParamValue)
         character(len=*) :: ParamName
		 integer :: ParamValue
		 
		 fname=trim(data_folder)//"/param"
         call h5open_f(err)
         rank=1
         data_dim1(1)=1
         call h5fopen_f(fname,H5F_ACC_RDWR_F,fid, err)
         call h5screate_simple_f(rank,data_dim1,dspace_id,err)     
     
         call h5dcreate_f(fid, ParamName, H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
         call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER, ParamValue, data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5sclose_f(dspace_id,err)
         call h5fclose_f(fid,err)
         call h5close_f(err)
	 end subroutine SaveParam_INT
	 
		 
 !---------------------------------------------------------------------------------------------------------------------------------
 !  Save total energy in fields and particles 
 !---------------------------------------------------------------------------------------------------------------------------------	 
	 
     subroutine save_total_energy
       real(dbpsn), dimension(:), allocatable :: energy, energy_this
	   
	   allocate(energy(Nflvr+2),energy_this(Nflvr+2))
       energy=0
       call CalcEnergy(energy_this)
       call MPI_REDUCE(energy_this,energy,Nflvr+2,MPI_DOUBLE_PRECISION,mpi_sum,0,MPI_COMM_WORLD)
	   
       if(proc.eq.0) then 
	       write(fname,'(I0)') t
	       fname=trim(data_folder)//"/energy_"//trim(fname)
	       call h5open_f(err)
	       rank=1
	       data_dim1(1)=1
	       call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
	       call h5screate_simple_f(rank,data_dim1,dspace_id,err)  
		  
		   call WriteReal(fid,dspace_id,'Benergy',energy(1))
		   call WriteReal(fid,dspace_id,'Eenergy',energy(2))
        
           call h5sclose_f(dspace_id,err)
	   
           data_dim1(1)=Nflvr
		   call h5screate_simple_f(rank,data_dim1,dspace_id,err)  
		   call WriteArrReal1(FileID=fid,sizex=Nflvr,var=energy(3:Nflvr+2),varname='PrtlEnergy')
	       call h5sclose_f(dspace_id,err)
		   
	       call h5fclose_f(fid,err)
	       call h5close_f(err)
         
       end if
	   deallocate(energy,energy_this)
     end subroutine save_total_energy
     
	 
	 
     subroutine save_performance_data
          INTEGER(HID_T) :: dspace_this

          write(str1,'(I0)') t
          !write(str2,'(I0)') proc
          fname=trim(data_folder)//"/prfm/"//trim(str1)
          rank=1
          call h5open_f(err)
          call h5pcreate_f(H5P_FILE_ACCESS_F, plist, err)
          call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD%MPI_VAL, MPI_INFO_NULL%MPI_VAL , err)
          call h5fcreate_f(fname, H5F_ACC_TRUNC_F, fid, err, access_prp = plist)
          call h5pclose_f(plist, err)

          data_dim1(1)=nproc
          call h5screate_simple_f(rank, data_dim1, dspace_id, err)
          data_dim1(1)=1
          call h5screate_simple_f(rank, data_dim1, memspace, err)
          
          offset1(1)=proc
        
          call h5pcreate_f(H5P_DATASET_XFER_F, plist, err) 
          call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, err)   
		  
		  call WriteReal_Collective_SP(fid,'tot',real(exec_time(31)),dspace_id,memspace,plist,offset1)
		  
		  call WriteReal_Collective_SP(fid,'movdep',real(exec_time(4)),dspace_id,memspace,plist,offset1)
          
		  call WriteReal_Collective_SP(fid,'movfld',real(exec_time(18)+exec_time(2)+exec_time(5)+exec_time(7)),dspace_id,memspace,plist,offset1)
		
		  call WriteReal_Collective_SP(fid,'prtlexchange',real(exec_time(12)),dspace_id,memspace,plist,offset1)
          
		  call WriteReal_Collective_SP(fid,'fldexchange',real(exec_time(18)),dspace_id,memspace,plist,offset1)
		  
          call WriteINT_Collective(fid,'np',np,dspace_id,memspace,plist,offset1)
		                      
          call h5sclose_f(dspace_id, err)     
          call h5pclose_f(plist, err)
          call h5sclose_f(memspace, err)
          call h5fclose_f(fid,err)
          call h5close_f(err)     
     end subroutine save_performance_data
	 
     
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

     call h5dcreate_f(fid,'tagp',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,tagp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'procp',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,procp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)

     call h5dcreate_f(fid,'qmp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,qmp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'wtp',h5psn,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,h5psn,wtp(1:used_prtl_arr_size),data_dim1,err)
     call h5dclose_f(dset_id,err)


     
	call h5sclose_f(dspace_id,err)	 
     call h5fclose_f(fid,err)
end subroutine save_particles_restart

subroutine save_param_restart
     INTEGER(HID_T) :: h5psn
	 integer :: n, i,j
	 integer, dimension(:), allocatable :: kproc_count
	 real(dbpsn), dimension(6) :: BC_pos_fld, BC_pos_prtl
	 real(psn), dimension(6) :: BC_speed
	 character(len=16):: vname
      integer, dimension(:), allocatable :: int_arr
	 
	 write(str1,'(I0)') t
     write(str2,'(I0)') proc
     fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
     rank=1
     data_dim1(1)=1
     call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
	 
	 !-----  scalar variables ---------! 
     call h5screate_simple_f(rank,data_dim1,dspace_id,err)     

     call WriteINT(fid,dspace_id,'mx',mx)
     call WriteINT(fid,dspace_id,'my',my)
     call WriteINT(fid,dspace_id,'mz',mz)
     call WriteINT(fid,dspace_id,'prtl_arr_size',prtl_arr_size)
     call WriteINT(fid,dspace_id,'used_prtl_arr_size',used_prtl_arr_size)
	 
     call WriteINT(fid,dspace_id,'np',np)
     call WriteINT(fid,dspace_id,'Nflvr',Nflvr)
	 
     call WriteINT(fid,dspace_id,'proc',proc)
	 
	 
	 !subdomain grid, borders, and proc. layout
     call WriteINT(fid,dspace_id,'nSubDomainsX',nSubDomainsX)
     call WriteINT(fid,dspace_id,'nSubDomainsY',nSubDomainsY)
     call WriteINT(fid,dspace_id,'nSubDomainsZ',nSubDomainsZ)
	 call WriteINT(fid,dspace_id,'isizeProcGrid',isizeProcGrid)
	 call WriteINT(fid,dspace_id,'jsizeProcGrid',jsizeProcGrid)
	 call WriteINT(fid,dspace_id,'iproc',iproc) 
	 call WriteINT(fid,dspace_id,'jproc',jproc) 
	 call WriteINT(fid,dspace_id,'kproc',kproc) 
	 call WriteINT(fid,dspace_id,'procxind',procxind)
	 call WriteINT(fid,dspace_id,'procyind',procyind)
	 call WriteINT(fid,dspace_id,'proczind',proczind)
	 call WriteINT(fid,dspace_id,'indepLBaxis',indepLBaxis)
	 
	
	 call h5sclose_f(dspace_id,err)
	 
	 !--------- 1D arrays ----------------------------------------
	 
	 call WriteINT1(fid,ProcGrid(iproc,jproc)%count+1,'borders',ProcGrid(iproc,jproc)%borders(0:ProcGrid(iproc,jproc)%count))
	 call WriteINT1(fid,nSubDomainsX+1,'xborders',xborders(0:nSubDomainsX))
	 call WriteINT1(fid,nSubDomainsY+1,'yborders',yborders(0:nSubDomainsY))
	 call WriteINT1(fid,nSubDomainsZ+1,'zborders',zborders(0:nSubDomainsZ))
	 call WriteINT1(fid,6,'box_bounds',box_bounds)
	 
	 !-------- Proc Grid -----------------------------------------
	 if(proc.eq.0) then !only master saves the proc grid data
		 
		 !save proc count
		 allocate(kproc_count(isizeProcGrid*jsizeProcGrid))
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				 kproc_count(j*isizeProcGrid + i+1) = ProcGrid(i,j)%count
			 end do 
		 end do 
		 call WriteINT1(fid,isizeProcGrid*jsizeProcGrid,'kproc_count',kproc_count)
		 deallocate(kproc_count)
		 
		 ! save proc list for every k-dimension of the ProcGrid array
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				  if(ProcGrid(i,j)%count.gt.0) then 
		               write(vname,'(I0)') j*isizeProcGrid + i+1
		               vname="procs"//trim(vname)
					   call WriteINT1(fid,ProcGrid(i,j)%count,vname,ProcGrid(i,j)%procs(0:ProcGrid(i,j)%count-1))
				  end if
			 end do 
		 end do 
		 
	 end if 
	 
	 
	 !save BC vars	 
	 do n=1,6
		 BC_Pos_Fld(n) = bc_face(n)%pos_fld
		 BC_Pos_Prtl(n) = bc_face(n)%pos_prtl
		 BC_speed(n) = bc_face(n)%speed 
	 end do 
	 call WriteArrReal1(FileID=fid,sizex=6,var=BC_Pos_Fld,varname='BC_Pos_Fld')
	 call WriteArrReal1(FileID=fid,sizex=6,var=BC_Pos_Prtl,varname='BC_Pos_Prtl')
	 call WriteArrReal1(FileID=fid,sizex=6,var=BC_speed,varname='BC_speed')
	 
	 
	 !Flavour attributes for all species
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
     call h5dcreate_f(fid,'FlvrSaveFldData',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,FlvrSaveFldData,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrType',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,FlvrType,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'FlvrSaveRatio',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,FlvrSaveRatio,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     call h5dcreate_f(fid,'CurrentTagID',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagID,data_dim1,err)
     call h5dclose_f(dset_id,err)
     call h5dcreate_f(fid,'CurrentTagProcID',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagProcID,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     
     !save other flavor properties 
     allocate(int_arr(Nflvr))
     int_arr = [(flvr_prpt(i)%ionization , i=1, Nflvr)]
     call h5dcreate_f(fid,'ionization_type',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
     call h5dclose_f(dset_id,err) 

     int_arr = [(flvr_prpt(i)%ionization_elc_flv , i=1, Nflvr)]
     call h5dcreate_f(fid,'ionization_elc_flv',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
     call h5dclose_f(dset_id,err) 

     int_arr = [(flvr_prpt(i)%Z_nucleus , i=1, Nflvr)]
     call h5dcreate_f(fid,'ionization_z_nucleus',H5T_NATIVE_INTEGER,dspace_id,dset_id,err)
     call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
     call h5dclose_f(dset_id,err) 
     deallocate(int_arr)

     
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

    call WriteINT(fid,dspace_id,'restart_time',t)
    call h5sclose_f(dspace_id,err)
	call h5fclose_f(fid,err)	
end subroutine save_restart_time

!============================================End of restart data save subroutines ====================================================

      

end module savedata