module reload
     use parameters
     use vars
     use hdf5
	 use prtl_tag
	 use mem_prtl
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

     subroutine RestartInitParam

          fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
          call h5open_f(err)
		  
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
		  call HDF5readINT(fid,'mx',mx)
		  call HDF5readINT(fid,'my',my)
		  call HDF5readINT(fid,'mz',mz)
		  call HDF5readINT(fid,'prtl_arr_size',prtl_arr_size)
		  
		  call InitPrtlArr(prtl_arr_size) !allocate prtl memory
		  
		  call HDF5readINT(fid,'used_prtl_arr_size',used_prtl_arr_size)
		  call HDF5readINT(fid,'np',np)
		 
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
		  integer, dimension(:), allocatable :: int_arr
		  integer :: i

          call h5open_f(err)
		  
          fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
          data_dim1(1)=1
          call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)		  	  
          call h5dopen_f(fid,'Nflvr', dset_id, err)		  
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,Nflvr,data_dim1,err)
          call h5dclose_f(dset_id,err) 
		  
          allocate(flvrqm(Nflvr),FlvrCharge(Nflvr),FlvrSaveFldData(Nflvr),FlvrType(Nflvr),FlvrSaveRatio(Nflvr))
          allocate(CurrentTagID(Nflvr),CurrentTagProcID(Nflvr))
		  allocate(flvr_prpt(NFlvr))

          data_dim1(1)=Nflvr
          call h5dopen_f(fid,'flvrqm', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,flvrqm,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrCharge', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,FlvrCharge,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrSaveFldData', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,FlvrSaveFldData,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrType', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,FlvrType,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'FlvrSaveRatio', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,FlvrSaveRatio,data_dim1,err)
          call h5dclose_f(dset_id,err) 
          call h5dopen_f(fid,'CurrentTagID', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagID,data_dim1,err)
          call h5dclose_f(dset_id,err)
		  call h5dopen_f(fid,'CurrentTagProcID', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,CurrentTagProcID,data_dim1,err)
          call h5dclose_f(dset_id,err)

		  allocate(int_arr(NFlvr))
		  call h5dopen_f(fid,'ionization_type', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
          call h5dclose_f(dset_id,err)
		  forall (i=1:Nflvr) flvr_prpt(i)%ionization = int_arr(i)	
		  
		  call h5dopen_f(fid,'ionization_elc_flv', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
          call h5dclose_f(dset_id,err)
		  forall (i=1:Nflvr) flvr_prpt(i)%ionization_elc_flv = int_arr(i)	
	

		  call h5dopen_f(fid,'ionization_z_nucleus', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,int_arr,data_dim1,err)
          call h5dclose_f(dset_id,err)
		  forall (i=1:Nflvr) flvr_prpt(i)%Z_nucleus = int_arr(i)
		  deallocate(int_arr)	
	
		  		  

		  
		  
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

		  call h5dopen_f(fid,'tagp', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,tagp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
		  call h5dopen_f(fid,'procp', dset_id, err)
          call h5dread_f(dset_id,H5T_NATIVE_INTEGER,procp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)

		  call h5dopen_f(fid,'qmp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,qmp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
		  call h5dopen_f(fid,'wtp', dset_id, err)
          call h5dread_f(dset_id,h5real_psn,wtp(1:used_prtl_arr_size),data_dim1,err)
          call h5dclose_f(dset_id,err)
		  call h5fclose_f(fid,err)



          call h5close_f(err)
     end subroutine RestartInitPrtl
	 
	 subroutine RestartAllVars
		 INTEGER(HSIZE_T), dimension(1) :: data_dim1
		 integer :: n, i,j
		 integer, dimension(:), allocatable :: kproc_count
		 real(dbpsn), dimension(6) :: BC_pos_fld, BC_pos_prtl
		 real(psn), dimension(6) :: BC_speed
		 character(len=16):: vname

         fname=trim(data_folder)//"/restart"//"/param_"//trim(str2)//"_"//trim(str1)
         call h5open_f(err)
         call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)

	 
		 call HDF5readINT(fid,'indepLBaxis',indepLBaxis)

		 call HDF5readINT(fid,'procxind',procxind)
		 call HDF5readINT(fid,'procyind',procyind)
		 call HDF5readINT(fid,'proczind',proczind)
		 call HDF5readINT(fid,'nSubDomainsX',nSubDomainsX)
		 call HDF5readINT(fid,'nSubDomainsY',nSubDomainsY)
		 call HDF5readINT(fid,'nSubDomainsZ',nSubDomainsZ)
		 
		 !x,y,z borders
         call h5dopen_f(fid,'xborders', dset_id, err)
		 data_dim1(1)=nSubDomainsX+1
		 allocate(xborders(0:nSubDomainsX))
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,xborders(0:nSubDomainsX),data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5dopen_f(fid,'yborders', dset_id, err)
		 data_dim1(1)=nSubDomainsY+1
		 allocate(yborders(0:nSubDomainsY))
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,yborders(0:nSubDomainsY),data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
         call h5dopen_f(fid,'zborders', dset_id, err)
		 allocate(zborders(0:nSubDomainsZ))
		 data_dim1(1)=nSubDomainsZ+1
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,zborders(0:nSubDomainsZ),data_dim1,err)
         call h5dclose_f(dset_id,err)
		 
		 call HDF5readINT(fid,'iproc',iproc)
		 call HDF5readINT(fid,'jproc',jproc)
		 call HDF5readINT(fid,'kproc',kproc)
		 call HDF5readINT(fid,'isizeProcGrid',isizeProcGrid)
		 call HDF5readINT(fid,'jsizeProcGrid',jsizeProcGrid)
		 
		 allocate(ProcGrid(0:isizeProcGrid-1,0:jsizeProcGrid-1))

		 !1D proc count is read by the master proc. and broadcasted to other proc.
		 allocate(kproc_count(isizeProcGrid*jsizeProcGrid))
		 
		 if(proc.eq.0) then 

			 data_dim1(1)=isizeProcGrid*jsizeProcGrid
			 call h5dopen_f(fid,'kproc_count', dset_id, err)
			 call h5dread_f(dset_id,H5T_NATIVE_INTEGER,kproc_count,data_dim1,err)
			 call h5dclose_f(dset_id,err)
			 
		 end if 
		 
		 call MPI_Bcast( kproc_count ,isizeProcGrid*jsizeProcGrid, MPI_INTEGER, 0, MPI_COMM_WORLD)
		 
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				 ProcGrid(i,j)%count = kproc_count(j*isizeProcGrid + i+1)
			 end do
		 end do
		
		 deallocate(kproc_count)
		 
		 !read/broadcast proc list for every k-dimension of the ProcGrid array
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				 if(ProcGrid(i,j)%count.gt.0) then 
					 
					 allocate(ProcGrid(i,j)%procs(  0:ProcGrid(i,j)%count-1 ))
					 
					 if(proc.eq.0) then 

						 data_dim1(1)=ProcGrid(i,j)%count
  		                 write(vname,'(I0)') j*isizeProcGrid + i+1
  		                 vname="procs"//trim(vname)
						 call h5dopen_f(fid,vname, dset_id, err)
						 call h5dread_f(dset_id,H5T_NATIVE_INTEGER, ProcGrid(i,j)%procs(  0:ProcGrid(i,j)%count-1 ), data_dim1, err)
						 call h5dclose_f(dset_id,err)
			 
					 end if 
					 
					 call MPI_Bcast( ProcGrid(i,j)%procs(  0:ProcGrid(i,j)%count-1 ) , ProcGrid(i,j)%count , MPI_INTEGER, 0, MPI_COMM_WORLD)
					 
				 end if
				 
			 end do 
		 end do
		 

		 allocate(ProcGrid(iproc,jproc)%borders(  0:ProcGrid(iproc,jproc)%count ))		 
         call h5dopen_f(fid,'borders', dset_id, err)
		 data_dim1(1)=ProcGrid(iproc,jproc)%count+1
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,ProcGrid(iproc,jproc)%borders(0:ProcGrid(iproc,jproc)%count),data_dim1,err)
         call h5dclose_f(dset_id,err)
		 	 
         call h5dopen_f(fid,'box_bounds', dset_id, err)
		 data_dim1(1)=6
         call h5dread_f(dset_id,H5T_NATIVE_INTEGER,box_bounds,data_dim1,err)
         call h5dclose_f(dset_id,err)
		 		 
		 call HDF5readRealDP_1D(fid,'BC_Pos_Fld',BC_pos_fld,6)
		 call HDF5readRealDP_1D(fid,'BC_Pos_Prtl',BC_pos_prtl,6)
		 call HDF5readReal_1D(fid,'BC_speed',BC_speed,6)
		 do n=1,6
			 bc_face(n)%pos_fld = BC_pos_fld(n)
			 bc_face(n)%pos_prtl = BC_pos_prtl(n)
			 bc_face(n)%speed = BC_speed(n)
		 end do

		 
         
		 call h5fclose_f(fid,err)
		 call h5close_f(err)
		 
	 end subroutine RestartAllVars
     
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
	      real(psn)                      :: var
		  integer                        :: err
	      data_dim1(1)=1
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,h5real_psn,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readReal
	 
	 subroutine HDF5readReal_1D(fid,varname,var,size)
	      INTEGER(HID_T)                 :: fid,dset_id 
	      INTEGER(HSIZE_T), dimension(1) :: data_dim1
	      character(len=*)               :: varname
		  integer                        :: err, size
		  real(psn), dimension(size)     :: var

	      data_dim1(1)=size
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,h5real_psn,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readReal_1D
	 
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
	 
	 subroutine HDF5readRealDP_1D(fid,varname,var,size)
	      INTEGER(HID_T)                 :: fid,dset_id 
	      INTEGER(HSIZE_T), dimension(1) :: data_dim1
	      character(len=*)               :: varname
		  integer                        :: size,err
		  real(dbpsn), dimension(size)   :: var

	      data_dim1(1)=size
		  call h5dopen_f(fid,varname, dset_id, err)
	      call h5dread_f(dset_id,h5real_dbpsn,var,data_dim1,err)
	      call h5dclose_f(dset_id,err) 
	 end subroutine HDF5readRealDP_1D
	 
end module reload