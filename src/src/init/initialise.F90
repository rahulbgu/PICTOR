module initialise
     use parameters
     use vars
	 use prtl_stats
	 use reload 
	 use mem_prtl
	 use mem_fld
	 use prtl_tag
	 use subdomains
	 use init_prtl
#ifdef OPEN_MP
     use omp_lib 
#endif 	 
     implicit none
contains
     ! this subroutine intializes location of all particles 
     subroutine InitAll		  

          call AllocateFldVars
		call InitFlds_default

		call AllocatePrtlVars 

     end subroutine InitAll
	  
	 subroutine InitOpenMP		 
#ifdef OPEN_MP
!$OMP PARALLEL
         ThreadID=OMP_GET_THREAD_NUM()
         if(ThreadID.eq.0) Nthreads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#else 
         ThreadID=0
         Nthreads=1
#endif 
		 if(proc.eq.0) print*,Nthreads,'CPU threads are active within each MPI proc'
		 
	 end subroutine InitOpenMP
	 
	 subroutine InitPrtlScan
         call InitScanPrtlArr !set np, tag flv=1,2, and check for errors, if possible   
	 end subroutine InitPrtlScan    
	 
	 subroutine InitPrtlArrSize 
		  integer :: Nelc 

          Nelc = Nelc_uniform()         
	  		  
		prtl_arr_size= min( int(3.0*Nelc) , 10000000 ) ! max. initial prtl size is 10M  
#ifdef gpu
          prtl_arr_size=gpu_prtl_arr_size
#endif 
	 end subroutine InitPrtlArrSize
	 
	 
     subroutine AllocateFldVars
          allocate(Ex(mx,my,mz),Ey(mx,my,mz),Ez(mx,my,mz))
          allocate(Bx(mx,my,mz),By(mx,my,mz),Bz(mx,my,mz))
		  
		call InitAuxFld
		  
     end subroutine AllocateFldVars
	 

     subroutine AllocatePrtlVars
		Nflvr=2 !by default the number of flv is 2, it is increased on demand in SetQbyM(help_setup.F90)
          allocate(flvrqm(Nflvr),FlvrCharge(Nflvr),FlvrSaveFldData(Nflvr),FlvrType(Nflvr),FlvrSaveRatio(Nflvr))
          allocate(CurrentTagID(Nflvr),CurrentTagProcID(Nflvr))
          allocate(flvr_prpt(Nflvr))

          !Default values 
          flvrqm(1)=qmi
          flvrqm(2)=qme
		FlvrCharge(1)=1.0_psn
		FlvrCharge(2)=-1.0_psn
          FlvrSaveFldData(1)=1
          FlvrSaveFldData(2)=1
          FlvrType(1)=0
          FlvrType(2)=0
          FlvrSaveRatio(1)=psave_ratio
          FlvrSaveRatio(2)=psave_ratio  
          CurrentTagID(1) = 0 
          CurrentTagID(2) = 0 
          CurrentTagProcID(1) = proc + 1
          CurrentTagProcID(2) = proc + 1 
           
	  
          call InitPrtlArr(prtl_arr_size)		  
		  
     end subroutine AllocatePrtlVars
    
	 
	 subroutine InitBoxDefaultBounds
		 box_bounds(1) = 0 
		 box_bounds(2) = nx
		 box_bounds(3) = 0 
		 box_bounds(4) = ny
		 box_bounds(5) = 0 
		 box_bounds(6) = nz
#ifdef twoD
		 box_bounds(5) = 0 
		 box_bounds(6) = 1
#endif		 
	 end subroutine InitBoxDefaultBounds
	 
	 subroutine InitBC_Face
		  integer :: n
		  do n=1,6
			  if(mod(n,2).eq.0) then 
				  bc_face(n)%pos_fld=huge(1)
				  bc_face(n)%pos_prtl=huge(1)
			  else
				  bc_face(n)%pos_fld=-huge(1)
				  bc_face(n)%pos_prtl=-huge(1) 
			  end if

			  bc_face(n)%type_fld='prdc' !default BC is periodic
			  bc_face(n)%type_prtl='prdc' !default BC is periodic
			  bc_face(n)%speed=0
			  bc_face(n)%attn_thickness=0 
		  end do
		  
	 end subroutine InitBC_Face

	 
! By default the intial elemectromagnetic fields is set to 0, but they can be reset in the setup file  
subroutine InitFlds_default
     integer :: i,j,k
     Ex=0.0_psn; Ey=0.0_psn; Ez=0.0_psn; Bx=0.0_psn; By=0.0_psn; Bz=0.0_psn
     Bx_ext0=0.0_psn; By_ext0=0.0_psn; Bz_ext0=0.0_psn
end subroutine InitFlds_default


subroutine InitScanPrtlArr
     integer :: n,count
	 real(psn) :: r1
	 real(psn) :: tag_fraction
     !count the total number of (active) particles
     count=0
     do n=1,prtl_arr_size
          if(flvp(n).ne.0) then 
			   count=count+1
			   used_prtl_arr_size=n
		  end if 
     end do 
     np=count
	 
end subroutine InitScanPrtlArr


subroutine InitRandomNumberSeed
     !this subrotines tries to get unique random seed for all proc
     integer :: seed_size
     integer, allocatable, dimension(:) :: seed_arr
     real(kind=8) :: r1
     integer :: i, clock
	 
	 !call system_clock(clock)
     call random_seed(SIZE=seed_size)
     allocate(seed_arr(seed_size))
	 
	 seed_arr(1) = 2930
     do i=2,seed_size
         seed_arr(i)= mod( 19*seed_arr(i-1) + 7*proc + 3457, 11399 )
     end do
     call random_seed(PUT=seed_arr)
     call random_number(r1)
    
	 
	 do i=1,1000
		 call random_number(r1)
	 end do
	 
!     call random_seed(PUT=seed_arr)
     
end subroutine InitRandomNumberSeed

          
          
end module initialise