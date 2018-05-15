module initialise
     use parameters
     use vars
	 use prtl_stats
	 use reload 
     use setup
#ifdef OPEN_MP
     use omp_lib 
#endif 	 
     implicit none
contains
     ! this subroutine intializes location of all particles 
     subroutine InitAll
          call InitRandomNumberSeed

#ifdef OPEN_MP
!$OMP PARALLEL
ThreadID=OMP_GET_THREAD_NUM()
if(ThreadID.eq.0) Nthreads = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
#else 
ThreadID=0
Nthreads=1
#endif 


if(proc.eq.0) print*,Nthreads,' CPU threads are active within each MPI proc'
		  
          if(restart) then 
               call RestartInitParam
          else
               call InitParam
          end if
          
          call InitPrtlBoundaries
          call InitScale
          call AllocateFldVars
          if(.not.restart) call AllocatePrtlVars !must be different for restart 
           if(restart) then
               call RestartInitFlds			   
               call RestartInitPrtl
			   call RestartAllVars
               call InitUser ! Init user must be designed for restart
           else
               call InitFlds_default
               call Initboundaries_default 
               call InitProcXYZind 
               call InitProcGrid 
               call InitUser
               call InitFldSaveDataLimit_default !limit to save fld data     
               call InitScanPrtlArr !set np, tag flv=1,2, and check for errors if possible   
			   call InitScanTestPrtlArr    
          end if
		  call InitPrtlStatVars
          call InitCurrent
          call InitTransferVarsSize !Warning :: load from restart file, if restarted
          call AllocateTransferVars
          if(.not.restart) call InitOverride
           !print*,'prtl_arr_size is',prtl_arr_size,'np',np	
     end subroutine InitAll

     subroutine InitScale
          ompe=c/c_ompe
          qi=(g0*ompe**2)/(epc*(1.0_psn+me/mi))
          qe=-qi
          qme=-1.0_psn
          qmi=1.0_psn/mi
          masse=abs(qe)*me
          massi=qi*mi
          !print*,'charge.mass is',masse,massi
     end subroutine InitScale

! #if defined(Hybrid1) || defined(Hybrid2) || defined(Hybrid3)
!      subroutine InitScale
!           ompe=c/c_ompe
!           qi=(g0*ompe**2)/ppc
!           qe=-qi
!           qme=0
!           qmi=1
!           masse=0
!           massi=qi
!           !print*,'charge.mass is',masse,massi
!      end subroutine InitScale
! #endif
    
     subroutine InitCurrent
          Jx=0
          Jy=0
          Jz=0
     end subroutine InitCurrent
     subroutine InitTransferVarsSize
          loutp_size=10*my*mz*epc*2
          routp_size=10*my*mz*epc*2
          rinp_size=10*my*mz*epc*2
          linp_size=10*my*mz*epc*2
          toutp_size=10*mx*mz*epc*2
          boutp_size=10*mx*mz*epc*2
          tinp_size=10*mx*mz*epc*2
          binp_size=10*mx*mz*epc*2
          uoutp_size=10*mx*my*epc*2
          doutp_size=10*mx*my*epc*2
          uinp_size=10*mx*my*epc*2
          dinp_size=10*mx*my*epc*2
		  linp_count=0
		  rinp_count=0
		  binp_count=0
		  tinp_count=0
		  uinp_count=0
		  dinp_count=0
		  lintp_count=0
		  rintp_count=0
		  bintp_count=0
		  tintp_count=0
		  uintp_count=0
		  dintp_count=0
     end subroutine InitTransferVarsSize

     subroutine Initboundaries_default
          integer :: i
          do i=0,nSubDomainsX
               xborders(i)=i*(mx-5)
          end do
          do i=0,nSubDomainsY
               yborders(i)=i*(my-5)
          end do
#ifdef twoD
               zborders(0)=0
			   zborders(1)=1
#else                          
          do i=0,nSubDomainsZ
               zborders(i)=i*(mz-5)
          end do
#endif               		  
     end subroutine Initboundaries_default
     subroutine InitProcXYZind
          integer :: i
#ifdef twoD          
          do i=0,nSubDomainsX*nSubDomainsY-1
               procxind(i)=i-nSubDomainsX*int(i/nSubDomainsX)
               procyind(i)=int(i/nSubDomainsX)
               proczind(i)=0
          end do
#else
          do i=0,nSubDomainsX*nSubDomainsY*nSubDomainsZ-1
               procxind(i)=i-nSubDomainsX*int(i/nSubDomainsX)
               procyind(i)=int(i/nSubDomainsX)-nSubDomainsY*int(i/(nSubDomainsX*nSubDomainsY))
               proczind(i)=int(i/(nSubDomainsX*nSubDomainsY))
          end do
#endif           
     end subroutine InitProcXYZind
     subroutine InitProcGrid
         integer :: i,j,k
#ifdef twoD
        do k=0,0
#else                 
        do k=0,nSubDomainsZ-1
#endif               
              do j=0,nSubDomainsY-1
                  do i=0,nSubDomainsX-1
                        proc_grid(i,j,k)=i+j*nSubDomainsX+k*nSubDomainsX*nSubDomainsY
                  end do
              end do
         end do
     end subroutine InitProcGrid
     subroutine AllocateFldVars
          allocate(Ex(mx,my,mz),Ey(mx,my,mz),Ez(mx,my,mz))
          allocate(Bx(mx,my,mz),By(mx,my,mz),Bz(mx,my,mz))
          allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz))
          allocate(F0(mx,my,mz))
          if(nMoverEMfilter.gt.0) then 
               allocate(FilteredEx(mx,my,mz),FilteredEy(mx,my,mz),FilteredEz(mx,my,mz))
          end if 
		  
     end subroutine AllocateFldVars

     subroutine AllocatePrtlVars
		  Nflvr=2 !by default the number of flv is 2, it is increased on demand in SetQbyM(help_setup.F90)
          allocate(flvrqm(Nflvr),FlvrSaveFldData(Nflvr),FlvrType(Nflvr),FlvrSplitPrtl(Nflvr))
          allocate(CurrentTagID(Nflvr),TagCounter(Nflvr))
          !Default values 
          flvrqm(1)=qmi
          flvrqm(2)=qme
          FlvrSaveFldData(1)=1
          FlvrSaveFldData(2)=1
          FlvrType(1)=0
          FlvrType(2)=0
          FlvrSplitPrtl(1)=0
          FlvrSplitPrtl(2)=0   
	 	  CurrentTagID(1)=1+proc*NtagProcLen !initialise the value of current tag 
	 	  CurrentTagID(2)=1+proc*NtagProcLen  
		  TagCounter(1)=1
		  TagCounter(2)=1
#ifndef twoD
          Nelc=epc*(mx-5)*(my-5)*(mz-5) ! Initial Number of electrons  
#else
          Nelc=epc*(mx-5)*(my-5)
#endif
          prtl_arr_size=int(3.0*Nelc)
		  
#ifdef gpu 
          prtl_arr_size=prtl_arr_size+gpu_prtl_arr_size 
#endif
#ifdef GPU_EXCLUSIVE
          prtl_arr_size=gpu_prtl_arr_size 
#endif 	  
		  !Allocate Prtl Arr
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
		  used_prtl_arr_size=0
		  prtl_random_insert_index=1
		  np=0 
		  allocate(SortedPrtlCountYZ(my,mz))
		  SortedPrtlCountYZ=0
		  !Allocate test prtl array size 
		  test_prtl_arr_size=Nelc/10
#ifdef gpu
          test_prtl_arr_size=test_prtl_arr_size+gpu_test_prtl_arr_size 
#endif
#ifdef GPU_EXCLUSIVE
          test_prtl_arr_size=gpu_test_prtl_arr_size
#endif 		  
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
		  used_test_prtl_arr_size=0
		  test_prtl_random_insert_index=1
		  ntp=0  
     end subroutine AllocatePrtlVars

     subroutine AllocateTransferVars     
          !arrays used to send recieve particles outliers     
          allocate(loutp(loutp_size),linp(linp_size),routp(routp_size),rinp(rinp_size))           
          allocate(toutp(toutp_size),tinp(tinp_size),boutp(boutp_size),binp(binp_size))           
          
          allocate(buff_lJx(3,my,mz),buff_rJx(3,my,mz),buff_lJy(3,my,mz),buff_rJy(3,my,mz),buff_lJz(3,my,mz),buff_rJz(3,my,mz))
          allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz))
#ifndef twoD
        allocate(uoutp(uoutp_size),uinp(uinp_size),doutp(doutp_size),dinp(dinp_size))           
        allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))
#endif           
     end subroutine AllocateTransferVars
     
     subroutine InitPrtlBoundaries
          xmin=3
          xmax=mx-2
          ymin=3
          ymax=my-2
#ifdef twoD
          zmin=1
          zmax=2
#else          
          zmin=3
          zmax=mz-2
#endif
          xlen=xmax-xmin
          ylen=ymax-ymin
          zlen=zmax-zmin
     end subroutine InitPrtlBoundaries

     subroutine InitParam	 
           mx=nx/nSubDomainsX+5
           my=ny/nSubDomainsY+5
#ifdef twoD
          mz=1 
#else
          mz=nz/nSubDomainsZ+5
#endif
           tstart=1     
           
     end subroutine InitParam
! By default the intial elemectromagnetic fields is set to 0, but they can be reset in the setup file  
subroutine InitFlds_default
     integer :: i,j,k
     do i=1,mx
          do j=1,my
               do k=1,mz
                    Ex(i,j,k)=0.0_psn
                    Ey(i,j,k)=0.0_psn
                    Ez(i,j,k)=0.0_psn
                    Bx(i,j,k)=0.0_psn
                    By(i,j,k)=0.0_psn
                    Bz(i,j,k)=0.0_psn
               end do
          end do
     end do
     Bx_ext0=0.0_psn
     By_ext0=0.0_psn
     Bz_ext0=0.0_psn
end subroutine InitFlds_default

subroutine InitFldSaveDataLimit_default
     fdataxi_box=xborders(0)
     fdataxf_box=xborders(nSubDomainsX)
     fdatayi_box=yborders(0)
     fdatayf_box=yborders(nSubDomainsY)
#ifndef twoD
     fdatazi_box=zborders(0)
     fdatazf_box=zborders(nSubDomainsZ)
#else     
     fdatazi_box=0
     fdatazf_box=0
#endif
     
end subroutine InitFldSaveDataLimit_default

subroutine InitScanPrtlArr
     integer :: n,count
     !count the total number of (active) particles
     count=0
     do n=1,prtl_arr_size
          if(qp(n).ne.0) then 
			   count=count+1
			   used_prtl_arr_size=n
		  end if 
     end do 
     np=count
	 prtl_random_insert_index=used_prtl_arr_size+1
	 
	 !Tag Some Selected Particles 
	 !CurrentTagID(1)=1+proc*NtagProcLen !initialise the value of current tag 
	 !CurrentTagID(2)=1+proc*NtagProcLen
!      if(psave_ratio.gt.0) then !tag the flv=1,2 particles by default
!           do n=1,Nelc,psave_ratio !This is not corret, electron and ions could be anywhere
!                if(qp(n).ne.0) then
! 				   tagp(n)=CurrentTagID(1)
! 				   if(mod(CurrentTagID(1),NtagProcLen).eq.0) CurrentTagID(1)=CurrentTagID(1)+NtagProcLen*(nproc-1)
! 			       CurrentTagID(1)=CurrentTagID(1)+1
! 		       end if
!           end do
!
!           do n=1,Nelc,psave_ratio
!                if(qp(n+Nelc).ne.0) then
! 				   tagp(n+Nelc)=CurrentTagID(2)
!  				   if(mod(CurrentTagID(2),NtagProcLen).eq.0) CurrentTagID(2)=CurrentTagID(2)+NtagProcLen*(nproc-1)
!  			       CurrentTagID(2)=CurrentTagID(2)+1
! 				end if
!           end do
!      end if

end subroutine InitScanPrtlArr

subroutine InitScanTestPrtlArr
    integer :: n,count
    count=0
    do n=1,test_prtl_arr_size
         if(qtp(n).ne.0) then
			 count=count+1
			 used_test_prtl_arr_size=n
		 end if 
    end do 
	ntp=count
    test_prtl_random_insert_index=used_test_prtl_arr_size+1
end subroutine InitScanTestPrtlArr




subroutine InitRandomNumberSeed
     !this subrotines tries to get unique random seed for all proc
     integer :: seed_size
     integer, allocatable, dimension(:) :: seed_arr
     real(kind=8) :: r1
     integer :: i
     call random_seed(SIZE=seed_size)
     allocate(seed_arr(seed_size))
     do i=1,seed_size
         seed_arr(i)=19*proc+i*207
     end do
     call random_seed(PUT=seed_arr)
     call random_number(r1)
     do i=1,int(1520*r1)
         call random_number(r1)
     end do
     seed_arr=0
     do i=1,seed_size
         seed_arr(i)=seed_arr(i)+int(100000000*r1)
          call random_number(r1)
     end do
    call random_seed(PUT=seed_arr)
     
end subroutine InitRandomNumberSeed

          
          
end module initialise