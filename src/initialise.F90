module initialise
     use parameters
     use vars
	 use prtl_stats
	 use reload 
	 use memory
#ifdef OPEN_MP
     use omp_lib 
#endif 	 
     implicit none
contains
     ! this subroutine intializes location of all particles 
     subroutine InitAll		  

          call AllocateFldVars
          call InitCurrent
          call InitTransferVars 

		  call AllocatePrtlVars 
          call InitFlds_default

     end subroutine InitAll
	 
	 subroutine InitDomainSkelton	 	 
		  tstart=1     
		  call InitNeighborsMPI
          call InitProcXYZind 
          call InitProcGrid     
          call InitFldSaveDataLimit_default !limit to save fld data   	 
	 end subroutine InitDomainSkelton
	 
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
         call InitScanPrtlArr !set np, tag flv=1,2, and check for errors if possible   
	     call InitScanTestPrtlArr    
	 end subroutine InitPrtlScan

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
	 
	 subroutine InitPrtlArrSize        
#ifndef twoD
          Nelc=epc*(mx-5)*(my-5)*(mz-5) ! Initial Number of electrons  
#else
          Nelc=epc*(mx-5)*(my-5)
#endif
          prtl_arr_size= min ( int(3.0*Nelc) , 10000000 ) ! min. initial prtl size is 10M  

#ifdef gpu
          prtl_arr_size=gpu_prtl_arr_size
#endif 
	
	 end subroutine InitPrtlArrSize
	 
     subroutine InitTransferVars
		  linp_count=0; rinp_count=0; binp_count=0; tinp_count=0; uinp_count=0; dinp_count=0
		  lintp_count=0; rintp_count=0; bintp_count=0; tintp_count=0; uintp_count=0; dintp_count=0
		                  
          allocate(buff_lJx(3,my,mz),buff_rJx(3,my,mz),buff_lJy(3,my,mz),buff_rJy(3,my,mz),buff_lJz(3,my,mz),buff_rJz(3,my,mz))
          allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz))
		  buff_lJx=0; buff_rJx=0; buff_lJy=0; buff_rJy=0; buff_lJz=0; buff_rJz=0;
		  buff_bJx=0; buff_tJx=0; buff_bJy=0; buff_tJy=0; buff_bJz=0; buff_tJz=0;
#ifndef twoD        
          allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))
		  buff_dJx=0; buff_uJx=0; buff_dJy=0; buff_uJy=0; buff_dJz=0; buff_uJz=0;
#endif   
     end subroutine InitTransferVars

     subroutine Initboundaries
          integer :: i
          do i=0,nSubDomainsX
               xborders(i)=i*(nx/nSubDomainsX)
          end do
		  xborders(nSubDomainsX)=nx
 		  
          do i=0,nSubDomainsY
               yborders(i)=i*(ny/nSubDomainsY)
          end do
		  yborders(nSubDomainsY)=ny
#ifdef twoD
               zborders(0)=0
			   zborders(1)=1
#else                          
          do i=0,nSubDomainsZ
               zborders(i)=i*(nz/nSubDomainsZ)
          end do
		  zborders(nSubDomainsZ)=nz
#endif               		  
     end subroutine Initboundaries
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
	 
	 subroutine InitNeighborsMPI
          lproc=proc-1
          rproc=proc+1
          if(modulo(proc,nSubDomainsX).eq.0) lproc=proc+nSubDomainsX-1
          if(modulo(proc+1,nSubDomainsX).eq.0) rproc=proc-(nSubDomainsX-1)
          
          tproc=proc+nSubDomainsX
          bproc=proc-nSubDomainsX
#ifdef twoD
          if(proc.ge.nSubDomainsX*(nSubDomainsY-1)) tproc=proc-nSubDomainsX*(nSubDomainsY-1)
          if(proc.lt.nSubDomainsX) bproc=proc+nSubDomainsX*(nSubDomainsY-1) 
#else
          if((proc-nSubDomainsX*nSubDomainsY*int(proc/(nSubDomainsX*nSubDomainsY))).ge.nSubDomainsX*(nSubDomainsY-1)) tproc=proc-nSubDomainsX*(nSubDomainsY-1)
          if((proc-nSubDomainsX*nSubDomainsY*int(proc/(nSubDomainsX*nSubDomainsY))).lt.nSubDomainsX) bproc=proc+nSubDomainsX*(nSubDomainsY-1)  
#endif               

#ifndef twoD
          uproc=proc+nSubDomainsX*nSubDomainsY
          dproc=proc-nSubDomainsX*nSubDomainsY
          if(proc.ge.nSubDomainsX*nSubDomainsY*(nSubDomainsZ-1)) uproc=proc-nSubDomainsX*nSubDomainsY*(nSubDomainsZ-1)
          if(proc.lt.nSubDomainsX*nSubDomainsY) dproc=proc+nSubDomainsX*nSubDomainsY*(nSubDomainsZ-1)          
#endif                		 
	 end subroutine InitNeighborsMPI
	 
     subroutine AllocateFldVars
          allocate(Ex(mx,my,mz),Ey(mx,my,mz),Ez(mx,my,mz))
          allocate(Bx(mx,my,mz),By(mx,my,mz),Bz(mx,my,mz))
          allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz))
          allocate(F0(mx,my,mz)) ! dummy field variable used in various places for temp. copy etc. 
          if(nMoverEMfilter.gt.0) then 
               allocate(FilteredEx(mx,my,mz),FilteredEy(mx,my,mz),FilteredEz(mx,my,mz))
			   FilteredEx=0;FilteredEy=0;FilteredEz=0;
          end if 		  
     end subroutine AllocateFldVars
	 
	 subroutine ReshapeAuxFld
	     deallocate(Jx,Jy,Jz,F0)
	     allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz),F0(mx,my,mz))
         call InitCurrent
	     
		 deallocate(buff_tJx,buff_tJy,buff_tJz,buff_bJx,buff_bJy,buff_bJz)
	     allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz)) 	 
#ifndef twoD
	     deallocate(buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz)
		 allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))    
#endif  		 
	 end subroutine ReshapeAuxFld
	 
	 subroutine InitMoveDeposit
		 
#ifdef twoD	
          allocate(VecEx(4,mx*my),VecEy(4,mx*my),VecEz(4,mx*my))
          allocate(VecBx(4,mx*my),VecBy(4,mx*my),VecBz(4,mx*my))
          allocate(VecJ(8,mx*my,Nthreads))
#else
          allocate(VecEx(8,mx*my*mz),VecEy(8,mx*my*mz),VecEz(8,mx*my*mz))	
          allocate(VecBx(8,mx*my*mz),VecBy(8,mx*my*mz),VecBz(8,mx*my*mz))	
          allocate(VecJ(12,mx*my*mz,Nthreads))
		  	
#endif
			VecEx=0
			VecEy=0
			VecEz=0
			VecBx=0
			VecBy=0
			VecBz=0
			
	 end subroutine InitMoveDeposit 
	 
	 subroutine ReshapeShortMoverFldArr
	     !Now update  arrays used in mover 
	      deallocate(VecEx,VecEy,VecEz,VecBx,VecBy,VecBz,VecJ)
		  
		  call InitMoveDeposit

	 end subroutine ReshapeShortMoverFldArr

     subroutine AllocatePrtlVars
		  Nflvr=2 !by default the number of flv is 2, it is increased on demand in SetQbyM(help_setup.F90)
          allocate(flvrqm(Nflvr),FlvrCharge(Nflvr),FlvrSaveFldData(Nflvr),FlvrType(Nflvr),FlvrSpare(Nflvr))
          allocate(CurrentTagID(Nflvr),TagCounter(Nflvr))
          !Default values 
          flvrqm(1)=qmi
          flvrqm(2)=qme
		  FlvrCharge(1)=1.0_psn
		  FlvrCharge(2)=-1.0_psn
          FlvrSaveFldData(1)=1
          FlvrSaveFldData(2)=1
          FlvrType(1)=0
          FlvrType(2)=0
          FlvrSpare(1)=0
          FlvrSpare(2)=0   
	 	  CurrentTagID(1)=1+proc*NtagProcLen !initialise the value of current tag 
	 	  CurrentTagID(2)=1+proc*NtagProcLen  
		  TagCounter(1)=1
		  TagCounter(2)=1
	  
          call InitPrtlArr(prtl_arr_size)		  
		  
		  !Allocate test prtl array size 
		  test_prtl_arr_size=1024!Nelc/10
#ifdef gpu
          test_prtl_arr_size=1024! Note: test prtls will be removed  
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
     
     subroutine InitPrtlBoundaries
          xmin=3.0_psn
          xmax=mx-2
          ymin=3.0_psn
          ymax=my-2
#ifdef twoD
          zmin=1.0_psn
          zmax=2.0_psn
#else          
          zmin=3.0_psn
          zmax=mz-2
#endif
          xlen=xmax-xmin
          ylen=ymax-ymin
          zlen=zmax-zmin
		  

     end subroutine InitPrtlBoundaries
	 
	 subroutine InitBCpos
		  !BC varaibles
	      BC_Xmin_Prtl=-1; BC_Xmax_Prtl=-1; BC_Ymin_Prtl=-1; BC_Ymax_Prtl=-1; BC_Zmin_Prtl=-1; BC_Zmax_Prtl=-1;
	 	  BC_Xmin_Fld=-1; BC_Xmax_Fld=-1; BC_Ymin_Fld=-1; BC_Ymax_Fld=-1; BC_Zmin_Fld=-1; BC_Zmax_Fld=-1;
	 end subroutine InitBCpos

     subroutine InitDomainSize	 
          mx=xborders(procxind(proc)+1)-xborders(procxind(proc))+5
          my=yborders(procyind(proc)+1)-yborders(procyind(proc))+5
#ifdef twoD
          mz=1 
#else
          mz=zborders(proczind(proc)+1)-zborders(proczind(proc))+5
#endif     
     end subroutine InitDomainSize
! By default the intial elemectromagnetic fields is set to 0, but they can be reset in the setup file  
subroutine InitFlds_default
     integer :: i,j,k
     Ex=0.0_psn; Ey=0.0_psn; Ez=0.0_psn; Bx=0.0_psn; By=0.0_psn; Bz=0.0_psn
     Bx_ext0=0.0_psn; By_ext0=0.0_psn; Bz_ext0=0.0_psn
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
	 real(psn) :: r1
	 real(psn) :: tag_fraction
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
	 
	 !Tag ions and electrons by default
	 CurrentTagID(1)=1+proc*NtagProcLen !initialise the value of current tag
	 CurrentTagID(2)=1+proc*NtagProcLen
	 
	 tag_fraction=1.0_psn/psave_ratio
     if(psave_ratio.gt.0) then !tag the flv=1,2 particles by default
          do n=1,prtl_arr_size
               if(qp(n).ne.0) then
				   call random_number(r1)
				   if(r1.lt.tag_fraction) then
					    if(flvp(n).eq.1) then   
				            tagp(n)=CurrentTagID(1)
				            if(mod(CurrentTagID(1),NtagProcLen).eq.0) CurrentTagID(1)=CurrentTagID(1)+NtagProcLen*(nproc-1)
			                CurrentTagID(1)=CurrentTagID(1)+1
						end if 
					    if(flvp(n).eq.2) then   
				            tagp(n)=CurrentTagID(2)
				            if(mod(CurrentTagID(2),NtagProcLen).eq.0) CurrentTagID(2)=CurrentTagID(2)+NtagProcLen*(nproc-1)
			                CurrentTagID(2)=CurrentTagID(2)+1
						end if 						
			       end if 
		       end if
          end do
     end if

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

subroutine InitPrtlTransferInOutArr
	integer :: n
	loutp_size=0; routp_size=0; toutp_size=0; boutp_size=0; uoutp_size=0; doutp_size=0;
	do n=1,prtl_arr_size
		if(qp(n).ne.0) then
			if(xp(n).lt.xmin+1) loutp_size=loutp_size+1
			if(xp(n).gt.xmax-1) routp_size=routp_size+1
			if(yp(n).lt.ymin+1) boutp_size=boutp_size+1
			if(yp(n).gt.ymax-1) toutp_size=toutp_size+1
#ifndef twoD			
			if(zp(n).lt.zmin+1) doutp_size=doutp_size+1
			if(zp(n).gt.zmax-1) uoutp_size=uoutp_size+1
#endif			
		end if 
	end do 
	
	do n=1,test_prtl_arr_size
		if(qtp(n).ne.0) then
			if(xtp(n).lt.xmin+1) loutp_size=loutp_size+1
			if(xtp(n).gt.xmax-1) routp_size=routp_size+1
			if(ytp(n).lt.ymin+1) boutp_size=boutp_size+1
			if(ytp(n).gt.ymax-1) toutp_size=toutp_size+1
#ifndef twoD			
			if(ztp(n).lt.zmin+1) doutp_size=doutp_size+1
			if(ztp(n).gt.zmax-1) uoutp_size=uoutp_size+1
#endif			
		end if 
	end do 
	
	loutp_size=loutp_size+outp_arr_block_size; routp_size=routp_size+outp_arr_block_size; toutp_size=toutp_size+outp_arr_block_size; boutp_size=boutp_size+outp_arr_block_size; uoutp_size=uoutp_size+outp_arr_block_size; doutp_size=doutp_size+outp_arr_block_size;
#ifdef twoD	
    allocate(loutp(loutp_size),routp(routp_size),toutp(toutp_size),boutp(boutp_size))
#else
    allocate(loutp(loutp_size),routp(routp_size),toutp(toutp_size),boutp(boutp_size),uoutp(uoutp_size),doutp(doutp_size)) 
#endif	 

    linp_size=0; rinp_size=0; tinp_size=0; binp_size=0; uinp_size=0; dinp_size=0           
#ifdef twoD	
    allocate(linp(linp_size),rinp(rinp_size),tinp(tinp_size),binp(binp_size))
#else
    allocate(linp(linp_size),rinp(rinp_size),tinp(tinp_size),binp(binp_size),uinp(uinp_size),dinp(dinp_size)) 
#endif	
end subroutine InitPrtlTransferInOutArr

!---------------------------------------------------------
! Enforce the boundary conditions defined in the setup file
!---------------------------------------------------------
subroutine ComplyBC 
		  integer :: n 
		  real(psn) :: xmax_global,xmin_global,ymin_global,ymax_global,zmin_global,zmax_global
		  real(psn) :: xmax_local,xmin_local,ymin_local,ymax_local,zmin_local,zmax_local

		  !----------  Particles --------------!
		  
		  xmax_global=BC_Xmax_Prtl
		  xmin_global=BC_Xmin_Prtl
		  ymax_global=BC_Ymax_Prtl
		  ymin_global=BC_Ymin_Prtl
		  zmax_global=BC_Zmax_Prtl
		  zmin_global=BC_Zmin_Prtl
		  
		  if(BC_Xmax_Prtl.lt.0) xmax_global=xborders(nSubDomainsX)+1
		  if(BC_Xmin_Prtl.lt.0) xmin_global=xborders(0)-1
		  if(BC_Ymax_Prtl.lt.0) ymax_global=yborders(nSubDomainsY)+1
		  if(BC_Ymin_Prtl.lt.0) ymin_global=yborders(0)-1
		  xmax_local=xmax_global-xborders(procxind(proc))+3
		  xmin_local=xmin_global-xborders(procxind(proc))+3
		  ymax_local=ymax_global-yborders(procyind(proc))+3
		  ymin_local=ymin_global-yborders(procyind(proc))+3		  		  
#ifdef twoD		  
          zmax_local=2
          zmin_local=-1
#else 	  
		  if(BC_Zmax_Prtl.lt.0) zmax_global=zborders(nSubDomainsZ)+1
		  if(BC_Zmin_Prtl.lt.0) zmin_global=zborders(0)-1
		  zmax_local=zmax_global-zborders(proczind(proc))+3
		  zmin_local=zmin_global-zborders(proczind(proc))+3
#endif 		  


          do n=1,prtl_arr_size
               if(qp(n).ne.0) then
				   if((xp(n).lt.xmin_local).or.(xp(n).gt.xmax_local)) call DeletePrtl(n)
				   if((yp(n).lt.ymin_local).or.(yp(n).gt.ymax_local)) call DeletePrtl(n)
				   if((zp(n).lt.zmin_local).or.(zp(n).gt.zmax_local)) call DeletePrtl(n)		 
			   end if
		  end do 
		  
          do n=1,test_prtl_arr_size
               if(qtp(n).ne.0) then
				   if((xtp(n).lt.xmin_local).or.(xtp(n).gt.xmax_local)) call DeleteTestPrtl(n)
				   if((ytp(n).lt.ymin_local).or.(ytp(n).gt.ymax_local)) call DeleteTestPrtl(n)
				   if((ztp(n).lt.zmin_local).or.(ztp(n).gt.zmax_local)) call DeleteTestPrtl(n)	
			   end if
		  end do 
		  
		  !----------  Fields --------------!
		  if(BC_Xmin_Fld_Type.eq.'cond') call SetFldPairToZero(Ey,Ez,1, int(BC_Xmin_Fld)-xborders(procxind(proc))+3, 1,my, 1,mz)
		  
end subroutine ComplyBC

subroutine SetFldPairToZero(Fld1,Fld2,i1,i2,j1,j2,k1,k2)
	real(psn), dimension(mx,my,mz) :: Fld1,Fld2
	integer :: i1,i2,j1,j2,k1,k2
	integer :: i,j,k
	
	if(i1.gt.mx .or. i2.lt.1) return
	if(j1.gt.my .or. j2.lt.1) return
	if(k1.gt.mz .or. k2.lt.1) return
	
	do k = max(k1,1),min(k2,mz)
		do j = max(j1,1),min(j2,my)
			do i = max(i1,1),min(i2,mx)
				Fld1(i,j,k) = 0.0_psn
				Fld2(i,j,k) = 0.0_psn
			end do 
		end do 
	end do
end subroutine SetFldPairToZero




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