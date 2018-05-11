module setup
     use parameters
     use vars
	 !use fields, only: smoothen_current_subdomain
     use help_setup
	 use memory
	 use savedata
	 use loadbalance
	 use deposit
     implicit none 
	 logical :: ParticleSplitting,InjectTurb,RadiationBoundaryCondition,TrimUpstream,SlowDownInjector,InitReflected,InjTestPrtl,TwoStream,CleanUpstreamPrtl,InitThermalBath
	 integer :: BoxSizeIncrement,current_filter
	 integer, parameter :: GammaTableLength=10000
	 integer, parameter :: InjCountTableSize=1000
 	 integer, parameter :: RefCells=c_ompe*8*8
	 real(dbpsn) :: xinj
	 real(psn)   :: time_this
     real(psn):: DriftUpstream,dxUpstream,dxInjector,TrimSpeed,SlowDownSpeed,PUIratio
	 integer  :: Xshock,xcond,delX_PrtlSplit,PrtlWall,InjBuffZone,RadiationBoundary,TrimPeriod,TrimAfterTimeStep,SlowDownAfter,TPratio,TPtag_ratio,SlowDownRefPrtl,StopWallScatter
     integer :: alast=1,untagged_pairs=1
	 real(psn):: Temp_ion, Temp_elc, TeTiRatio
	 real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion,GammaPUI,Gamma_PDF_PUI
	 !real(psn), dimension(:), allocatable :: InjCountTable1,InjCountTable2,InjCountTable3

	 integer, parameter :: Nmodes=7 ! number of Alfven modes 
	 real(psn), dimension(Nmodes,3):: vecK, vecKXB, vecKXBXB, vecKXdB
	 real(psn), dimension(Nmodes)  :: AmpK, PhaseK 
	 
	 real(psn) :: Alfven_speed, PlasmaBeta, Btheta,Bphi,BextMag
	 real(psn) :: driftEx,driftEy,driftEz
contains 
	
subroutine InitUser
	    
	!------------------------------------------------------------
    !Physical Parameters
	!------------------------------------------------------------
		Alfven_speed=0.01_psn
		PlasmaBeta=0.1
		DriftUpstream=-0.05!-g0
		dxInjector=abs(1.5_psn*DriftUpstream)
        TeTiRatio=1.0
		PUIratio=0.25!fraction of pickup ions 
		TwoStream=.false.
		StopWallScatter=1!particle's momentum is isotropized 1 cell ahead of the wall 
		ParticleSplitting=.false.
		InjectTurb=.false.
		RadiationBoundaryCondition=.false. ! True: Radiation Boundary Condition; false: Expanding Box
		TrimUpstream=.false. !Upstream is cut short by resetting injector's location
		SlowDownInjector=.true.
		SlowDownRefPrtl=1000 ! Particles reflected from the shock that reach the injector are slowed down after this many time steps
		InitReflected=.false.
		InjTestPrtl=.false.
		CleanUpstreamPrtl=.false. !very useful for parallel shocks, the routine needs to be set up for different problems
		InitThermalBath=.false.
		

		Btheta=90*(pi/180.0_psn) !Angle (in radians) between magnetic field and z-axis
		Bphi= 80*(pi/180.0_psn)! Angle between x-y projection of magnetic field and x-axis   
		
		!------------------------------------------------------------
	    !Setup Parameters (all positions are in global grid-cell cordiantes)
		!------------------------------------------------------------	
		xcond=4 !global x-cord of the conductor 
		PrtlWall=xcond+6! global x-cord of the Wall that reflects particles 
	    xinj=PrtlWall+32!+128!intial position, in global cordinates 
		BoxSizeIncrement=128
		current_filter=0!number of times current is filtered (edges are treated differently w.r.t. periodic box )
	    delX_PrtlSplit=100
		InjBuffZone=128
		RadiationBoundary=xborders(nSubDomainsX)-10!xinj+64
		TrimPeriod=1300!0	
		TrimAfterTimeStep=10000
		TrimSpeed=abs(DriftUpstream*0.7)
		SlowDownAfter=10000
		SlowDownSpeed=abs(0.5*DriftUpstream)!Injector speed is reduced to this speed
		TPratio=128
		TPtag_ratio=128
		if(TwoStream) xinj=34! 
		time_this=t
			
		!============= Derived Quantitites =================================================
		if(abs(DriftUpstream).gt.1.0_psn) then 
			dxUpstream=c*sqrt((DriftUpstream-1.0_psn)*(DriftUpstream+1.0_psn))/abs(DriftUpstream)
		else 
			dxUpstream=c*abs(DriftUpstream)
		end if
		
        Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
        Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
        
		call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
        call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
		
		if(PUIratio.gt.0) then
			call InitNewFlvr(FlvID=3,QbyM=qmi,Type=0,SaveFld=1,Split=0)
			!call Init_KappaDist_Table(GammaTableLength,GammaPUI,Gamma_PDF_PUI,Temp_ion*100,4)
			call Init_PUI_ShellDist_Table(GammaTableLength,GammaPUI,Gamma_PDF_PUI,2.0*abs(DriftUpstream))
		end if 
		
		if(InitThermalBath) call InsertThermalBathPrtl
		!------------------------------------------------------------
	    !Derived quantities
		!------------------------------------------------------------
        !set the background magentic field 
        BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
        Bz_ext0=BextMag*cos(Btheta)
        Bx_ext0=BextMag*sin(Btheta)*cos(Bphi)
		By_ext0=BextMag*sin(Btheta)*sin(Bphi)
		
        driftEx=0.0
        driftEy=DriftUpstream*Bz_ext0
        driftEz=-DriftUpstream*By_ext0
		
        if((curr_filter.gt.0).and.(current_filter.gt.0)) print*,'Warning!! Two different current filtering schemes are used simultaneously'  
		!There is a reflecting wall at the very left 
		if(.not.restart) then 
			if(TwoStream) then 
				 xborders=xborders-nx/2
	 			 call SetUpstreamFldRight(0-xborders(procxind(proc))+3)
				 call SetUpstreamFldLeft(0-xborders(procxind(proc))+3)
			else 
				call SetUpstreamFldRight(PrtlWall-xborders(procxind(proc))+3)
			end if  
	
			if(InjTestPrtl) call InitTestPrtl
		else			
            call ReadRestartDataSetup
		end if  
		
end subroutine InitUser


subroutine InjectNewPrtl
	 implicit none 
	 integer  :: procxind_this
     real(psn):: x1,x2,xinj_local
	 integer  :: xrad_local
	 
	 procxind_this=procxind(proc)
	 
	 !----First Inject Left Moving Particles ----
	 
	 !x1 and x2 determined boundary of the region in local cordinates where new upstream plasma is injected 
	 x1=max(xinj,real(xborders(procxind_this)))
	 x2=min(xinj+dxInjector,real(xborders(procxind_this+1)))
	 
	 x1=x1-xborders(procxind_this)+3 !change to local cordinate 
	 x2=x2-xborders(procxind_this)+3 
	 
	 xinj_local=max(0.5_psn,xinj-xborders(procxind_this)+3)
	 xrad_local=RadiationBoundary-xborders(procxind(proc))+3
	 !xrad_local=max(1,floor(RadiationBoundary)-xborders(procxind_this)+3)
	 
	 !call ClearInjectionRegionRight(xinj_local+1)
	 call ReflectNewPrtlRight(xinj_local)  ! Is it needed if I am using steady injector?
	 
	 !call AttenuateFld(xrad_local)
! 	 if(RadiationBoundaryCondition) then
! 		 call SetUpstreamFldRight(xrad_local+1) !only the ghost zones should be reset
! 	 else
! 		 xrad_local=xborders(nSubDomainsX)-8-xborders(procxind(proc))+3
! 		 call SetUpstreamFldRight(xrad_local)
! 	 end if
	 
	 
	  call HomogeneousDriftingPrtlEqualWeightPUI_SteadyInjector(DriftUpstream,abs(3.0*DriftUpstream)) 
	  call HomogeneousDriftingPrtlEqualWeightPUI(x1,x2,DriftUpstream) 
	 !call HomogeneousDriftingPrtl_SteadyInjector(DriftUpstream,abs(3*DriftUpstream)) 
	 !call HomogeneousDriftingPrtl(x1,x2,DriftUpstream) 
	 
     if(InjTestPrtl) call InjectTestPrtl(x1,x2,DriftUpstream) 
	 
	 !----Now Right Moving Particles ----
	 if(TwoStream) then 
		 
		 x1=max(-xinj-dxInjector,real(xborders(procxind_this)))
		 x2=min(-xinj,real(xborders(procxind_this+1)))
	 
		 x1=x1-xborders(procxind_this)+3 !change to local cordinate 
		 x2=x2-xborders(procxind_this)+3 
	 
		 xinj_local=min(mx-0.5,-xinj-xborders(procxind_this)+3)
		 xrad_local=-RadiationBoundary-xborders(procxind(proc))+3 !At the moment there is no radiation bounary at the left side 
	 
		 !call ClearInjectionRegionLeft(xinj_local-1)
		 !call ReflectNewPrtlLeft(xinj_local)
		 
! 		 if(RadiationBoundaryCondition) then
! 			 call SetUpstreamFldLeft(xrad_local-1) !only the ghost zones should be reset
! 		 else
! 			 xrad_local=xborders(0)+8-xborders(procxind(proc))+3
! 			 call SetUpstreamFldLeft(xrad_local)
! 		 end if
	 
		 !	 if(InjectTurb) call InitEMfieldAlfvenWave(x1) !At the moment turublence is not inject in this beam 
		 call HomogeneousDriftingPrtlEqualWeightPUI_SteadyInjector(-DriftUpstream,abs(3.0*DriftUpstream)) 
		 call HomogeneousDriftingPrtlEqualWeightPUI(x1,x2,-DriftUpstream) 
		 !call HomogeneousDriftingPrtl_SteadyInjector(-DriftUpstream,abs(3*DriftUpstream)) 
		 !call HomogeneousDriftingPrtl(x1,x2,-DriftUpstream) 
	     if(InjTestPrtl) call InjectTestPrtl(x1,x2,-DriftUpstream) 
     end if 
	 	 
	 
	 xinj=xinj+dxInjector ! update location of the injector 
	 !RadiationBoundary=RadiationBoundary+dxInjector
	 if(proc.eq.0) print*,'xinj',xinj,'xrad',RadiationBoundary
	 		  
end subroutine InjectNewPrtl

!This subsourtine creates all particles of same weight
subroutine HomogeneousDriftingPrtlEqualWeightPUI(x1,x2,drift)
	implicit none 
	real(psn) :: x1,x2,drift
	real(psn) :: dx
    real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,r1,tag_fraction
    real(psn) :: ubeta,vbeta,wbeta,beta_drift
    real(psn) :: MeanCount1,MeanCount2   
	real(psn):: MeanCount1Prev=-1,MeanCount2Prev=-1
	real(psn), dimension(InjCountTableSize) :: InjCountTable1,InjCountTable2
	integer   :: Nelc_New,i,tag
	
	save MeanCount1Prev,MeanCount2Prev,InjCountTable1,InjCountTable2
	
	tag_fraction=1.0_psn/psave_ratio
	CurrentTagID=1+NtagProcLen*proc
    TagCounter=1
	
	
	dx=max(0.0_psn,x2-x1)  
#ifdef twoD
   	MeanCount1=(1.0_psn-PUIratio)*dx*epc*(my-5)
	MeanCount2=PUIratio*dx*epc*(my-5)
#else 
	MeanCount1=(1.0_psn-PUIratio)*dx*epc*(my-5)*(mz-5)
	MeanCount2=PUIratio*dx*epc*(my-5)*(mz-5)
#endif 	

	if(MeanCount1Prev.ne.MeanCount1) then 	
		call InitPoissonDist(InjCountTableSize,InjCountTable1,MeanCount1)
		MeanCount1Prev=MeanCount1
	end if 
 
     call GetInjPrtlCount(MeanCount1,Nelc_New,InjCountTableSize,InjCountTable1)

	 if(Nelc_New.gt.0) then
		 do i=1,Nelc_New
			!Electrons  
            call GenerateTag(tag,1,tag_fraction)
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)			
		    call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,time_this)
			!ions
			call GenerateTag(tag,2,tag_fraction)
			call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)
		 end do
	 end if
     
 	 if(MeanCount2Prev.ne.MeanCount2) then 	
 		 call InitPoissonDist(InjCountTableSize,InjCountTable2,MeanCount2)
 		 MeanCount2Prev=MeanCount2
 	 end if
     
     call GetInjPrtlCount(MeanCount2,Nelc_New,InjCountTableSize,InjCountTable2)	
	  
	 if(Nelc_New.gt.0) then
		 do i=1,Nelc_New
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)			
			!ions
			call GenerateTag(tag,2,tag_fraction)
			call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)

			call GenerateTag(tag,3,tag_fraction*10)
			call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaPUI,Gamma_PDF_PUI)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,3,time_this) 
		 end do
	 end if
end subroutine HomogeneousDriftingPrtlEqualWeightPUI


! The following subroutines assumes a stationary injector and lets particles cross into the simulation box 
subroutine HomogeneousDriftingPrtlEqualWeightPUI_SteadyInjector(drift,max_speed)
	implicit none 
	real(psn) :: drift,max_speed
	real(psn) :: x1,x2,xinj1,xinj2!x1,x2 define domain when plasma is actually injected
	real(psn) :: dx
    real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,r1,tag_fraction,ginv
    real(psn) :: ubeta,vbeta,wbeta,beta_drift
    real(psn) :: MeanCount1,MeanCount2   
	real(psn) :: MeanCount1Prev=-1,MeanCount2Prev=-1
	real(psn), dimension(InjCountTableSize) :: InjCountTable1,InjCountTable2
	integer   :: Nelc_New,i,tag
	integer   :: procxind_this
	
	save MeanCount1Prev,MeanCount2Prev,InjCountTable1,InjCountTable2
	
	tag_fraction=1.0_psn/psave_ratio
	CurrentTagID=1+NtagProcLen*proc
    TagCounter=1
	
    procxind_this=procxind(proc)
	dx=c*max_speed
	if(drift.lt.0) then 
		x1=max(xinj-dx,real(xborders(procxind_this)))
   	    x2=min(xinj,real(xborders(procxind_this+1)))
		xinj1=xinj !domain of the plasma in the upstream which can cross the injector
		xinj2=xinj+dx
	else 
	    x1=max(-xinj,real(xborders(procxind_this))) ! the second beam drifting to the right
	    x2=min(-xinj+dx,real(xborders(procxind_this+1)))
	    xinj1=-xinj-dx
		xinj2=-xinj
	end if 
    
	if(x1.gt.x2) return !Nothing to inject  	
    x1=x1-xborders(procxind_this)+3 !change to local cordinate 
    x2=x2-xborders(procxind_this)+3 	
    xinj1=xinj1-xborders(procxind_this)+3 !change to local cordinate 
    xinj2=xinj2-xborders(procxind_this)+3 	
	
	!print*,'proc',x1,x2,xinj1,xinj2,xinj,t
	  
#ifdef twoD
   	MeanCount1=(1.0_psn-PUIratio)*dx*epc*(my-5)
	MeanCount2=PUIratio*dx*epc*(my-5)
#else
	MeanCount1=(1.0_psn-PUIratio)*dx*epc*(my-5)*(mz-5)
	MeanCount2=PUIratio*dx*epc*(my-5)*(mz-5)
#endif 	

	   if(MeanCount1Prev.ne.MeanCount1) then 	
		call InitPoissonDist(InjCountTableSize,InjCountTable1,MeanCount1)
		MeanCount1Prev=MeanCount1
	    end if 

 
     call GetInjPrtlCount(MeanCount1,Nelc_New,InjCountTableSize,InjCountTable1)
	 if(Nelc_New.gt.0) then
		 do i=1,Nelc_New
			!Ions
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,xinj1,xinj2,ymin,ymax,zmin,zmax)					
		    call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
			
			ginv=1.0_psn/sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
			xlocal=xlocal+ugamma*ginv*c !move particles and see if it happend to be in the simualtion box
  
            if(xlocal.gt.x1.and.xlocal.lt.x2) then     
	            call GenerateTag(tag,1,tag_fraction)	
				call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,time_this)
				!Electrons
				call GenerateTag(tag,2,tag_fraction)
				call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
				call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)
		    end if 
		 end do
	 end if
     
 	 if(MeanCount2Prev.ne.MeanCount2) then 	
 		 call InitPoissonDist(InjCountTableSize,InjCountTable2,MeanCount2)
 		 MeanCount2Prev=MeanCount2
	 end if
     
     call GetInjPrtlCount(MeanCount2,Nelc_New,InjCountTableSize,InjCountTable2)	
	  
	 if(Nelc_New.gt.0) then
		 do i=1,Nelc_New
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,xinj1,xinj2,ymin,ymax,zmin,zmax)		
			call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaPUI,Gamma_PDF_PUI)
			
			ginv=1.0_psn/sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
			xlocal=xlocal+ugamma*ginv*c
	        if(xlocal.gt.x1.and.xlocal.lt.x2) then     
			
				!pick up ions 
				call GenerateTag(tag,3,tag_fraction*10)
				call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,3,time_this) 
			
			
				!electrons
				call GenerateTag(tag,2,tag_fraction)
				call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
				call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)
		    end if
		 end do
	 end if
end subroutine HomogeneousDriftingPrtlEqualWeightPUI_SteadyInjector


subroutine GenerateTag(tag,FlvID,fraction)
	integer :: tag, FlvID 
	real(psn) :: fraction
	real(psn) :: r1
	call random_number(r1)
	if(r1.lt.fraction) then 
		call GetTagID(tag,FlvID,1)
	else 
		tag=0
	end if 	
end subroutine GenerateTag

subroutine InitTestPrtl
	call InitNewFlvr(FlvID=4,QbyM=qmi,Type=-11,SaveFld=1,Split=0)
	call InitNewFlvr(FlvID=5,QbyM=qmi,Type=-11,SaveFld=1,Split=0)
end subroutine InitTestPrtl
subroutine InjectTestPrtl(x1,x2,drift)
	implicit none 
	real(psn) :: x1,x2,drift
    real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,tag_fraction
	real(psn) :: dx,MeanCount3,MeanCount3Prev=0
	real(psn), dimension(InjCountTableSize) :: InjCountTable3
	integer   :: Nelc_New,i,tag
    
	tag_fraction=1.0_psn/TPtag_ratio
	dx=max(0.0_psn,x2-x1)  
#ifdef twoD
	MeanCount3=dx*epc*(my-5)/real(TPratio)
#else 
	MeanCount3=dx*epc*(my-5)*(mz-5)/real(TPratio)
#endif 

	if(MeanCount3Prev.ne.MeanCount3) then 	
		call InitPoissonDist(InjCountTableSize,InjCountTable3,MeanCount3)
		MeanCount3Prev=MeanCount3
	end if 
	
    call GetInjPrtlCount(MeanCount3,Nelc_New,InjCountTableSize,InjCountTable3)
    if(Nelc_New.gt.0) then
		do i=1,Nelc_New
            call GenerateTag(tag,4,tag_fraction)
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)			
		    call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,4,time_this)
			
            call GenerateTag(tag,5,tag_fraction)
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)			
		    call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaPUI,Gamma_PDF_PUI)
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,5,time_this) 		
		end do 
	end if 		
end subroutine InjectTestPrtl

subroutine GetInjPrtlCount(MeanCount,Count,TableSize,Table)
	real(psn) :: MeanCount
	integer   :: Count,TableSize
	real(psn), dimension(TableSize) :: Table 
	if(MeanCount.eq.0) then
		Count=0
		return
	end if 
	Count=nint(MeanCount)
	if(Count.lt.InjCountTableSize) call GetIntPoissonDist(TableSize,Table,Count) ! Possion distribution  		  
end subroutine GetInjPrtlCount


subroutine ClearInjectionRegionRight(x1)
	real(psn) :: x1
	integer :: n
	if(x1.le.mx) then 
		do n=1,used_prtl_arr_size !reflect particles from the injector, relection is done in the lab frame
			if((qp(n).ne.0).and.(xp(n).ge.x1)) up(n)=-up(n)	
		end do 
		do n=1,used_test_prtl_arr_size !reflect particles from the injector,maybe reflecton should be done in upstream frame
			if((qtp(n).ne.0).and.(xtp(n).ge.x1)) utp(n)=-utp(n)
		end do	
	end if
end subroutine ClearInjectionRegionRight

subroutine ClearInjectionRegionLeft(x1)
	real(psn) :: x1
	integer :: n
	if(x1.ge.1) then 
		do n=1,used_prtl_arr_size !reflect particles
			if((qp(n).ne.0).and.(xp(n).le.x1)) up(n)=-up(n)
		end do 
		do n=1,used_test_prtl_arr_size !reflect particles
			if((qtp(n).ne.0).and.(xtp(n).le.x1)) utp(n)=-utp(n)
		end do 
	end if
end subroutine ClearInjectionRegionLeft

subroutine ReflectNewPrtlLeft(x1)
	real(psn) :: x1
	integer :: n
	if(x1.ge.1) then 
		do n=1,used_prtl_arr_size !reflect any particle to the left of the injector 
			if((qp(n).ne.0).and.(xp(n).le.x1)) call ReflectPrtlMovingFrameX(-DriftUpstream,up(n),vp(n),wp(n))
		end do 
		do n=1,used_test_prtl_arr_size !reflect any particle to the left of the injector 
			if((qtp(n).ne.0).and.(xtp(n).le.x1)) call ReflectPrtlMovingFrameX(-DriftUpstream,utp(n),vtp(n),wtp(n))
		end do 
	end if
end subroutine ReflectNewPrtlLeft

subroutine ReflectNewPrtlRight(x1)
	real(psn) :: x1
	real(psn) :: x0,y0,z0
	integer :: n
	if(x1.le.mx) then 
		do n=1,used_prtl_arr_size !remove any particle to the right of the injector 
			if((qp(n).ne.0).and.(xp(n).ge.x1)) then 
				call ReflectPrtlMovingFrameX(DriftUpstream,up(n),vp(n),wp(n))	
				!x0=xp(n)
				!y0=yp(n)
				!z0=zp(n)
				!xp(n)=2*x1-x0
				!call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
			end if 		
		end do 
		do n=1,used_test_prtl_arr_size !remove any particle to the right of the injector 
			if((qtp(n).ne.0).and.(xtp(n).ge.x1)) then 
				xtp(n)=2*x1-xtp(n)
				call ReflectPrtlMovingFrameX(DriftUpstream,utp(n),vtp(n),wtp(n))
			end if 			
		end do 
	end if
end subroutine ReflectNewPrtlRight

subroutine EnlargeBox
	 !Enlarging right side of the box  
	 if((xinj+InjBuffZone).gt.xborders(nSubDomainsX)) then 
		!print*,'Increasing Size of the Box',t
	    if(procxind(proc).eq.nSubDomainsX-1) then
            call EnlargeAllFldArr(BoxSizeIncrement)
            call ReshapeAuxVars
			call SetUpstreamFldRight(mx-BoxSizeIncrement)
	        call ReshapeShortMoverFldArr
	    end if
   	    xborders(nSubDomainsX)=xborders(nSubDomainsX)+BoxSizeIncrement
		fdataxf_box=xborders(nSubDomainsX)
	 end if
	 
	 if(TwoStream) then !Now Enlarging left side of the box
		 if((-xinj-InjBuffZone).lt.xborders(0)) then 
		    if(procxind(proc).eq.0) then
	            call EnlargeAllFldArr(-BoxSizeIncrement)
	            call ReshapeAuxVars
				call SetUpstreamFldLeft(BoxSizeIncrement+3)
		        call ReshapeShortMoverFldArr
				xp(:)=xp(:)+BoxSizeIncrement
				xtp(:)=xtp(:)+BoxSizeIncrement
		    end if
	   	    xborders(0)=xborders(0)-BoxSizeIncrement
			fdataxi_box=xborders(0)
	     end if
	 end if 
end subroutine EnlargeBox

subroutine EnlargeAllFldArr(increment)
     integer :: increment
     call EnlargeFldArr(Ex,increment)
     call EnlargeFldArr(Ey,increment)
     call EnlargeFldArr(Ez,increment)
     call EnlargeFldArr(Bx,increment)
     call EnlargeFldArr(By,increment)
     call EnlargeFldArr(Bz,increment)
	 if(nMoverEMfilter.gt.0) then 
	     call EnlargeFldArr(filteredEx,increment)
	     call EnlargeFldArr(filteredEy,increment)
	     call EnlargeFldArr(filteredEz,increment)
	 end if
end subroutine EnlargeAllFldArr
subroutine EnlargeFldArr(Fld,increment)
	 integer     :: increment
	 real(psn), dimension(:,:,:), allocatable :: Fld
	 real(psn), dimension(:,:,:),allocatable :: FldTemp
     allocate(FldTemp(mx+abs(increment),my,mz))  !negative increment is interpreted as enlargment to the left side of the box
	 FldTemp=0.0_psn   
	 if(increment.gt.0) then 
		 FldTemp(1:mx,:,:)=Fld(1:mx,:,:)
	 else 
		 FldTemp(abs(increment)+1:mx+abs(increment),:,:)=Fld(1:mx,:,:)
	 end if 
	 deallocate(Fld)
     call move_alloc(FldTemp,Fld)
end subroutine EnlargeFldArr
subroutine ReshapeAuxVars
     !update all other relevant auxiliary arrays and varaiables local to this domain 
	 mx=mx+BoxSizeIncrement
     xmax=mx-2
     xlen=xmax-3

     deallocate(Jx,Jy,Jz,F0)
     allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz),F0(mx,my,mz))
     Jx=0.0_psn!reset the current 
     Jy=0.0_psn
     Jz=0.0_psn 
     deallocate(buff_tJx,buff_tJy,buff_tJz,buff_bJx,buff_bJy,buff_bJz)
     allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz)) 	 
#ifndef twoD
     deallocate(buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz)
	 allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))    
#endif      
end subroutine ReshapeAuxVars

subroutine ConductingWall
    integer :: n
	integer :: j,k,xcond_local,PrtlWall_local
	real(psn) :: r1,r2
	real(psn) :: mag,cos_theta,sin_theta,cos_phi,sin_phi,phi,mag_max
	
	xcond_local=xcond-xborders(procxind(proc))+3
	PrtlWall_local=PrtlWall-xborders(procxind(proc))+3
		
	if(PrtlWall_local.lt.0) return 	
	
	
!     if(t.gt.StopWallScatter) then
! 		do n=1,used_prtl_arr_size
! 			if((xp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) up(n)=-up(n)
! 		end do
! 		do n=1,used_test_prtl_arr_size
! 			if((xtp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) utp(n)=-utp(n)
! 		end do
! 	else
		!mag_max=(0.2+t/StopWallScatter)*abs(DriftUpstream)/sqrt((1.0_psn-DriftUpstream)*(DriftUpstream+1.0_psn))
		if(t.gt.StopWallScatter) then
			do n=1,used_prtl_arr_size
				if((xp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) then !particle's wall is few cells ahead of the conductor
					xp(n)=PrtlWall_local+(PrtlWall_local-xp(n))
					up(n)=-up(n)
				end if
			end do
			do n=1,used_test_prtl_arr_size
				if((xtp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) then !particle's wall is few cells ahead of the conductor
					xtp(n)=PrtlWall_local+(PrtlWall_local-xtp(n))
					utp(n)=-utp(n)
				end if
			end do
		else 		
			do n=1,used_prtl_arr_size
				if(xp(n).lt.PrtlWall_local) then 
					call random_number(r1)
					cos_theta=r1
					sin_theta=sin(acos(cos_theta))
					call random_number(r2)
					phi=2*r2-1
					!cos_phi=cos(phi)
					!sin_phi=sin(phi)
					mag=sqrt(up(n)**2+vp(n)**2)
					!mag=min(mag,mag_max)
					up(n)= mag*cos_theta
					vp(n)= mag*sin_theta*sign(1.0,phi)!*cos_phi
					!wp(p)= mag*sin_theta*sin_phi
			    end if 
			end do 
			do n=1,used_test_prtl_arr_size
			    if(xtp(n).lt.PrtlWall_local) then 
					call random_number(r1)
					cos_theta=r1
					sin_theta=sin(acos(cos_theta))
					call random_number(r2)
					phi=2*r2-1
					!cos_phi=cos(phi)
					!sin_phi=sin(phi)
					mag=sqrt(utp(n)**2+vtp(n)**2)
					!mag=min(mag,mag_max)
					utp(n)= mag*cos_theta
					vtp(n)= mag*sin_theta*sign(1.0,phi)!*cos_phi
					!wtp(n)= mag*sin_theta*sin_phi 
			    end if
			end do
		end if 	
	
	if(xcond_local.ge.1) then 
	     Ey(1:min(mx,xcond_local),:,:)=0.0_psn
	     Ez(1:min(mx,xcond_local),:,:)=0.0_psn
	end if

end subroutine ConductingWall


subroutine WallScatterPrtl(x1)
	integer, intent(in) :: x1
	integer :: x1_local,x2_local
	integer :: n 
	real(psn) :: r1,r2
	real(psn) :: mag,cos_theta,sin_theta,cos_phi,sin_phi,phi
	
	x1_local=x1-xborders(procxind(proc))+3
	x2_local=x1+1-xborders(procxind(proc))+3
	if((x1_local.gt.mx-2).or.(x2_local.lt.3)) return 
	
	do n=1,used_prtl_arr_size
		if((xp(n).gt.x1_local).and.(xp(n).lt.x2_local)) then 
			call random_number(r1)
			cos_theta=2*r1-1
			sin_theta=sin(acos(cos_theta))
			call random_number(r2)
			phi=2*pi*r2
			cos_phi=cos(phi)
			sin_phi=sin(phi)
			mag=sqrt(up(n)**2+vp(n)**2+wp(n)**2)
			up(n)= mag*cos_theta
			vp(n)= mag*sin_theta*cos_phi
			wp(n)= mag*sin_theta*sin_phi
	    end if 
	end do 
	do n=1,used_test_prtl_arr_size
	    if((xtp(n).gt.x1_local).and.(xtp(n).lt.x2_local)) then 
			call random_number(r1)
			cos_theta=2*r1-1
			sin_theta=sin(acos(cos_theta))
			call random_number(r2)
			phi=2*pi*r2
			cos_phi=cos(phi)
			sin_phi=sin(phi)
			mag=sqrt(utp(n)**2+vtp(n)**2+wtp(n)**2)
			utp(n)= mag*cos_theta
			vtp(n)= mag*sin_theta*cos_phi
			wtp(n)= mag*sin_theta*sin_phi
	    end if
	end do 	
end subroutine WallScatterPrtl


subroutine AttenuateFld(x2)
	integer :: x1,x2
	integer :: Atten_Ncell=40!number of cells for attenuation
	real(psn) :: AttenFactor=0.9_psn
    integer :: i,j,k,imin,imax
	
	x1=x2-Atten_Ncell 
	imax=min(x2,mx)
	imin=max(1,x1)
    if(x1.le.mx) then
		 do k=1,mz
	          do j=1,my
	               do i=imin,imax                   
	                    Ex(i,j,k)=(Ex(i,j,k)-driftEx)*AttenFactor+driftEx		
	                    Ey(i,j,k)=(Ey(i,j,k)-driftEy)*AttenFactor+driftEy
	                    Ez(i,j,k)=(Ez(i,j,k)-driftEz)*AttenFactor+driftEz   
	                    Bx(i,j,k)=Bx(i,j,k)*AttenFactor		
	                    By(i,j,k)=By(i,j,k)*AttenFactor
	                    Bz(i,j,k)=Bz(i,j,k)*AttenFactor                   
	               end do
	          end do
	     end do
    end if		
end subroutine AttenuateFld

subroutine SetUpstreamFldRight(x1)!to ensure that the fld right of the injector is always uniform drift value 
	integer :: x1
    integer :: i,j,k,imin
	
	imin=max(1,x1)
	if(x1.le.mx) then
		 do k=1,mz
	          do j=1,my
	               do i=imin,mx                   
	                    Ex(i,j,k)=driftEx
	                    Ey(i,j,k)=driftEy
	                    Ez(i,j,k)=driftEz        
	               end do
	          end do
	     end do
	end if 
end subroutine SetUpstreamFldRight

subroutine SetUpstreamFldLeft(x1)!to ensure that the fld right of the injector is always uniform drift value 
	integer :: x1
    integer :: i,j,k,imax	
	imax=min(mx,x1)
	if(x1.ge.1) then
		
		 do k=1,mz
	          do j=1,my
	               do i=1,imax     
	                    Ex(i,j,k)=-driftEx
	                    Ey(i,j,k)=-driftEy
	                    Ez(i,j,k)=-driftEz       
	               end do
	          end do
	     end do
	end if 
end subroutine SetUpstreamFldLeft

!The following two subroutines implements absorbing boundary conditions 
subroutine FldBoundaryRight(x1)
	integer :: x1
	integer, parameter :: Ncell=8 !Number of cells in the dissipating zone 
    integer :: i,j,k,imin
	real(psn) :: AttenFactor
	
	imin=max(1,x1)
	if(x1.le.mx) then
		 do k=1,mz
	          do j=1,my
	               do i=imin,mx  
					    if(i.ge.x1+Ncell) then
		                    Ex(i,j,k)=driftEx
		                    Ey(i,j,k)=driftEy
		                    Ez(i,j,k)=driftEz        
					    else 
						   	AttenFactor=(real(abs(x1-i))/real(Ncell))**3							
	   	                    Ex(i,j,k)=Ex(i,j,k)-(Ex(i,j,k)-driftEx)*AttenFactor
	   	                    Ey(i,j,k)=Ey(i,j,k)-(Ey(i,j,k)-driftEy)*AttenFactor
	   	                    Ez(i,j,k)=Ez(i,j,k)-(Ez(i,j,k)-driftEz)*AttenFactor   
					   end if 
	               end do
	          end do
	     end do
	end if 	
end subroutine FldBoundaryRight

subroutine FldBoundaryLeft(x1)
	integer :: x1
	integer, parameter :: Ncell=8 !Number of cells in the dissipating zone 
    integer :: i,j,k,imax
	real(psn) :: AttenFactor	
	
	imax=min(mx,x1)	
	if(x1.ge.1) then	
		 do k=1,mz
	          do j=1,my
	               do i=1,imax  
				    if(i.le.x1-Ncell) then                				      
	                    Ex(i,j,k)=-driftEx
	                    Ey(i,j,k)=-driftEy
	                    Ez(i,j,k)=-driftEz  
					else 
					    AttenFactor=(real(abs(x1-i))/real(Ncell))**3
   	                    Ex(i,j,k)=Ex(i,j,k)-(Ex(i,j,k)+driftEx)*AttenFactor ! note that the drift is negative here 
   	                    Ey(i,j,k)=Ey(i,j,k)-(Ey(i,j,k)+driftEy)*AttenFactor
   	                    Ez(i,j,k)=Ez(i,j,k)-(Ez(i,j,k)+driftEz)*AttenFactor   
					end if       
	               end do
	          end do
	     end do
	end if 
end subroutine FldBoundaryLeft


subroutine RadiationBoundaryX(ind)
	integer :: ind
	integer :: i,j,k
	real(psn) :: c1,c2,c3
	c1=2.0_psn*fldc/(1.0_psn+fldc)
    c2=0.41421356237_psn ! (1-tan^2(pi/8))/2
	c3=0.5_psn*(1.0_psn-c2)*c1
	
#ifdef twoD
	do k=1,1
		do j=2,my-1
			By(ind,j,k)=By(ind,j,k)&
			            +c1*( By(ind-1,j,k)-By(ind,j,k) + c2*(Bx(ind,j,k)-Bx(ind,j-1,k)))&
						-c3*((Ex(ind,j+1,k)-Ex(ind,j,k)) +(Ex(ind-1,j,k)-Ex(ind-1,j+1,k)) )&
						-fldc*(Ez(ind,j,k)-Ez(ind-1,j,k))
					
			Bz(ind,j,k)=Bz(ind,j,k)&
			            +c1*( Bz(ind-1,j,k)-Bz(ind,j,k) )&
						+fldc*(Ey(ind,j,k)-Ey(ind-1,j,k))-fldc*(Ex(ind-1,j+1,k)-Ex(ind-1,j,k))						
		end do 
	end do
#else	
	do k=2,mz-1
		do j=2,my-1
			By(ind,j,k)=By(ind,j,k)&
			            +c1*( By(ind-1,j,k)-By(ind,j,k) + c2*(Bx(ind,j,k)-Bx(ind,j-1,k)))&
						-c3*((Ex(ind,j+1,k)-Ex(ind,j,k)) +(Ex(ind-1,j,k)-Ex(ind-1,j+1,k)) )&
						+fldc*(Ex(ind-1,j,k+1)-Ex(ind-1,j,k))-fldc*(Ez(ind,j,k)-Ez(ind-1,j,k))
						
			Bz(ind,j,k)=Bz(ind,j,k)&
			            +c1*( Bz(ind-1,j,k)-Bz(ind,j,k) +c2*(Bx(ind,j,k)-Bx(ind,j,k-1)))&
						-c3*( (Ex(ind,j,k+1)-Ex(ind,j,k)) + (Ex(ind-1,j,k+1)-Ex(ind-1,j,k)) )&
						+fldc*(Ey(ind,j,k)-Ey(ind-1,j,k))-fldc*(Ex(ind-1,j+1,k)-Ex(ind-1,j,k))						
		end do 
	end do
#endif   
		
end subroutine RadiationBoundaryX

subroutine LargeBoxInitLoadBalance
	integer :: i,large_domain
	
	large_domain=(nx-32*(nSubDomainsX/2))/real(nSubDomainsX/2)
	if(nx.gt.32*nSubDomainsX) then
		xborders_new(0)=xborders(0)
		do i=1,nSubDomainsX-1
			if(i.le.nSubDomainsX/2) then
				 xborders_new(i)=xborders_new(i-1)+32
			 else
				 xborders_new(i)=xborders_new(i-1)+large_domain
			 end if
		end do
		xborders_new(nSubDomainsX)=xborders(nSubDomainsX)
		call SetNewXBorders
        call ReorderPrtlArr
	    call ReorderTestPrtlArr
    end if
end subroutine LargeBoxInitLoadBalance

subroutine TrimUpstreamLength
	integer :: TrimLen
	real(psn) :: xinj_local
	if(t.lt.TrimAfterTimeStep) return
	if(modulo(t,TrimPeriod).eq.0) then 
	
	    TrimLen=floor(TrimSpeed*c*TrimPeriod)
		xinj=xinj-TrimLen
		xinj_local=max(0.5_psn,xinj-xborders(procxind(proc))+3)
		call ClearInjectionRegionRight(xinj_local)
    end if 
end subroutine TrimUpstreamLength

subroutine SlowInjector
	if(t.lt.SlowDownAfter) return
	dxInjector=SlowDownSpeed+0.99*abs(dxInjector-SlowDownSpeed)
    !call SlowDownParticles !This subroutine does not seem to work with pick up ions  
end subroutine SlowInjector

subroutine SlowDownParticles
	integer, parameter :: DecCells=64!number of cells used to apply deceleration
	real(psn)    :: XstartDec_local,gbeta_drift
	integer      :: n,Tmax
		
	if(abs(DriftUpstream).gt.1.0_psn) then 
		gbeta_drift=-sqrt((DriftUpstream-1.0_psn)*(DriftUpstream+1.0_psn))
	else 
		gbeta_drift=DriftUpstream/sqrt((1.0_psn-DriftUpstream)*(DriftUpstream+1.0_psn))
	end if
	XStartDec_local=xinj-xborders(procxind(proc))+3-DecCells
	Tmax=t-floor(2*real(DecCells)/abs(DriftUpstream)) !assuming a cold drifting plasma
	if(XStartDec_local.gt.1) then
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
		    if((var1p(n).lt.Tmax).and.(xp(n).gt.XStartDec_local)) then ! to make sure that the particle is old enough
				up(n)=gbeta_drift+0.5_psn*(up(n)-gbeta_drift)
                vp(n)=0.5_psn*vp(n)
				wp(n)=0.5_psn*wp(n)
			end if	 
		end do 
		do n=1,used_test_prtl_arr_size
			if(qtp(n).eq.0) cycle
		    if((var1tp(n).lt.Tmax).and.(xtp(n).gt.XStartDec_local)) then ! to make sure that the particle is old enough
				utp(n)=gbeta_drift+0.5_psn*(utp(n)-gbeta_drift)
                vtp(n)=0.5_psn*vtp(n)
				wtp(n)=0.5_psn*wtp(n)
			end if	 
		end do 
	end if 	
end subroutine SlowDownParticles


!The following subroutine is used to intialize a nearly cold gas before the wall 
subroutine InsertThermalBathPrtl
	 real(psn) :: x1,x2  
     real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,TempThis
 	 integer   :: i,Nelc_New,tag
	 
	
	 x1=PrtlWall-xborders(procxind(proc))+3 
	 x2=xinj-12-xborders(procxind(proc))+3
	 x1=max(x1,xmin)
	 x2=min(x2,xmax)
	 
#ifdef twoD 
     Nelc_New=2*nint((x2-x1)*epc*(my-5))!First compute how many new particles are to be added 
#else 	 
     Nelc_New=2*nint((x2-x1)*epc*(my-5)*(mz-5))
#endif

    TempThis=0.1*(DriftUpstream**2) 
    if(Nelc_New.gt.0) then 	
         do i=1,Nelc_New	 
	         
			 !electrons 
		     call GetTagID(tag,1)
			 call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)
			 call GetVelGamma_MaxwellNonRel_StaticPDF(0.0_psn,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
			 call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,time_this)	
    
			 !ions 
			 call GetTagID(tag,2)
			 !call GetVelGamma_MaxwellNonRel(0.0_psn,Temp_elc,ugamma,vgamma,wgamma)
             call GetVelGamma_MaxwellNonRel_StaticPDF(0.0_psn,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
			 call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)
		 end do 
	 end if 	 	 	  
end subroutine InsertThermalBathPrtl

subroutine SaveRestartDataSetup
	if((restart_save_period.gt.0).and.(modulo(t,restart_save_period).eq.0)) then
		if(proc.ne.0) return
	    write(str1,'(I0)') t
        fname=trim(data_folder)//"/restart"//"/SetupVar_"//trim(str1)
        call h5open_f(err)
        call h5fcreate_f(fname,H5F_ACC_TRUNC_F,fid, err)
        rank=1
        data_dim1(1)=1
        call h5screate_simple_f(rank,data_dim1,dspace_id,err)
	    call h5dcreate_f(fid,'xinj',H5T_NATIVE_DOUBLE,dspace_id,dset_id,err)
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,xinj,data_dim1,err)
		call h5dclose_f(dset_id,err)
        call h5sclose_f(dspace_id,err)
        call h5fclose_f(fid,err)
        call h5close_f(err)
    end if 
end subroutine SaveRestartDataSetup
subroutine ReadRestartDataSetup
    write(str1,'(I0)') restart_time
    fname=trim(data_folder)//"/restart"//"/SetupVar_"//trim(str1)
    call h5open_f(err)   
    call h5fopen_f(fname,H5F_ACC_RDONLY_F,fid, err)
    data_dim1(1)=1
    call h5dopen_f(fid,'xinj', dset_id, err)
    call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,xinj,data_dim1,err)
    call h5dclose_f(dset_id,err)
	call h5fclose_f(fid,err)  
    call h5close_f(err)  
end subroutine ReadRestartDataSetup
 

!================================================================================
!The following subroutines are meant for customization and must be part of the setup module,
!even though they are not used and the default subroutines are used instead
!Note:: in order to enable use of the following subroutines the corresponding custom 
!variables must be change to .true.
!================================================================================
subroutine InitOverride !override the default intial conditons
	!xborders(0)=xwall
end subroutine InitOverride  

!The following subroutine is called after moving particles and depositing the current on grid 
! Any additional step which are setup specific should be implelemnted here 
subroutine PostMovDep
	integer :: PrtlWall_local,xrad_local
	real(psn) :: xinj_local
	
	xinj_local=max(0.5_psn,xinj-xborders(procxind(proc))+3)
	!call ReflectNewPrtlRight(xinj_local)
	
	PrtlWall_local=PrtlWall-xborders(procxind(proc))+3
	xrad_local=RadiationBoundary-xborders(procxind(proc))+3
	if(t.lt.StopWallScatter) return
	if(.not.TwoStream) then !Current Adjustments at the Prtl Wall 
		if(PrtlWall_local.ge.3.and.PrtlWall_local.le.mx-2) then
			!To  ensure that the current is deposited on right place for reflected particles
			Jx(PrtlWall_local,:,:)=Jx(PrtlWall_local,:,:)-Jx(PrtlWall_local-1,:,:)
		 	Jy(PrtlWall_local,:,:)=Jy(PrtlWall_local,:,:)+Jy(PrtlWall_local-1,:,:)
		 	Jz(PrtlWall_local,:,:)=Jz(PrtlWall_local,:,:)+Jz(PrtlWall_local-1,:,:)
			Jx(PrtlWall_local-1,:,:)=0.0_psn
			Jy(PrtlWall_local-1,:,:)=0.0_psn
			Jz(PrtlWall_local-1,:,:)=0.0_psn
		end if
    end if 
    
if(RadiationBoundaryCondition.and.(xrad_local.le.mx).and.(xrad_local.ge.2)) call RadiationBoundaryX(xrad_local)
end subroutine PostMovDep

subroutine PreAddCurrent
	integer :: PrtlWall_local,xinj_local,xinj_left_local
	xinj_local=floor(xinj)-xborders(procxind(proc))+3
	PrtlWall_local=PrtlWall-xborders(procxind(proc))+3
!     if(TwoStream) then
! 		xinj_left_local=floor(-xinj)-xborders(procxind(proc))+3
! 		call smoothen_current_subdomain(current_filter,xinj_left_local,xinj_local)
! 	else
! 		call smoothen_current_subdomain(current_filter,PrtlWall_local,xinj_local)
! 	end if 
end subroutine PreAddCurrent

subroutine Finalsubroutines !this subroutine is called at the end of every time step
     !define here all additional actions that is to be perform at the end of every time step
     time_this=t ! =var1 which is used to store injection time of particles 
	 if(.not.TwoStream) call ConductingWall
	 !if(WallScatter) call WallScatterPrtl(PrtlWall+1)
	 
	 call InjectNewPrtl
	 if(TwoStream) then !Load balancing 
		 call LoadBalanceShockTwoStream(xinj)
	 else 
		 call LoadBalanceShock(xinj,PrtlWall+4) !+4 is to ensure that proc. with wall is underloaded for extra work
	 end if 
	 if((t.eq.1).and.(.not.TwoStream)) call LargeBoxInitLoadBalance
	 if(TrimUpstream) call TrimUpstreamLength
	 if(SlowDownInjector) call SlowInjector !Not Implemented for two Streams
	 if((SlowDownRefPrtl.gt.0).and.(t.gt.SlowDownRefPrtl)) call SlowDownParticles
	 !if(CleanUpstreamPrtl) call RemoveWallPrtlFromUpstream
	 
	 !call SetUpstreamFldRight(xborders(nSubDomainsX)-8-xborders(procxind(proc))+3)
	 !if(TwoStream) call SetUpstreamFldLeft(xborders(0)+64-xborders(procxind(proc))+3)
	 
	 call FldBoundaryRight(xborders(nSubDomainsX)-32-xborders(procxind(proc))+3) ! Currently 8 cells are in the boundary layer
	 if(TwoStream) call FldBoundaryLeft(xborders(0)+32-xborders(procxind(proc))+3)
	 
	 call EnlargeBox
	 
	 
	 call SplitParticlesShock
     if(modulo(t,spec_save_period).eq.0) then 
	 	 call FindShockXpos_BmagPeak(Xshock)
		 if(g0.gt.1) then 
		     call SaveGammaSpecInSubDomain(Xshock-delX_PrtlSplit,Xshock+delX_PrtlSplit,&
		     yborders(0),yborders(nSubDomainsY),zborders(0),zborders(size(zborders)-1),'1') ! eng. spec at the shock	 
		 else 
		     call SaveSpeedSpecInSubDomain(Xshock-delX_PrtlSplit,Xshock+delX_PrtlSplit,&
		     yborders(0),yborders(nSubDomainsY),zborders(0),zborders(size(zborders)-1),'1') ! eng. spec at the shock 
			 
		     call SaveSpeedSpecInSubDomain(PrtlWall,PrtlWall+100,&
		     yborders(0),yborders(nSubDomainsY),zborders(0),zborders(size(zborders)-1),'2') ! to save energy data	 
		 end if 
     end if 
	 
	 !if(t.gt.100000) prtl_save_period=800
	 !if(t.gt.100000) prtl_save_period=6400
	 call SaveRestartDataSetup
end subroutine Finalsubroutines
      

!================================================================================
!    Some Auxilary subroutines useful for the shock problem                     
!================================================================================
!Note: currently it does not work, needs to be updated for SOA
subroutine SplitParticlesShock
	integer :: x1,x2 	
	if(ParticleSplitting.eqv..false.) return
	call FindShockXpos_BmagPeak(Xshock)
	if(proc.eq.0) print*,'Current location of the shock is:',Xshock
		
! 	if(g0.gt.1) then
! 		call SplitPrtlInSubDomain_Rel(Xshock,Xshock+delX_PrtlSplit,yborders(0),yborders(nSubDomainsY),zborders(0),zborders(size(zborders)-1))
! 	else
! 		call SplitPrtlInSubDomain_NonRel(Xshock,Xshock+delX_PrtlSplit,yborders(0),yborders(nSubDomainsY),zborders(0),zborders(size(zborders)-1))
! 	end if
end subroutine SplitParticlesShock

end module setup