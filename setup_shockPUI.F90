module setup
     use parameters
     use vars
     use help_setup
	 use memory
	 use savedata
	 use loadbalance
     implicit none 
contains 

subroutine InitUser
	!------------------------------------------------------------
    !Physical Parameters
	!------------------------------------------------------------
	Alfven_speed=0.2_psn ! Alfven Speed (in units of c)
	PlasmaBeta=0.2 ! upstream ion plasma beta
	DriftUpstream=-0.05!Upstream drift speed 
	dxInjector=abs(0.0_psn*DriftUpstream) ! displacment of the injector per time step
    TeTiRatio=1.0 ! Electron to Ion Temperature
	PUIratio=0.5!fraction of pickup ions 
	Btheta=90*(pi/180.0_psn) !Angle (in radians) between magnetic field and z-axis
	Bphi= 90*(pi/180.0_psn)! Angle between x-y projection of magnetic field and x-axis   
	!------------------------------------------------------------
    !Setup Parameters (all positions are in global grid-cell cordiantes)
	!------------------------------------------------------------	
	xcond=4 !global x-cord of the conductor reflecting EM waves
	PrtlWall=xcond+128! global x-cord of the Wall that reflects particles 
    xinj=PrtlWall+10000!+128!intial position of the injector, in global cordinates 
	BoxSizeIncrement=128*2 !The box-size is incremented by this amount when needed
	
	!------------------------------------------------------------
    !Derived quantities
	!------------------------------------------------------------	
    Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
    Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
	call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
    call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
	
	if(abs(DriftUpstream).gt.1.0_psn) then 
		dxUpstream=c*sqrt((DriftUpstream-1.0_psn)*(DriftUpstream+1.0_psn))/abs(DriftUpstream)
	else 
		dxUpstream=c*abs(DriftUpstream)
	end if
    BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)   ! External magentic field  
    Bz_ext0=BextMag*cos(Btheta)
    Bx_ext0=BextMag*sin(Btheta)*cos(Bphi)
	By_ext0=BextMag*sin(Btheta)*sin(Bphi)
    driftEx=0.0
    driftEy=DriftUpstream*Bz_ext0
    driftEz=-DriftUpstream*By_ext0
	
end subroutine InitUser



subroutine InjectNewPrtl
	 implicit none 
	 integer  :: procxind_this
     real(psn):: x1,x2,xinj_local
	 
	 procxind_this=procxind(proc)
	 
	 
	 !x1 and x2 determined boundary of the region in local cordinates where new upstream plasma is injected 
	 x1=max(xinj,real(xborders(procxind_this)))
	 x2=min(xinj+dxInjector,real(xborders(procxind_this+1)))
	 
	 x1=x1-xborders(procxind_this)+3 !change to local cordinate 
	 x2=x2-xborders(procxind_this)+3 
	 
	 xinj_local=max(0.5_psn,xinj-xborders(procxind_this)+3)
	 
	 call ReflectNewPrtlRight(xinj_local)  ! Is it needed if I am using steady injector?

	 call HomogeneousDriftingPrtl(x1,x2,DriftUpstream) 	 	 
	 
	 xinj=xinj+dxInjector ! update location of the injector 
	 if(proc.eq.0) print*,'Injector is at:',xinj
	 		  
end subroutine InjectNewPrtl

subroutine HomogeneousDriftingPrtl(x1,x2,drift)
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
end subroutine HomogeneousDriftingPrtl



end module setup 