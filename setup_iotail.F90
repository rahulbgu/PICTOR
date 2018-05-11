module setup
     use parameters
     use vars
	 use help_setup
real(psn) :: MeanCount1
integer, parameter :: InjCountTableSize=1000
integer, parameter :: GammaTableLength=10000
real(psn), dimension(InjCountTableSize) :: InjCountTable1
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
real(psn) :: Temp_ion,Temp_elc
real(psn) :: Alfven_speed,PlasmaBeta,Btheta,Bphi,BextMag

 contains 

subroutine InitUser
	    
    !Physical parameters 
	Alfven_speed=0.0003125_psn
	PlasmaBeta=0.2
	
    !set the background magentic field 
    BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
    Bz_ext0=BextMag*cos(Btheta)
    Bx_ext0=BextMag*sin(Btheta)*cos(Bphi)
	By_ext0=BextMag*sin(Btheta)*sin(Bphi)
	
    Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
    Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
    
	call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
    call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
	
#ifdef twoD 
     MeanCount1=0.5_psn*epc*(my-5)!First compute how many new particles are to be added 
#else 	 
     MeanCount1=0.5_psn*epc*(my-5)*(mz-5)
#endif
	call InitPoissonDist(InjCountTableSize,InjCountTable1,MeanCount1)
      
end subroutine InitUser

subroutine InjectNewPrtl
	 implicit none 
	 integer  :: procxind_this
     real(psn):: x1,x2,xinj_local
	 integer  :: xinj
	 
	 xinj=10
	 procxind_this=procxind(proc)
	 
	 !----First Inject Left Moving Particles ----
	 
	 !x1 and x2 determined boundary of the region in local cordinates where new upstream plasma is injected 
	 x1=max(xinj-0.5,real(xborders(procxind_this)))
	 x2=min(xinj+dxInjector,real(xborders(procxind_this+1)))
	 
	 x1=x1-xborders(procxind_this)+3 !change to local cordinate 
	 x2=x2-xborders(procxind_this)+3 
	 
	 call ShearedDriftPrtl(x1,x2) 	 		  
end subroutine InjectNewPrtl

subroutine ShearedDriftPrtl(x1,x2)
	real(psn) :: dx
	real(psn) :: x1,x2
	real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
	real(psn) :: drift
	integer   :: Nelc_New,i,tag
	
	
	tag=0
	drift=0.1
	
	call GetInjPrtlCount(MeanCount1,Nelc_New,InjCountTableSize,InjCountTable1)
    if(Nelc_New.gt.0) then
		call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,ymin,ymax,zmin,zmax)
		call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,0)
		
		call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,0)
	end if 
end subroutine ShearedDriftPrtl

subroutine 







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



	 
subroutine InitOverride !override the default intial conditons
	!xborders(0)=xwall
end subroutine InitOverride  

!The following subroutine is called after moving particles and depositing the current on grid 
! Any additional step which are setup specific should be implelemnted here 
subroutine PostMovDep

end subroutine PostMovDep

subroutine PreAddCurrent
end subroutine PreAddCurrent

subroutine Finalsubroutines !this subroutine is called at the end of every time step	 
	 call InjectNewPrtl
end subroutine Finalsubroutines	 


end module setup 