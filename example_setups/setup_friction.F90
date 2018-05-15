module setup
     use parameters
     use vars
	 use help_setup
	 implicit none
	 
integer, parameter :: GammaTableLength=10000
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
real(psn) :: Temp_ion,Temp_elc
real(psn) :: Alfven_speed,PlasmaBeta,Btheta,Bphi,BextMag,TeTiRatio
real(psn) :: xfast,yfast,zfast,xfast_local,yfast_local,zfast_local

 contains 

subroutine InitUser
	real(psn) :: xlocal,ylocal,zlocal,ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,beta      
	integer :: i
    !Physical parameters 
	Alfven_speed=0.0_psn
	PlasmaBeta=0.1
	Btheta=(pi/180)*2.0
	Bphi=0.0

	TeTiRatio=1
	xfast=10
	yfast=ny/2
	zfast=0.5_psn
	

	
    !set the background magentic field 
    BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
    Bz_ext0=BextMag*cos(Btheta)
    Bx_ext0=BextMag*sin(Btheta)*cos(Bphi)
	By_ext0=BextMag*sin(Btheta)*sin(Bphi)
	
    Temp_ion=0.01!PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
    Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
    
	call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
    call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)

	!background ions
	do i=1,Nelc
       call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)
       call GetVelGamma_MaxwellJuttner(0.0,Temp_ion,ugamma,vgamma,wgamma)
	   call InsertParticleAt(i,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
   end do
   
   !background electrons 		  
   do i=1,Nelc
	    call GetVelGamma_MaxwellJuttner(0.0,Temp_elc,ugamma,vgamma,wgamma)  
        call InsertParticleAt(i+Nelc,xp(i),yp(i),zp(i),ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn) 
   end do   
      
   !Insert a Fast Moving Ion
   call InitNewFlvr(FlvID=3,QbyM=qmi,Type=0,SaveFld=1,Split=0)
   xfast_local=xfast-xborders(procxind(proc))+3
   yfast_local=yfast-yborders(procyind(proc))+3
   zfast_local=zfast-zborders(proczind(proc))+3
   
#ifndef twoD   
   if((xfast_local.ge.3.and.xfast_local.lt.mx-2).and.(yfast_local.ge.3.and.yfast_local.lt.my-2).and.(zfast_local.ge.3.and.zfast_local.lt.mz-2)) then 
#else
   if((xfast_local.ge.3.and.xfast_local.lt.mx-2).and.(yfast_local.ge.3.and.yfast_local.lt.my-2)) then 
#endif	   
       beta=0.8_psn
	   ugamma=beta/sqrt((1.0_psn-beta)*(1.0_psn+beta))
	   vgamma=0.0
	   wgamma=0.0
	   call InsertParticleAt(2*Nelc+1,xfast_local,yfast_local,zfast_local,ugamma,vgamma,wgamma,1.0_psn,1,3,0.0_psn) 	
	   call GetVelGamma_MaxwellJuttner(0.0,Temp_elc,ugamma,vgamma,wgamma) 
	   call InsertParticleAt(2*Nelc+2,xfast_local,yfast_local,zfast_local,ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn) !balance the charge	 	    
   end if    	  
end subroutine InitUser


	 
subroutine InitOverride !override the default intial conditons
end subroutine InitOverride  

!The following subroutine is called after moving particles and depositing the current on grid 
! Any additional step which are setup specific should be implelemnted here 
subroutine PostMovDep

end subroutine PostMovDep

subroutine PreAddCurrent
end subroutine PreAddCurrent

subroutine Finalsubroutines !this subroutine is called at the end of every time step	

end subroutine Finalsubroutines	 


end module setup 