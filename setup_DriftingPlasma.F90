!this module set up the pic simulation for alfven turbulence cascade 
!SetupType 1: Standing Alfven waves,i.e. initially dE=0, dB=0
!SetupType 2: Shear Alfven wave 

module setup

use parameters
use vars
use help_setup
!use loadbalance
implicit none 

real(psn) :: Temp_ion,Temp_elc,Alfven_speed,PlasmaBeta
real(psn) :: BextMag,Btheta,TeTiRatio,drift
integer, parameter :: GammaTableLength=10000
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
contains 
     
subroutine InitPosVel
          real(psn) :: ugamma,vgamma,wgamma,drift
          real(psn) :: xlocal,ylocal,zlocal
		  integer   :: i,tag
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------          
          Alfven_speed=0.01 !in units of c 
          PlasmaBeta=0.00!3
		  drift=0.1
          TeTiRatio=1.0     
          BTheta=90*(pi/180) !Angle (radian) between magnetic field and z-axis  
          !--------------------------------------------------
                           
          !set the background magentic field 
          BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
          Bz_ext0=BextMag*cos(Btheta)
          Bx_ext0=BextMag*sin(Btheta)
          !thermal spread in the velocity, 
          Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
          Temp_elc=Temp_ion*TeTiRatio*(mi/me) !Temp_elc=kT_e/m_ec^2
		  
  		  call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
          call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
          
!--------------------!!create ion-electron pairs!!-------------------------! 
		   tag=0
           do i=1,Nelc
               call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)   
			   call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
			   call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,t)
			   call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)
			   call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,t)
          end do                 
!-----------------------------------------------------------------------------!
end subroutine InitPosVel

     
!================================================================================
!The following subroutines are meant for customisation and must be part of the setup module,
!even though they are not used and the default subroutines are used instead
!Note:: in order to enable use of the following subroutines the corresponding custom 
!variables must be change to .true.
!================================================================================
subroutine InitOverride!override the default intial conditons, such as initial fields, intial current, external magnetic field 
end subroutine InitOverride  

subroutine PostMovDep
end subroutine PostMovDep

subroutine Finalsubroutines !this subroutine is deleted, see shock_setup 
     !define here all additional actions that is to be perform at the end of every time step
     !call LoadBalanceHomogeneous
end subroutine Finalsubroutines


end module setup