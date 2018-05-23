!Weibel Instability: counter-streaming relativisitc beams 
module setup

use parameters
use vars
use help_setup
implicit none 

real(psn) :: Temp_ion,Temp_elc,Alfven_speed
real(psn) :: BextMag,Btheta,TeTiRatio
contains 
     
subroutine InitUser
          real(psn) :: ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,drift,gamma_beta
          real(psn) :: xlocal,ylocal,zlocal
		  real(psn) :: r1,r2,r3
		  integer   :: i
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------          
          Alfven_speed=0.0!in units of c 
          TeTiRatio=1.0 ! Temperature of electron/ Temperature of the ions  
          BTheta=90*(pi/180) !Angle (radian) between magnetic field and z-axis  
		  !--------------------------------------------------                         
          !set the background magentic field 
          BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
          Bz_ext0=BextMag*cos(Btheta)
          Bx_ext0=BextMag*sin(Btheta)
		  
          !Initial Temperatures of the ions and the electrons 
          !Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
          !Temp_ion=0.0!in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
		  !Temp_elc=Temp_ion*TeTiRatio*(mi/me) !Temp_elc=kT_e/m_ec^2          


!--------------------!!create ion-electron pairs!!-------------------------!           
vgamma=0.0
wgamma=0.0
gamma_beta=sqrt((1.0_psn+g0)*(g0-1.0_psn))
		  !ions		   
           do i=1,Nelc
			   
               call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)  ! defined in help_setup.F90  		   
			   if(i.le.Nelc/2) ugamma=gamma_beta
			   if(i.gt.Nelc/2) ugamma=-gamma_beta
               !call GetVelGamma_MaxwellNonRel(drift,Temp_ion,ugamma,vgamma,wgamma)
			   call InsertParticleAt(i,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
			   !call InsertParticleAt(array index,x,y,z,u,v,w,charge ,Tag ID,Flavour ID, a dummy multipurpose variable)                                                            
          end do
		  
          !electrons 		  
          do i=1,Nelc
			   if(i.le.Nelc/2) ugamma=-gamma_beta
			   if(i.gt.Nelc/2) ugamma=gamma_beta
               !call GetVelGamma_MaxwellNonRel(drift,Temp_elc,ugamma,vgamma,wgamma)
			   call InsertParticleAt(i+Nelc,xp(i),yp(i),zp(i),ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn) !electrons are initialized at the same location as theions
          end do 

end subroutine InitUser

     
!================================================================================
!The following subroutines are meant for customisation and must be part of the setup module
!================================================================================
!This subroutine is called just one right before entering the main interation loop. The default values of the vriables can be changed at this point.       
subroutine InitOverride

end subroutine InitOverride  

! Some setup sepcific actions can be defiend in the folliwng three subroutine
!This subroutine is called once after the end of every time step. Usage example: add some new particles at a boundary. 
subroutine Finalsubroutines 
     
end subroutine Finalsubroutines


!This subroutine is called right after moving the particles. Usage example: scatter particles
subroutine PostMovDep

end subroutine PostMovDep

!This subroutine is called right before including the current in the Maxwell's equations. Usage example: add some external (antena) current
subroutine PreAddCurrent

end subroutine PreAddCurrent


end module setup