!this module set up the pic simulation for alfven turbulence cascade 
!SetupType 1: Standing Alfven waves,i.e. initially dE=0, dB=0
!SetupType 2: Shear Alfven wave 

module setup

use parameters
use vars
use help_setup
implicit none 

contains 
     
subroutine InitPosVel
          real(psn) :: ugamma,vgamma,wgamma,Temp_ion,Temp_elc
          real(psn) :: xlocal,ylocal,zlocal
		  integer   :: i

                           
        Temp_ion=0.01_psn !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
        Temp_elc=0.01_psn*(mi/me)! Temp_elc=kT_e/m_ec^2
		call InitNewFlvr(FlvID=3,QbyM=qme,Type=0,SaveFld=1,Split=0)
		        
!--------------------!!create ion-electron pairs!!-------------------------! 
         !ions: normal weight 
          do i=1,Nelc
               call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)     
               call GetVelGamma_MaxwellNonRel(0.0_psn,Temp_ion,ugamma,vgamma,wgamma)
			   call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0)
			   
               call GetVelGamma_MaxwellNonRel(0.0_psn,Temp_elc,ugamma,vgamma,wgamma)
			   call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,0,2,0)
			   
			   call GetVelGamma_MaxwellNonRel(0.0_psn,Temp_elc,ugamma,vgamma,wgamma)
			   call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,0,3,0)
          end do
!           !electrons: light weight
!           do i=1,Nelc
! 			   call GetVelGamma_MaxwellJuttner(0.0_psn,Temp_elc,ugamma,vgamma,wgamma)
!                call InsertNewPrtl(p(i)%x,p(i)%y,p(i)%z,ugamma,vgamma,wgamma,-1.0_psn,0,2,0)
!           end do                  
!-----------------------------------------------------------------------------!     
          
end subroutine InitPosVel

     
!================================================================================
!The following subroutines are meant for customisation and must be part of the setup module,
!even though they are not used and the default subroutines are used instead
!Note:: in order to enable use of the following subroutines the corresponding custom 
!variables must be change to .true.
!================================================================================
subroutine InitOverride !override the default intial conditons, such as initial fields, intial current, external magnetic field 
end subroutine InitOverride  

subroutine Finalsubroutines !this subroutine is deleted, see shock_setup 
     !define here all additional actions that is to be perform at the end of every time step
     !call LoadBalanceHomogeneous
end subroutine Finalsubroutines
      
subroutine PostMovDep
end subroutine PostMovDep	  




end module setup