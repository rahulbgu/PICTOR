!this module set up the pic simulation for alfven turbulence cascade 
!SetupType 1: Add backgorund current as external current source 
!SetupType 2: Cosmic rays are like particles (Not implemented yet)

module setup

use parameters
use vars
use help_setup
!use loadbalance
implicit none 

real(psn) :: Temp_ion,Temp_elc
real(psn) :: Alfven_speed 
real(psn) :: PlasmaBeta
real(psn) :: BextMag
real(psn) :: JxCR,JyCR,JzCR,fCR,Zeta,DriftCR
!Variables to define Initial spectrum 
integer   :: SetupType


contains
     
subroutine InitUser
          real(psn) :: ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,beta_drift,TeTiRatio,Btheta
          real(psn) :: xlocal,ylocal,zlocal
          integer   :: i,ntp,TPtag_ratio
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------
          SetupType=1!see the header for detail
          
          fCR=0.1_psn!fraction of cosmic rays
          Alfven_speed=0.01_psn!Alfven speed
          !Zeta=10 !4*pi*Jcr/kB
          PlasmaBeta=0.1_psn!3
          TeTiRatio=1.0_psn
          BTheta=90*(pi/180)!Angle (radian) between magnetic field and z-axis                      
          ntp=int(Nelc/256) ! number of Test Particles are TP_ratio times less than number of simulation particles 
          TPtag_ratio=4
          !--------------------------------------------------
          BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)          
          !set the background magentic field 
          Bz_ext0=BextMag*cos(Btheta)
          Bx_ext0=BextMag*sin(Btheta)
          !thermal spread in the velocity, 
          Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5_psn !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
          Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
          JxCR=fCR*qi*epc*c
          JyCR=0.0_psn!fCR*qi*epc*c!0
          JzCR=0.0_psn
!--------------------!!create ion-electron pairs!!-------------------------! 
          ubeta=-fCR*(me/(mi+me))
		  do i=1,Nelc!ions 
                call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)     
                call GetVelGamma_MaxwellNonRel(ubeta,Temp_ion,ugamma,vgamma,wgamma)
                call InsertParticleAt(i,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
          end do
		  
          ubeta=fCR*(mi/(mi+me))          
          do i=1,Nelc!electrons 
                call GetVelGamma_MaxwellNonRel(ubeta,Temp_elc,ugamma,vgamma,wgamma)          
                call InsertParticleAt(i+Nelc,xp(i),yp(i),zp(i),ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn)                     
          end do 
!=============================================================================!     
          
end subroutine InitUser


     
!================================================================================
!The following subroutines are meant for customisation and must be part of the setup module,
!even though they are not used and the default subroutines are used instead
!Note:: in order to enable use of the following subroutines the corresponding custom 
!variables must be change to .true.
!================================================================================
subroutine InitOverride !override the default intial conditons, such as initial fields, intial current, external magnetic field

end subroutine InitOverride  

subroutine Finalsubroutines !this subroutine is called at the end of every time step
     !define here all additional actions that is to be perform at the end of every time step

     !call LoadBalanceHomogeneous
end subroutine Finalsubroutines
      
subroutine PostMovDep
    if(SetupType.eq.1) then 
         call AddBackgroundCurrentCR
    end if 
end subroutine PostMovDep

subroutine PreAddCurrent
end subroutine PreAddCurrent


subroutine AddBackgroundCurrentCR
#ifndef twoD
       Jx(3:mx-3,3:my-3,3:mz-3)=Jx(3:mx-3,3:my-3,3:mz-3)+JxCR
       Jy(3:mx-3,3:my-3,3:mz-3)=Jy(3:mx-3,3:my-3,3:mz-3)+JyCR
       Jz(3:mx-3,3:my-3,3:mz-3)=Jz(3:mx-3,3:my-3,3:mz-3)+JzCR
#else
       Jx(3:mx-3,3:my-3,1)=Jx(3:mx-3,3:my-3,1)+JxCR
       Jy(3:mx-3,3:my-3,1)=Jy(3:mx-3,3:my-3,1)+JyCR
       Jz(3:mx-3,3:my-3,1)=Jz(3:mx-3,3:my-3,1)+JzCR
#endif
end subroutine AddBackgroundCurrentCR



end module setup
