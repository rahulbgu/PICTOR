!this module set up the pic simulation for alfven turbulence cascade 
!SetupType 1: Standing Alfven waves,i.e. initially dE=0, dB=0
!SetupType 2: Shear Alfven wave 

module setup

use parameters
use vars
use help_setup
!use loadbalance
use DerivedQnty
use interpolation
use savedata_routines, only : CalcPrtlChargeFlux
use fields, only : UpdateCurrentsAllEdges
implicit none 

real(psn) :: Temp_ion,Temp_elc
real(psn) :: Alfven_speed 
real(psn) :: PlasmaBeta
real(psn) :: BextMag,dbeta 
!Variables to define Initial spectrum 
integer            :: SetupType,DriveTurb
logical            :: SplitTestParticles
real(psn)          :: ScatteringRate
integer, parameter :: Nmodes=7 ! number of modes at the begining 
real(psn), dimension(Nmodes,3):: vecK, vecKXB, vecKXBXB, vecKXdB
real(psn), dimension(Nmodes)  :: AmpK, PhaseK, InitialAmpK 
real(psn) :: B2_Now,InitialB2 ! used in driving turublence   
!integer, parameter :: MFP_PDF_TableSize=100000
!real(psn), dimension(MFP_PDF_TableSize) :: MFP_PDF_Table

contains 
     
subroutine InitUser
          real(psn) :: ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,beta_drift,TeTiRatio,Btheta
          real(psn) :: xlocal,ylocal,zlocal
          integer   :: i,NTestPrtl,TPtag_ratio
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------
          SetupType=2 !see the header for detail
          DriveTurb=0
          SplitTestParticles=.false.
          
          Alfven_speed=0.2_psn!in units of c 
          PlasmaBeta=0.1_psn!3
          TeTiRatio=1.0_psn     
          BTheta=90*(pi/180) !Angle (radian) between magnetic field and z-axis  
		  ScatteringRate=0.01 ! randian/time-step
#ifdef twoD          
          vecK(1,:)=(/ 2*pi/nx, 2*pi/ny, 0.0_psn /)
          vecK(2,:)=(/ 2*pi/nx, 4*pi/ny, 0.0_psn /)
          vecK(3,:)=(/-2*pi/nx,-2*pi/ny, 0.0_psn /)
          vecK(4,:)=(/-2*pi/nx,-4*pi/ny, 0.0_psn /)
          vecK(5,:)=(/ 2*pi/nx, 0.0_psn    , 0.0_psn /)
          vecK(6,:)=(/ 2*pi/nx, 0.0_psn    , 0.0_psn /)
          vecK(7,:)=(/ 2*pi/nx, 0.0_psn    , 0.0_psn /)
#else           
          vecK(1,:)=(/  2*pi/nx, 2*pi/ny, 0.0_psn /)
          vecK(2,:)=(/  2*pi/nx, 4*pi/ny, 0.0_psn /)
          vecK(3,:)=(/ -4*pi/nx, 2*pi/ny, 0.0_psn /)
          vecK(4,:)=(/ -2*pi/nx, 0.0_psn    ,-2*pi/nz /)
          vecK(5,:)=(/ -2*pi/nx, 0.0_psn    ,-4*pi/nz /)
          vecK(6,:)=(/  4*pi/nx, 0.0_psn    ,-2*pi/nz /)
          vecK(7,:)=(/  2*pi/nx, 0.0_psn    , 0.0_psn /)
#endif          
      
          AmpK=(/ 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.0/)
		  AmpK=AmpK*0.0
          !AmpK=(/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
          !PhaseK=(/ 0.1*pi, 0.4*pi, 0.7*pi, 1.1*pi, 1.4*pi, 1.7*pi, 0.7*pi /)
		  PhaseK=(/ 0.1*pi, 0.8*pi, 1.2*pi, 1.5*pi, 0*pi, 0.4*pi, 1.7*pi /)
          !call RandomiseAlfvenWavePhase

          NTestPrtl=int(Nelc/(128*128*64)) ! number of Test Particles are TP_ratio times less than number of simulation particles 
          TPtag_ratio=1
          !--------------------------------------------------
                           

          !set the background magentic field 
          BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
          Bz_ext0=BextMag*cos(Btheta)
          Bx_ext0=BextMag*sin(Btheta)
          !thermal spread in the velocity, 
          Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
          Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
          call InitKXB !initialise unit vectors along kXB
          call InitKXBXB ! initialise unit vectors along VXB_ext, useful for determining direction of E 
          call InitKXdB ! unit vector along kXdB, useful in setting curlB=J
          
!--------------------!!create ion-electron pairs!!-------------------------! 
         !ions 
           do i=1,Nelc
               call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)      
               call GetDriftVelocity(xlocal,ylocal,zlocal,ubeta,vbeta,wbeta,qmi)
               beta_drift=sqrt(ubeta**2+vbeta**2+wbeta**2)
               call GetVelGamma_MaxwellJuttner(beta_drift,Temp_ion,ugamma,vgamma,wgamma,direction=(/ubeta,vbeta,wbeta/))
               call InsertParticleAt(i,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
          end do
          !electrons 
          do i=1,Nelc
               call GetDriftVelocity(xp(i),yp(i),zp(i),ubeta,vbeta,wbeta,qme)
               beta_drift=sqrt(ubeta**2+vbeta**2+wbeta**2)
               call GetVelGamma_MaxwellJuttner(beta_drift,Temp_elc,ugamma,vgamma,wgamma,direction=(/ubeta,vbeta,wbeta/))          
               call InsertParticleAt(i+Nelc,xp(i),yp(i),zp(i),ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn)                     
          end do 
                    
          call SetQbyM(1,qmi)
          call SetQbyM(2,qme)
          call SetFlvrSaveFldData(1,1)
          call SetFlvrSaveFldData(2,1)
          !Need to call a subroutine to tag particles 
		  !call InitPoissonDist(MFP_PDF_TableSize,MFP_PDF_Table,real(ScatterPeriod))
!-----------------------------------------------------------------------------!     


!=================Test Particles =============================================!

!--- used on TITAN ------

call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=3,QbyM=qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=4,QbyM=qme,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=5,QbyM=-qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=6,QbyM=-qme,SaveFldData=0)

call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=7,QbyM=qmi/4,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=8,QbyM=qmi/2,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=9,QbyM=2*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=10,QbyM=4*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=11,QbyM=8*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=12,QbyM=16*qmi,SaveFldData=0)

call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=13,QbyM=-qmi/4,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=14,QbyM=-qmi/2,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=15,QbyM=-2*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=16,QbyM=-4*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=17,QbyM=-8*qmi,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=18,QbyM=-16*qmi,SaveFldData=0)

call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=19,QbyM=qme/2,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=20,QbyM=qme/4,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=21,QbyM=qme/8,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=22,QbyM=qme/16,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=23,QbyM=qme/32,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=24,QbyM=qme/64,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=25,QbyM=qme/128,SaveFldData=0)

call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=26,QbyM=-qme/2,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=27,QbyM=-qme/4,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=28,QbyM=-qme/8,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=29,QbyM=-qme/16,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=30,QbyM=-qme/32,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=31,QbyM=-qme/64,SaveFldData=0)
call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=32,QbyM=-qme/128,SaveFldData=0)


!Another (OLDER) set of test particles 
!
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=3,QbyM=qmi/1.5,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=4,QbyM=qmi/2,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=5,QbyM=qmi/3,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=6,QbyM=qmi/4,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=7,QbyM=qmi/5,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=8,QbyM=qmi/6,SaveFldData=0)
!
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=8*Temp_ion,TagRatio=TPtag_ratio,FlvID=9,QbyM=qmi/1.5,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=8*Temp_ion,TagRatio=TPtag_ratio,FlvID=10,QbyM=qmi/2,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=8*Temp_ion,TagRatio=TPtag_ratio,FlvID=11,QbyM=qmi/4,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=8*Temp_ion,TagRatio=TPtag_ratio,FlvID=12,QbyM=qmi/6,SaveFldData=0)
!
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=0.125*Temp_ion,TagRatio=TPtag_ratio,FlvID=13,QbyM=qmi/1.5,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=0.125*Temp_ion,TagRatio=TPtag_ratio,FlvID=14,QbyM=qmi/2,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=0.125*Temp_ion,TagRatio=TPtag_ratio,FlvID=15,QbyM=qmi/4,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=0.125*Temp_ion,TagRatio=TPtag_ratio,FlvID=16,QbyM=qmi/6,SaveFldData=0)
!
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=17,QbyM=qmi,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=18,QbyM=qme,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=19,QbyM=qmi/2,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=20,QbyM=qmi/4,SaveFldData=0)
!
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_ion,TagRatio=TPtag_ratio,FlvID=21,QbyM=-qmi,SaveFldData=0)
! call InsertTestParticlesDriftMaxwellian(Nprtl=NTestPrtl,Temp=Temp_elc,TagRatio=TPtag_ratio,FlvID=22,QbyM=-qme,SaveFldData=0)



!=============================================================================!

          
end subroutine InitUser


subroutine InsertTestParticlesDriftMaxwellian(Nprtl,Temp,TagRatio,FlvID,QbyM,SaveFldData)
     integer,   intent(IN):: Nprtl,TagRatio,FlvID,SaveFldData
     real(psn), intent(IN):: Temp,QbyM   
	 real(psn)            :: qm_this
     real(xpsn)           :: xlocal
     real(ypsn)           :: ylocal
     real(zpsn)           :: zlocal
     real(psn)            :: ugamma,vgamma,wgamma
     real(psn)            :: ubeta,vbeta,wbeta,beta_drift
     integer              :: n
     integer              :: ptagged,ptagID
     ptagged= proc*(Nprtl/TagRatio)
     do n=1,Nprtl
 	        if(mod(n,TagRatio).eq.0) then 
	            ptagID=ptagged+1
	            ptagged=ptagged+1 
	        else
	            ptagID=0 
	        end if 
               call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)
			   if(QbyM.gt.0) then 
			       qm_this=max(QbyM,qmi) ! to stop heavy ions from getting very large polarization drift 
			   else 
				   qm_this=QbyM
				   if(abs(qm_this).gt.qmi) qm_this=-qmi ! to stop heavy ions from getting large drift speed
			   end if 
               call GetDriftVelocity(xlocal,ylocal,zlocal,ubeta,vbeta,wbeta,qm_this) ! qmi is used for even heavy ions to keep polarisation drift small 
               beta_drift=sqrt(ubeta**2+vbeta**2+wbeta**2)
               if(beta_drift.gt.1) then 
                    print*,'beta drfit exceeded speed of light beta_drfit:',beta_drift
                    STOP
               end if
               call GetVelGamma_MaxwellJuttner(beta_drift,Temp,ugamma,vgamma,wgamma,direction=(/ubeta,vbeta,wbeta/))
               call InsertNewTestPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,ptagID,FlvID,0.0_psn)
     end do 
     call SetQbyM(FlvID,QbyM)
     call SetFlvrSaveFldData(FlvID,SaveFldData)
     call SetFlvrType(FlvID,-1)
     FlvrSplitPrtl(FlvID)=0
     CurrentTagID(FlvID)=ptagged+1
end subroutine InsertTestParticlesDriftMaxwellian


!The following subroutine creates only shear Alfven waves 
subroutine GetDriftVelocity(xlocal,ylocal,zlocal,ubeta,vbeta,wbeta,qm)
     real(psn) :: ubeta,vbeta,wbeta
     real(psn) :: xlocal,ylocal,zlocal
     real(dbpsn)::xglobal,yglobal,zglobal
     real(psn) :: qm
     real(psn) :: upol_drift,vpol_drift,wpol_drift
      
    xglobal=xlocal-3+xborders(procxind(proc))
    yglobal=ylocal-3+yborders(procyind(proc))
#ifdef twoD     
    zglobal=zlocal
#else
    zglobal=zlocal-3+zborders(proczind(proc))
#endif

     
     !EXB drift 
     call GetEXBDrift(xglobal,yglobal,zglobal,ubeta,vbeta,wbeta)
     !Polarisation drift, due to time dependent electric field 
	 ! Can be large if the simulation box size is small
      if(SetupType.eq.2) then
            call GetCurrent(xglobal,yglobal,zglobal,upol_drift,vpol_drift,wpol_drift,qm)
            ubeta=ubeta+upol_drift
            vbeta=vbeta+vpol_drift
            wbeta=wbeta+wpol_drift
      end if          
end subroutine GetDriftVelocity

!set the initial value of electric and magnetic fields corresponding to superposition of all intial Alfven modes 
subroutine InitEMfieldAlfvenWave 
     integer :: i,j,k
     real(dbpsn):: xglobal,yglobal,zglobal
     real(psn)  :: Ax,Ay,Az
     do k=1,mz
          do j=1,my
               do i=1,mx
                    !Set the value of purturbed magnetic field 
                    call GetGlobalPosition_DP(real(i),real(j+0.5),real(k+0.5),xglobal,yglobal,zglobal)
                    call GetdelB(xglobal,yglobal,zglobal,Ax,Ay,Az)
                    Bx(i,j,k)=Bx(i,j,k)+Ax
                    
                     call GetGlobalPosition_DP(real(i+0.5),real(j),real(k+0.5),xglobal,yglobal,zglobal)
                     call GetdelB(xglobal,yglobal,zglobal,Ax,Ay,Az)
                     By(i,j,k)=By(i,j,k)+Ay
                    
                     call GetGlobalPosition_DP(real(i+0.5),real(j+0.5),real(k),xglobal,yglobal,zglobal)
                     call GetdelB(xglobal,yglobal,zglobal,Ax,Ay,Az)
                     Bz(i,j,k)=Bz(i,j,k)+Az
                    
                    !set the perturbed electric field 
                    call GetGlobalPosition_DP(real(i+0.5),real(j),real(k),xglobal,yglobal,zglobal)
                    call GetdelE(xglobal,yglobal,zglobal,Ax,Ay,Az)
                    Ex(i,j,k)=Ex(i,j,k)+Ax
                    
                    call GetGlobalPosition_DP(real(i),real(j+0.5),real(k),xglobal,yglobal,zglobal)
                    call GetdelE(xglobal,yglobal,zglobal,Ax,Ay,Az)
                    Ey(i,j,k)=Ey(i,j,k)+Ay
                    
                    call GetGlobalPosition_DP(real(i),real(j),real(k+0.5),xglobal,yglobal,zglobal)
                    call GetdelE(xglobal,yglobal,zglobal,Ax,Ay,Az)
                    Ez(i,j,k)=Ez(i,j,k)+Az
                    
               end do
          end do
     end do
end subroutine InitEMfieldAlfvenWave

!the following subroutine gives right drift to the electrons and ions such that the curlB=J
! subroutine GetPolarisationDrift
! end subroutine GetPolarisatinDrift
subroutine GetCurrent(xglobal,yglobal,zglobal,ubeta,vbeta,wbeta,qm)
     real(dbpsn)::xglobal,yglobal,zglobal
     real(psn) :: ubeta,vbeta,wbeta
     real(psn) :: qm
     real(psn) :: Ax,Ay,Az
     
     call SumAmpAllModesCurlB(xglobal,yglobal,zglobal,Ax,Ay,Az)
     if(qm.gt.0) then 
        ubeta=Ax*(qmi/qm)/((qi*epc)*(1+me/mi)) !(Important: qmi/qm is changed to 1, becasue heavy particles get relatively large drift)
        vbeta=Ay*(qmi/qm)/((qi*epc)*(1+me/mi))
        wbeta=Az*(qmi/qm)/((qi*epc)*(1+me/mi))
     else
        ubeta=Ax*abs(qme/qm)*(me/mi)/((qe*epc)*(1+me/mi))
        vbeta=Ay*abs(qme/qm)*(me/mi)/((qe*epc)*(1+me/mi))
        wbeta=Az*abs(qme/qm)*(me/mi)/((qe*epc)*(1+me/mi))
     end if
end subroutine GetCurrent

subroutine GetEXBDrift(x,y,z,Ampx,Ampy,Ampz)
     real(dbpsn), intent(IN) :: x,y,z
     real(dbpsn)             :: CosPhase
     real(psn), intent(INOUT):: Ampx,Ampy,Ampz
     integer                 :: i
     real(psn)               :: KdotB0,signKdotB0
     Ampx=0.0_psn
     Ampy=0.0_psn
     Ampz=0.0_psn
     do i=1,Nmodes
        KdotB0=vecK(i,1)*Bx_ext0+vecK(i,2)*By_ext0+vecK(i,3)*Bz_ext0
          signKdotB0=sign(1.0_psn,KdotB0)
          CosPhase=cos(vecK(i,1)*x+vecK(i,2)*y+vecK(i,3)*z+PhaseK(i))
          Ampx=Ampx-AmpK(i)*signKdotB0*Alfven_speed*vecKXB(i,1)*real(CosPhase)
          Ampy=Ampy-AmpK(i)*signKdotB0*Alfven_speed*vecKXB(i,2)*real(CosPhase)
          Ampz=Ampz-AmpK(i)*signKdotB0*Alfven_speed*vecKXB(i,3)*real(CosPhase)
     end do  
end subroutine GetEXBdrift

subroutine GetdelB(x,y,z,Ampx,Ampy,Ampz)
     real(dbpsn), intent(IN) :: x,y,z
     real(dbpsn)             :: CosPhase
     real(psn), intent(INOUT):: Ampx,Ampy,Ampz
     integer                 :: i
     Ampx=0.0
     Ampy=0.0
     Ampz=0.0
     do i=1,Nmodes
          CosPhase=cos(vecK(i,1)*x+vecK(i,2)*y+vecK(i,3)*z+PhaseK(i))
          Ampx=Ampx+AmpK(i)*BextMag*vecKXB(i,1)*real(CosPhase)
          Ampy=Ampy+AmpK(i)*BextMag*vecKXB(i,2)*real(CosPhase)
          Ampz=Ampz+AmpK(i)*BextMag*vecKXB(i,3)*real(CosPhase)     
     end do 
          
end subroutine GetdelB

subroutine GetdelE(x,y,z,Ampx,Ampy,Ampz)
     real(dbpsn), intent(IN) :: x,y,z
     real(dbpsn)             :: CosPhase
     real(psn), intent(INOUT):: Ampx,Ampy,Ampz
     integer                 :: i
     real(psn)               :: KdotB0,signKdotB0 
     Ampx=0.0
     Ampy=0.0
     Ampz=0.0
     do i=1,Nmodes     
          KdotB0=vecK(i,1)*Bx_ext0+vecK(i,2)*By_ext0+vecK(i,3)*Bz_ext0
          signKdotB0=sign(1.0_psn,KdotB0)
          CosPhase=cos(vecK(i,1)*x+vecK(i,2)*y+vecK(i,3)*z+PhaseK(i))
          Ampx=Ampx+AmpK(i)*signKdotB0*Alfven_speed*vecKXBXB(i,1)*BextMag*real(CosPhase)
          Ampy=Ampy+AmpK(i)*signKdotB0*Alfven_speed*vecKXBXB(i,2)*BextMag*real(CosPhase)
          Ampz=Ampz+AmpK(i)*signKdotB0*Alfven_speed*vecKXBXB(i,3)*BextMag*real(CosPhase)
     end do  
     !print*,'Amp',Ampx,Ampy,Ampz
end subroutine GetdelE

subroutine SumAmpAllModesCurlB(x,y,z,Ampx,Ampy,Ampz)
     real(dbpsn), intent(IN) :: x,y,z
     real(dbpsn)             :: SinPhase
     real(psn), intent(INOUT):: Ampx,Ampy,Ampz
     real(psn)               :: Kmag
     integer                 :: i
     Ampx=0.0
     Ampy=0.0
     Ampz=0.0
     do i=1,Nmodes
          SinPhase=-sin(vecK(i,1)*x+vecK(i,2)*y+vecK(i,3)*z+PhaseK(i))
          Kmag=sqrt(vecK(i,1)**2+vecK(i,2)**2+vecK(i,3)**2)
          Ampx=Ampx+AmpK(i)*Kmag*vecKXdB(i,1)*real(SinPhase)*BextMag
          Ampy=Ampy+AmpK(i)*Kmag*vecKXdB(i,2)*real(SinPhase)*BextMag
          Ampz=Ampz+AmpK(i)*Kmag*vecKXdB(i,3)*real(SinPhase)*BextMag
     end do 
end subroutine SumAmpAllModesCurlB




!------------------------------------------------------------------------------------------------
! Following subroutines are relavant for driving turbulence 
!------------------------------------------------------------------------------------------------
subroutine AddNewAlfvenWaves
     !call RandomiseAlfvenWavePhase
     call CalcSumVecFldSQ(Bx,By,Bz,B2_Now)
     if(t.eq.tstart) then 
          InitialB2=B2_Now
          InitialAmpK=AmpK
     end if
     
     !if(modulo(t,60).ne.0) return 
     
     if(B2_Now.lt.InitialB2) then
          AmpK=sqrt((InitialB2-B2_Now)/(2*InitialB2))*InitialAmpK 
          call InitEMfieldAlfvenWave 
          !call CalcPrtlVelocityField(4)
          call AddDriftToPrtl
     end if
     
end subroutine AddNewAlfvenWaves
subroutine AddDriftToPrtl
     integer ::n
     real(psn) :: driftx,drifty,driftz,DriftGamma,DriftBeta,gamma,projn
     do n=1,prtl_arr_size
          if(qp(n).eq.0) cycle
          !call InterpVectorField_CellCenters(xp(n),yp(n),zp(n),driftx0,drifty0,driftz0,Jx,Jy,Jz)
          call GetDriftVelocity(xp(n),yp(n),zp(n),driftx,drifty,driftz,flvrqm(flvp(n)))
          DriftBeta=sqrt(driftx**2+drifty**2+driftz**2)
          if(DriftBeta.ne.0) then 
               DriftGamma=1.0_psn/sqrt((1.0_psn-DriftBeta)*(1.0_psn+DriftBeta))
               gamma=sqrt(1.0+up(n)**2+vp(n)**2+wp(n)**2)     
               driftx=driftx/DriftBeta
               drifty=drifty/DriftBeta
               driftz=driftz/DriftBeta
               projn=up(n)*driftx+vp(n)*drifty+wp(n)*driftz
              up(n)=up(n)+(DriftGamma-1)*projn*driftx + DriftGamma*DriftBeta*gamma*driftx
              vp(n)=vp(n)+(DriftGamma-1)*projn*drifty + DriftGamma*DriftBeta*gamma*drifty
              wp(n)=wp(n)+(DriftGamma-1)*projn*driftz + DriftGamma*DriftBeta*gamma*driftz
         end if
     end do
end subroutine AddDriftToPrtl

subroutine RandomiseAlfvenWavePhase
     integer :: i
     real(psn) :: r1
     do i=1,Nmodes
          call random_number(r1)
          PhaseK(i)=2*pi*r1
     end do 
end subroutine RandomiseAlfvenWavePhase

!================================================================================================

subroutine InitKXB
     integer :: i
     do i=1,Nmodes
          !Get the unit vectors in the direction of kXB
          call GetUnitVectorAlongAXB(vecK(i,1),vecK(i,2),vecK(i,3),Bx_ext0,By_ext0,Bz_ext0,vecKXB(i,1),vecKXB(i,2),vecKXB(i,3))
     end do 
end subroutine InitKXB

subroutine InitKXBXB
     integer :: i
     do i=1,Nmodes
          call GetUnitVectorAlongAXB(vecKXB(i,1),vecKXB(i,2),vecKXB(i,3),Bx_ext0,By_ext0,Bz_ext0,vecKXBXB(i,1),vecKXBXB(i,2),vecKXBXB(i,3))
     end do 
end subroutine InitKXBXB

subroutine InitKXdB
     integer :: i
     do i=1,Nmodes
          call GetUnitVectorAlongAXB(vecK(i,1),vecK(i,2),vecK(i,3),vecKXB(i,1),vecKXB(i,2),vecKXB(i,3),vecKXdB(i,1),vecKXdB(i,2),vecKXdB(i,3))
    end do 
end subroutine InitKXdB

     
!================================================================================
!The following subroutines are meant for customisation and must be part of the setup module,
!even though they are not used and the default subroutines are used instead
!Note:: in order to enable use of the following subroutines the corresponding custom 
!variables must be change to .true.
!================================================================================
subroutine InitOverride !override the default intial conditons, such as initial fields, intial current, external magnetic field 
     if(SetupType.eq.2) call InitEMfieldAlfvenWave
end subroutine InitOverride  

subroutine Finalsubroutines !this subroutine is deleted, see shock_setup 
     !define here all additional actions that is to be perform at the end of every time step
     !call LoadBalanceHomogeneous
	 !if(ScatteringRate.gt.0) call CoulombScattering ! Scattering is turned off
     
	 !if(DriveTurb.eq.1) call AddNewAlfvenWaves

!      if(SplitTestParticles) then
!           if((spec_save_period.gt.0).and.(modulo(t,spec_save_period).eq.0)) call SplitHighEnergyTailSpeedSpec
!      end if
     
     if(t.gt.15000) fld_save_period=50
     if(t.gt.16000) fld_save_period=400
          

end subroutine Finalsubroutines
subroutine PostMovDep
end subroutine PostMovDep
subroutine PreAddCurrent
end subroutine PreAddCurrent

subroutine CoulombScattering
	integer :: n	
	do n=1,used_prtl_arr_size
		if(qp(n).eq.0) cycle 
		call SmallAngleScattering(up(n),vp(n),wp(n),ScatteringRate)
	end do 	
end subroutine CoulombScattering



! subroutine CoulombScattering
! 	integer :: n
! 	real(psn) :: pVx,pVy,pVz
! 	integer :: index
! 	real(psn) ::r1
! 	!Local turubulent velocity
! 	!call CalcPrtlChargeFlux(2) !only electrons are included
! 	!call UpdateCurrentsAllEdges
!
! 	do n=1,used_prtl_arr_size
! 		if(qp(n).eq.0) cycle
! 		!call InterpVecGridPoints(xp(n),yp(n),zp(n),pVx,pVy,pVz,Jx,Jy,Jz)
! 		if(var1p(n).le.0) then
! 			call HardScattering(up(n),vp(n),wp(n))
! 			call random_number(r1)
! 			call BinarySearch(MFP_PDF_TableSize,MFP_PDF_Table,r1,index)
! 			var1p(n)=MFP_PDF_Table(index)
! 		else
! 			var1p(n)=var1p(n)-1
! 		end if
! 	end do
! end subroutine CoulombScattering

end module setup