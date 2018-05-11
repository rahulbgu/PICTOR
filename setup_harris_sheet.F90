!Set up for magnetic reconnection 
module setup

use parameters
use vars
use help_setup
use memory
use loadbalance
implicit none 

integer, parameter :: InjCountTableSize=1000
integer, parameter :: GammaTableLength=10000
real(psn) :: Temp_ion,Temp_elc,Alfven_speed,PlasmaBeta
real(psn) :: BextMag,Btheta,TeTiRatio
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
integer :: TopEdge_Ypos,BottomEdge_Ypos,LeftEdge_Xpos,RightEdge_Xpos
real(psn) :: SheetThickness,BackgroundPrtlFraction,AmpB_Harris
contains 
     
subroutine InitUser
          real(psn) :: ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,drift
          real(psn) :: xlocal,ylocal,zlocal
		  real(psn) :: r1,r2,r3
		  integer   :: i
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------          
		  Alfven_speed=0.01_psn !in units of c 
          TeTiRatio=1.0_psn  ! Electron Temperature / Ion Temperature   
		  SheetThickness=10*c_ompe! Thickness of the current sheet; c_compe is the pnumber of cells per electron skin-depth
          BackgroundPrtlFraction=0.1 ! the background plasma density / peak density of current bearing particles in the current sheet
                           
          !set the value of magentic field amplitudes 
          AmpB_Harris=Alfven_speed*sqrt(epc*(massi+masse)*c*c)   		  
          !Bz_ext0=AmpB_Harris ! Guide Field 
          
		  !Position of the boundaries of the box in the global grid
		  TopEdge_Ypos=ny-12
		  BottomEdge_Ypos=12          
		  LeftEdge_Xpos=0
		  RightEdge_Xpos=nx!-12    
          !-------------------------------------------------
          ! Initialize speed distribution functions for the ions and electrons
          !-------------------------------------------------     	        
		  !Plasma Temperature 
		  Temp_ion=((AmpB_Harris**2)*0.5_psn)/(epc*massi*c*c*(1.0_psn+TeTiratio)) !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
		  Temp_elc=Temp_ion*TeTiRatio*(mi/me)!Temp_elc=kT_e/m_ec^2
          !look up table for the PDF 
		  call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
		  call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
          !-------------------------------------------------
          ! Initialize Particles and Fields
          !-------------------------------------------------  
		  call InitBackgroundPrtl
		  call InitDriftingPrtl
		  call InitMagFld
		  !call RuptureSheet

end subroutine InitUser

subroutine Finalsubroutines 
     !define here all setup-specific actions that are required at the end of every time step
     call LoadBalanceYRecn ! Shift boundaries to balance the load and to imporove the performance
	 call SetFldBoundaries ! Field Boundary Condition
	 call ReflPrtlTop(real(TopEdge_Ypos-yborders(procyind(proc))+3,psn)) ! Particle reflecting Boundary condtition at the top edge 
	 call ReflPrtlBottom(real(BottomEdge_Ypos-yborders(procyind(proc))+3,psn)) ! Particle reflecting Boundary condtition at the bottom edge 
end subroutine Finalsubroutines

!================================================================================
! set-up specific subroutines
!================================================================================
subroutine RuptureSheet
	integer :: n
	real(psn) :: x1,x2,y1,y2
	x1=nx/2.0_psn-xborders(procxind(proc))+3
	x2=2.0+ nx/2.0_psn-xborders(procxind(proc))+3
	y1=ny/2.0_psn-yborders(procyind(proc))+3
	y2=2.0+ ny/2.0_psn-yborders(procyind(proc))+3
	
		do n=1,used_prtl_arr_size !particles
			if((xp(n).gt.x1).and.(xp(n).lt.x2).and.(yp(n).gt.y1).and.(yp(n).lt.y2)) call DeletePrtl(n)
		end do 
end subroutine RuptureSheet 
subroutine VelPrtb(x,y,z,AmpV,vx,vy,vz)
	real(psn) , intent(IN) :: x,y,z,AmpV
	real(psn) , intent(OUT):: vx,vy,vz
	real(psn) :: kmin,kmax,spec_ind,dlogk,ki,phase
	integer :: Nmodes,i
	kmin=2*pi/nx
	kmax=2*pi/c_ompe
	dlogk=0.2_psn
	Nmodes=(log10(kmax)-log10(kmin))/dlogk+1
    vx=0.0_psn
	vy=0.0_psn
	vz=0.0_psn
	do i=1,Nmodes
		ki=kmin*10.0_psn**((i-1)*dlogk)
		phase=2*pi*real(i)/real(Nmodes)
		vy=vy+AmpV*cos(ki*x+phase)
	end do  
end subroutine VelPrtb


subroutine SetFldBoundaries

	call ReflFldBoundaryTop(ny-6-yborders(procyind(proc))+3)
		
	call ReflFldBoundaryBottom(6-yborders(procyind(proc))+3)
		
end subroutine SetFldBoundaries

subroutine InitMagFld
	integer :: i,j,k
	real(psn) :: yglobal
#ifdef twoD
    do k=1,1
#else 	
	do k=1,mz
#endif 
       do j=1,my 
		   do i=1,mx
			  yglobal=j+0.5_psn+yborders(procyind(proc))-3
		      Bx(i,j,k)=-AmpB_Harris*tanh( (yglobal-(ny/2.0_psn)) /SheetThickness)
		   end do 
       end do  
    end do 
end subroutine InitMagFld

subroutine ReflFldBoundaryBottom(y1)
	integer :: y1
	if(y1.lt.1) return
	Ex(:,1:min(my,y1),:)=0.0_psn
	Ez(:,1:min(my,y1),:)=0.0_psn
	!By(:,1:min(my,y1),:)=0.0_psn	
end subroutine ReflFldBoundaryBottom
subroutine ReflFldBoundaryTop(y1)
	integer :: y1
	if(y1.gt.my) return
	Ex(:,max(1,y1):my,:)=0.0_psn
	Ez(:,max(1,y1):my,:)=0.0_psn
	!By(:,max(1,y1):my,:)=0.0_psn
end subroutine ReflFldBoundaryTop



subroutine InitBackgroundPrtl
	 implicit none 
     real(psn):: x1,x2,y1,y2
	 real(psn) :: MeanCount1
	 real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
	 real(psn), dimension(InjCountTableSize) :: InjCountTable1
 	 integer   :: Nelc_New,i
	 
	 	 
	 x1=max(real(LeftEdge_Xpos,psn),real(xborders(procxind(proc)),psn))
	 x2=min(real(RightEdge_Xpos,psn),real(xborders(procxind(proc)+1),psn)) 
	 x1=x1-xborders(procxind(proc))+3 !change to local cordinate 
	 x2=x2-xborders(procxind(proc))+3 
	 
	 y1=max(real(BottomEdge_Ypos,psn),real(yborders(procyind(proc)),psn))
	 y2=min(real(TopEdge_Ypos,psn),real(yborders(procyind(proc)+1),psn))	 
	 y1=y1-yborders(procyind(proc))+3 !change to local cordinate 
	 y2=y2-yborders(procyind(proc))+3 
	 	 
#ifdef twoD
	 MeanCount1=(x2-x1)*(y2-y1)*epc*BackgroundPrtlFraction
#else 
     MeanCount1=(x2-x1)*(y2-y1)*epc*(mz-5)*BackgroundPrtlFraction
#endif
     call InitPoissonDist(InjCountTableSize,InjCountTable1,MeanCount1)
	 call GetInjPrtlCount(MeanCount1,Nelc_New,InjCountTableSize,InjCountTable1)

	 do i=1,Nelc_New
		!Ions  
		call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,y1,y2,zmin,zmax)			
	    call GetVelGamma_MaxwellNonRel_StaticPDF(0.0_psn,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion)
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
		!Electrons
		call GetVelGamma_MaxwellNonRel_StaticPDF(0.0_psn,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc)		
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn)
	 end do 	 	 	 		  
end subroutine InitBackgroundPrtl

subroutine InitDriftingPrtl
	 implicit none 
     real(psn):: x1,x2,y1,y2
	 real(psn) :: MeanCount1
	 real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,r1,yglobal,xglobal,zglobal
	 real(psn), dimension(InjCountTableSize) :: InjCountTable1
 	 integer   :: Nelc_New,i
	 real(psn) :: den
	 real(psn) :: drifti,drifte,drift,drifti_tot,drifte_tot,drift_prtb
	 real(psn) :: vx,vy,vz
	 real(psn) :: vx_prtb,vy_prtb,vz_prtb,AmpPrtb
	 
	 	 
	 x1=max(real(LeftEdge_Xpos,psn),real(xborders(procxind(proc)),psn))
	 x2=min(real(RightEdge_Xpos,psn),real(xborders(procxind(proc)+1),psn)) 
	 x1=x1-xborders(procxind(proc))+3 !change to local cordinate 
	 x2=x2-xborders(procxind(proc))+3 
	 
	 y1=max(real(BottomEdge_Ypos,psn),real(yborders(procyind(proc)),psn))
	 y2=min(real(TopEdge_Ypos,psn),real(yborders(procyind(proc)+1),psn))	 
	 y1=y1-yborders(procyind(proc))+3 !change to local cordinate 
	 y2=y2-yborders(procyind(proc))+3 
	 	 
#ifdef twoD
     MeanCount1=(x2-x1)*(y2-y1)*epc
#else 
     MeanCount1=(x2-x1)*(y2-y1)*epc*(mz-5)
#endif

     call InitPoissonDist(InjCountTableSize,InjCountTable1,MeanCount1)
	 call GetInjPrtlCount(MeanCount1,Nelc_New,InjCountTableSize,InjCountTable1)
     
	 drift=AmpB_Harris/(qi*epc*SheetThickness)
	 drifti=drift*(1.0_psn/(TeTiratio+1.0_psn)) 
	 drifte=-drift*(TeTiratio/(TeTiratio+1.0_psn))
	 AmpPrtb=drift*0.01
	 
	 do i=1,Nelc_New  
		!Ions  
		call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1,x2,y1,y2,zmin,zmax)	
		call random_number(r1)
		xglobal=xlocal+xborders(procxind(proc))-3
		yglobal=ylocal+yborders(procyind(proc))-3
		zglobal=zlocal+zborders(proczind(proc))-3
		
		den=1.0_psn/(cosh((yglobal-(ny/2.0_psn))/SheetThickness))**2
		if(r1.gt.den) cycle

		!drift_prtb=0.01*drift*sin(4*pi*(xglobal-LeftEdge_Xpos)/(RightEdge_Xpos-LeftEdge_Xpos)) !perturb the velocity of particles
		call VelPrtb(xglobal,yglobal,zglobal,AmpPrtb,vx_prtb,vy_prtb,vz_prtb)
		vx=0.0_psn+vx_prtb
		vy=0.0_psn+vy_prtb
		vz=drifti +vz_prtb
		drifti_tot=sqrt(vx**2+vy**2+vz**2)
	    call GetVelGamma_MaxwellNonRel_StaticPDF(drifti_tot,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion,(/vx, vy, vz/))
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,0,1,0.0_psn)
		!Electrons
		vx=0.0_psn+vx_prtb
		vy=0.0_psn+vy_prtb
		vz=drifte +vz_prtb
		drifte_tot=sqrt(vx**2+vy**2+vz**2)
		call GetVelGamma_MaxwellNonRel_StaticPDF(drifte_tot,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc,(/vx, vy, vz/))
		call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,0,2,0.0_psn)
	 end do 	 	 	 		  
end subroutine InitDriftingPrtl

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


subroutine ReflPrtlTop(y1)
	real(psn) :: y1
	integer :: n
	if(y1.lt.my-1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).gt.y1)) then 
				 vp(n)=-vp(n)
                 yp(n)=y1-(yp(n)-y1)
			end if 
     	end do
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(ytp(n).gt.y1)) then 
				vtp(n)=-vtp(n)
				ytp(n)=y1-(ytp(n)-y1)
			end if 
		end do 
	end if
end subroutine ReflPrtlTop
subroutine ReflPrtlBottom(y1)
	real(psn) :: y1
	integer :: n
	if(y1.gt.1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).lt.y1)) then 
				vp(n)=-vp(n)
				yp(n)=y1+(y1-yp(n))
			end if 
		end do 
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(ytp(n).lt.y1)) then 
				 vtp(n)=-vtp(n)
				 ytp(n)=y1+(y1-ytp(n))
			 end if 
		end do 
	end if
end subroutine ReflPrtlBottom

     
!================================================================================
!The following subroutines are meant for customization and must be part of each setup module
!================================================================================
subroutine InitOverride!override the default intial conditons, such as initial fields, intial current, external magnetic field 
end subroutine InitOverride  

subroutine PostMovDep
	integer :: PrtlWall_local
	!copy current in the top layers 
	PrtlWall_local=TopEdge_Ypos-yborders(procyind(proc))+3
	if(PrtlWall_local.ge.2.and.PrtlWall_local.le.my-2) then
		Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local+1,:)
	 	Jy(:,PrtlWall_local-1,:)=Jy(:,PrtlWall_local-1,:)-Jy(:,PrtlWall_local,:)
	 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local+1,:)
		Jx(:,PrtlWall_local+1,:)=0.0_psn
		Jy(:,PrtlWall_local,:)=0.0_psn
		Jz(:,PrtlWall_local+1,:)=0.0_psn
    end if
	
	PrtlWall_local=BottomEdge_Ypos-yborders(procyind(proc))+3
	if(PrtlWall_local.ge.3.and.PrtlWall_local.le.my-2) then
		!To  ensure that the current is deposited on right place for reflected particles
		Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local-1,:)
	 	Jy(:,PrtlWall_local,:)=Jy(:,PrtlWall_local,:)-Jy(:,PrtlWall_local-1,:)
	 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local-1,:)
		Jx(:,PrtlWall_local-1,:)=0.0_psn
		Jy(:,PrtlWall_local-1,:)=0.0_psn
		Jz(:,PrtlWall_local-1,:)=0.0_psn
	end if
end subroutine PostMovDep
subroutine PreAddCurrent
end subroutine PreAddCurrent



end module setup