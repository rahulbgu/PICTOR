!Set up for magnetic reconnection 
module setup

use parameters
use vars
use help_setup
use memory
!use loadbalance
implicit none 

integer, parameter :: InjCountTableSize=1000
integer, parameter :: GammaTableLength=10000
real(psn) :: Temp_ion,Temp_elc,Alfven_speed,PlasmaBeta,DriftSpeed
real(psn) :: BextMag,Btheta,TeTiRatio
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
integer :: TopEdge_Ypos,BottomEdge_Ypos,LeftEdge_Xpos,RightEdge_Xpos
real(psn) :: ExTop,EyTop,EzTop,BxTop,ByTop,BzTop
real(psn) :: ExBottom,EyBottom,EzBottom,BxBottom,ByBottom,BzBottom
contains 
     
subroutine InitUser
          real(psn) :: ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,drift
          real(psn) :: xlocal,ylocal,zlocal
		  real(psn) :: r1,r2,r3
		  integer   :: i
        
          !-------------------------------------------------
          ! PARAMETERS for this setup 
          !-------------------------------------------------          
          Alfven_speed=0.5 !in units of c 
          DriftSpeed=0.05
		  PlasmaBeta=0.1!3
          TeTiRatio=1.0     
          BTheta=90*(pi/180) !Angle (radian) between magnetic field and z-axis  
          
		  
		  !Dependent Varaibles ---------------------
                           
          !set the background magentic field 
          !BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
          !Bz_ext0=BextMag*cos(Btheta)
          !Bx_ext0=BextMag*sin(Btheta)
          
		  !Plasma Temperature 
          Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
          Temp_elc=Temp_ion*TeTiRatio*(mi/me)!Temp_elc=kT_e/m_ec^2
		  
		  BxTop=0!Alfven_speed*sqrt(epc*(massi+masse)*c*c)
		  ByTop=0
		  BzTop=0
		  BxBottom=0!-BxTop
		  ByBottom=0
		  ByBottom=0
		  
		  ExTop=0
		  EyTop=0
		  EzTop=0!-DriftSpeed*BxTop
		  ExBottom=0
		  EyBottom=0
		  EzBottom=0!DriftSpeed*BxBottom
		  
  		  call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
          call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)
		            
!--------------------!!Initialize Particles and Fields!!-------------------------! 
TopEdge_Ypos=ny-12
BottomEdge_Ypos=12          

LeftEdge_Xpos=12
RightEdge_Xpos=nx-12

call InitPrtl


call SetFldBoundaries
call ClearPrtlOutsideBox
               
!-----------------------------------------------------------------------------!
end subroutine InitUser

!================================================================================
! set-up specific subroutines
!================================================================================

subroutine SetFldBoundaries
	!!!call AttenuateFldBottom(BottomEdge_Ypos-2-yborders(procyind(proc))+3)
	!!!call AttenuateFldBottom(TopEdge_Ypos+2-yborders(procyind(proc))+3)
	
	!!call SetFldTop(ny-3-yborders(procyind(proc))+3)
	!!call SetFldTop(ny-2-yborders(procyind(proc))+3)
	!!call SetFldTop(ny-1-yborders(procyind(proc))+3)
	call SetFldTop(ny-6-yborders(procyind(proc))+3)
	!call SetFldTop(TopEdge_Ypos-yborders(procyind(proc))+3)
	
	
	!!call SetFldBottom(3-yborders(procyind(proc))+3)
	!!call SetFldBottom(2-yborders(procyind(proc))+3)
	!!call SetFldBottom(1-yborders(procyind(proc))+3)
	call SetFldBottom(6-yborders(procyind(proc))+3)
	
	!call SetFldBottom(BottomEdge_Ypos-yborders(procyind(proc))+3)
	
	
! 	call CopyFldLeft(LeftEdge_Xpos-1-xborders(procxind(proc))+3)
! 	call CopyFldLeft(LeftEdge_Xpos-2-xborders(procxind(proc))+3)
! 	call CopyFldLeft(LeftEdge_Xpos-3-xborders(procxind(proc))+3)
! 	call CopyFldLeft(LeftEdge_Xpos-4-xborders(procxind(proc))+3)
! 	call CopyFldLeft(LeftEdge_Xpos-5-xborders(procxind(proc))+3)
! 	call CopyFldLeft(LeftEdge_Xpos-6-xborders(procxind(proc))+3)
! 	call SetFldLeftZero(LeftEdge_Xpos-7-xborders(procxind(proc))+3)
!
!
! 	call CopyFldRight(RightEdge_Xpos-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+1-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+2-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+3-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+4-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+5-xborders(procxind(proc))+3)
! 	call CopyFldRight(RightEdge_Xpos+6-xborders(procxind(proc))+3)
! 	call SetFldRightZero(RightEdge_Xpos+7-xborders(procxind(proc))+3)
	
	
	
	
end subroutine SetFldBoundaries
subroutine UpdateTopBottomDriftFld
	real(psn) :: Amp
	Amp=t/1000.0_psn ! linearly increase the value of amplitude
	Amp=min(1.0_psn,Amp)
	
	BxTop=Amp*Alfven_speed*sqrt(epc*(massi+masse)*c*c)
	BxBottom=-BxTop
	EzTop=-DriftSpeed*BxTop
    EzBottom=DriftSpeed*BxBottom
end subroutine UpdateTopBottomDriftFld

subroutine ClearPrtlOutsideBox
    call ClearPrtlTop(real(TopEdge_Ypos,psn)-yborders(procyind(proc))+3)
	call ClearPrtlBottom(real(BottomEdge_Ypos,psn)-yborders(procyind(proc))+3)
	call ClearPrtlLeft(real(LeftEdge_Xpos,psn)-0.5_psn-xborders(procxind(proc))+3)
	call ClearPrtlRight(real(RightEdge_Xpos,psn)+0.5_psn-xborders(procxind(proc))+3)
end subroutine ClearPrtlOutsideBox

subroutine SetFldTop(y1)
	integer :: y1
	
	if(y1.gt.my) return
	Ex(:,max(1,y1):my,:)=ExTop
	Ey(:,max(1,y1):my,:)=EyTop
	Ez(:,max(1,y1):my,:)=EzTop
	
	Bx(:,max(1,y1):my,:)=BxTop
	By(:,max(1,y1):my,:)=ByTop
	Bz(:,max(1,y1):my,:)=BzTop
end subroutine SetFldTop
subroutine SetFldBottom(y1)
	integer :: y1
	
	if(y1.lt.1) return
	    Ex(:,1:min(my,y1),:)=ExBottom
	    Ez(:,1:min(my,y1),:)=EzBottom
	    By(:,1:min(my,y1),:)=ByBottom
    
	if(y1.le.1) return !to make sure that y1-1 remain a positive integer 
	Ey(:,1:min(my,y1-1),:)=EyBottom
	Bx(:,1:min(my,y1-1),:)=BxBottom
	Bz(:,1:min(my,y1-1),:)=BzBottom
end subroutine SetFldBottom


subroutine AttenuateFldBottom(y2)
	integer :: y1,y2
	integer :: Atten_Ncell=8!number of cells for attenuation
	real(psn) :: AttenFactor=0.5_psn
    integer :: i,j,k,jmin,jmax
	
	y1=y2-Atten_Ncell 
	jmax=min(y2,my)
	jmin=max(1,y1)
    if(y2.ge.1) then
		 do k=1,mz
	          do j=jmin,jmax
	               do i=1,mx                   
	                    Ex(i,j,k)=(Ex(i,j,k)-ExBottom)*AttenFactor+ExBottom		
	                    Ey(i,j,k)=(Ey(i,j,k)-EyBottom)*AttenFactor+EyBottom
	                    Ez(i,j,k)=(Ez(i,j,k)-EzBottom)*AttenFactor+EzBottom   
	                    Bx(i,j,k)=(Bx(i,j,k)-BxBottom)*AttenFactor+BxBottom		
	                    By(i,j,k)=(By(i,j,k)-ByBottom)*AttenFactor+ByBottom
	                    Bz(i,j,k)=(Bz(i,j,k)-BzBottom)*AttenFactor+BzBottom  
	               end do
	          end do
	     end do
    end if		
end subroutine AttenuateFldBottom

subroutine AttenuateFldTop(y1)
	integer :: y1,y2
	integer :: Atten_Ncell=8!number of cells for attenuation
	real(psn) :: AttenFactor=0.5_psn
    integer :: i,j,k,jmin,jmax
	
	y2=y1+Atten_Ncell 
	jmax=min(y2,my)
	jmin=max(1,y1)
    if(y1.le.my) then
		 do k=1,mz
	          do j=jmin,jmax
	               do i=1,mx                   
	                    Ex(i,j,k)=(Ex(i,j,k)-ExTop)*AttenFactor+ExTop		
	                    Ey(i,j,k)=(Ey(i,j,k)-EyTop)*AttenFactor+EyTop
	                    Ez(i,j,k)=(Ez(i,j,k)-EzTop)*AttenFactor+EzTop   
	                    Bx(i,j,k)=(Bx(i,j,k)-BxTop)*AttenFactor+BxTop		
	                    By(i,j,k)=(By(i,j,k)-ByTop)*AttenFactor+ByTop
	                    Bz(i,j,k)=(Bz(i,j,k)-BzTop)*AttenFactor+BzTop  
	               end do
	          end do
	     end do
    end if		
end subroutine AttenuateFldTop




subroutine CopyFldRight(x1)
	integer :: x1
	if((x1.ge.1).and.x1.lt.mx) then 
	    Ey(x1+1,:,:)=Ey(x1,:,:)
	    Ez(x1+1,:,:)=Ez(x1,:,:)
	    Bx(x1+1,:,:)=Bx(x1,:,:)	
	end if
	if((x1.gt.1).and.x1.le.mx) then 
	    Ex(x1,:,:)=Ex(x1-1,:,:)
	    By(x1,:,:)=By(x1-1,:,:)
	    Bz(x1,:,:)=Bz(x1-1,:,:)
	end if 	
end subroutine CopyFldRight

subroutine CopyFldLeft(x1)
	integer :: x1
	if((x1.ge.1).and.(x1.lt.mx)) then 
	    Ey(x1,:,:)=Ey(x1+1,:,:)
	    Ez(x1,:,:)=Ez(x1+1,:,:)
	    Ex(x1,:,:)=Ex(x1+1,:,:)
	    By(x1,:,:)=By(x1+1,:,:)
	    Bz(x1,:,:)=Bz(x1+1,:,:)
	    Bx(x1,:,:)=Bx(x1+1,:,:)	
	end if 	
end subroutine CopyFldLeft
subroutine SetFldLeftZero(x1)
	integer :: x1	
	if(x1.ge.1) then 
	    Ey(1:min(mx,x1),:,:)=0.0_psn
	    Ez(1:min(mx,x1),:,:)=0.0_psn
	    Ex(1:min(mx,x1),:,:)=0.0_psn
	    By(1:min(mx,x1),:,:)=0.0_psn
	    Bz(1:min(mx,x1),:,:)=0.0_psn
	    Bx(1:min(mx,x1),:,:)=0.0_psn
	end if 
end subroutine SetFldLeftZero
subroutine SetFldRightZero(x1)
	integer :: x1
	if(x1.lt.mx) then
	    Ey(max(1,x1+1):mx,:,:)=0.0_psn
	    Ez(max(1,x1+1):mx,:,:)=0.0_psn
	    Bx(max(1,x1+1):mx,:,:)=0.0_psn	
	end if
	if(x1.le.mx) then 
	    Ex(max(1,x1):mx,:,:)=0.0_psn
	    By(max(1,x1):mx,:,:)=0.0_psn
	    Bz(max(1,x1):mx,:,:)=0.0_psn
	end if 	
end subroutine SetFldRightZero

subroutine InitPrtl
	real(psn):: y1,y2
    
	y1=max(real(BottomEdge_Ypos,psn),real(yborders(procyind(proc)),psn))
    y2=min(real(TopEdge_Ypos,psn),real(yborders(procyind(proc)+1),psn))
    
	y1=y1-yborders(procyind(proc))+3 !change to local cordinate 
    y2=y2-yborders(procyind(proc))+3 
	call HomogeneousDriftingPrtl(y1,y2,0.0_psn) 
end subroutine InitPrtl

subroutine InjectNewPrtl
	 implicit none 
	 integer  :: procyind_this
     real(psn):: y1,y2
	 integer  :: xrad_local
	 
	 procyind_this=procyind(proc)
	 
	 !----First Inject Particles at the Top Edge ----
	 
	 !x1 and x2 determined boundary of the region in local cordinates where new upstream plasma is injected 
	 y1=max(real(TopEdge_Ypos),real(yborders(procyind_this)))
	 y2=min(real(TopEdge_Ypos)+0.5,real(yborders(procyind_this+1)))
	 
	 y1=y1-yborders(procyind_this)+3 !change to local cordinate 
	 y2=y2-yborders(procyind_this)+3 
	 
	 call HomogeneousDriftingPrtl(y1,y2,-DriftSpeed) 
	 
	 
	 !----Now inject particles from the bottom edge ----
	 y1=max(real(BottomEdge_Ypos)-0.5,real(yborders(procyind_this)))
	 y2=min(real(BottomEdge_Ypos),real(yborders(procyind_this+1)))
	
	 y1=y1-yborders(procyind_this)+3 !change to local cordinate 
	 y2=y2-yborders(procyind_this)+3 
	 
	 call HomogeneousDriftingPrtl(y1,y2,DriftSpeed) 
	 	 		  
end subroutine InjectNewPrtl

!This subsourtine creates all particles of same weight
subroutine HomogeneousDriftingPrtl(y1,y2,drift)
	implicit none 
	real(psn) :: y1,y2,drift
	real(psn) :: dy
    real(psn) :: xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,r1,tag_fraction
    real(psn) :: ubeta,vbeta,wbeta,beta_drift
    real(psn) :: MeanCount1,MeanCount2   
	real(psn):: MeanCount1Prev=-1,MeanCount2Prev=-1
	real(psn), dimension(InjCountTableSize) :: InjCountTable1,InjCountTable2
	integer   :: Nelc_New,i,tag
	real(psn) :: time_this
	real(psn) :: x1,x2 
	
	save MeanCount1Prev,MeanCount2Prev,InjCountTable1,InjCountTable2
	
	tag_fraction=1.0_psn/psave_ratio
	CurrentTagID=1+NtagProcLen*proc
    TagCounter=1
	time_this=t
	
	x1=max(xmin,real(LeftEdge_Xpos,psn)-xborders(procxind(proc))+3)
	x2=min(xmax,real(RightEdge_Xpos,psn)-xborders(procxind(proc))+3)
	
	dy=max(0.0_psn,y2-y1)  
#ifdef twoD
   	MeanCount1=dy*epc*(mx-5)
#else 
	MeanCount1=dy*epc*(mx-5)*(mz-5)
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
			call GetRandomPositionHomogeneousInSubDomain(xlocal,ylocal,zlocal,x1+4,x2-4,y1,y2,zmin,zmax)			
		    call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaIon,Gamma_PDF_Ion,(/0.0_psn, 1.0_psn , 0.0_psn/))
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,time_this)
			!ions
			call GenerateTag(tag,2,tag_fraction)
			call GetVelGamma_MaxwellNonRel_StaticPDF(drift,ugamma,vgamma,wgamma,GammaTableLength,GammaElc,Gamma_PDF_Elc,(/0.0_psn, 1.0_psn , 0.0_psn/))		
			call InsertNewPrtl(xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,time_this)
		 end do
	end if	  
end subroutine HomogeneousDriftingPrtl
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

subroutine CopyPrtlLeft(x1,dx)
	real(psn) :: x1,dx
	integer :: n
	
	if(x1.gt.1) then
		do n=1,used_prtl_arr_size 
			if((qp(n).ne.0).and.(xp(n).gt.x1).and.(xp(n).lt.x1+dx)) then 
                call InsertNewPrtl(2*x1-xp(n),yp(n),zp(n),up(n),vp(n),wp(n),qp(n),tagp(n),flvp(n),var1p(n))
			end if
		end do 
	end if  
end subroutine CopyPrtlLeft
subroutine CopyPrtlRight(x1,dx)
	real(psn) :: x1,dx
	integer :: n
	if(x1.gt.1) then
		do n=1,used_prtl_arr_size 
			if((qp(n).ne.0).and.(xp(n).gt.x1).and.(xp(n).lt.x1+dx)) then 
                call InsertNewPrtl(2*x1-xp(n),yp(n),zp(n),up(n),vp(n),wp(n),qp(n),tagp(n),flvp(n),var1p(n))
			end if
		end do 
	end if  
end subroutine CopyPrtlRight

subroutine ClearPrtlTop(y1)
	real(psn) :: y1
	integer :: n
	if(y1.lt.my-1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).gt.y1)) call DeletePrtl(n)
     	end do 
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(ytp(n).gt.y1)) call DeleteTestPrtl(n)
		end do 
	end if
end subroutine ClearPrtlTop

subroutine ClearPrtlBottom(y1)
	real(psn) :: y1
	integer :: n
	if(y1.gt.1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).lt.y1)) call DeletePrtl(n)
		end do 
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(ytp(n).lt.y1)) call DeleteTestPrtl(n)
		end do 
	end if
end subroutine ClearPrtlBottom

subroutine ClearPrtlLeft(x1)
	real(psn) :: x1
	integer :: n
	if(x1.gt.1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(xp(n).lt.x1)) call DeletePrtl(n)
		end do 
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(xtp(n).lt.x1)) call DeleteTestPrtl(n)
		end do 
	end if
end subroutine ClearPrtlLeft
subroutine ClearPrtlRight(x1)
	real(psn) :: x1
	integer :: n
	if(x1.lt.mx-1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(xp(n).gt.x1)) call DeletePrtl(n)

		end do 
		do n=1,used_test_prtl_arr_size !test particles
			if((qtp(n).ne.0).and.(xtp(n).gt.x1)) call DeleteTestPrtl(n)
		end do 
	end if
end subroutine ClearPrtlRight

     
!================================================================================
!The following subroutines are meant for customization and must be part of each setup module
!================================================================================
subroutine InitOverride!override the default intial conditons, such as initial fields, intial current, external magnetic field 
end subroutine InitOverride  

subroutine Finalsubroutines !this subroutine is deleted, see shock_setup 
     !define here all additional actions that is to be perform at the end of every time step
     !call LoadBalanceHomogeneous
	 call SetFldBoundaries
	 call ClearPrtlOutsideBox
	 call InjectNewPrtl
	 call UpdateTopBottomDriftFld
end subroutine Finalsubroutines

subroutine PostMovDep
end subroutine PostMovDep
subroutine PreAddCurrent
end subroutine PreAddCurrent


! The following subroutine was being written to purturb harris current sheet, status: unfinished

! subroutine InitBfldPrtb ! Initialise Magnetic perturubation
! 	integer :: i,j,k
! 	real(psn) :: yglobal,xglobal,AmpB_prtb,lx,ly
!
! 	AmpB_prtb=0.1*AmpB_Harris ! Amplitude of the purturbation
!
! 	lx=2*(RightEdge_Xpos-LeftEdge_Xpos)
! 	ly=2*(TopEdge_Ypos-BottomEdge_Ypos)
! #ifdef twoD
!     do k=1,1
! #else
! 	do k=1,mz
! #endif
!        do j=1,my
! 		   do i=1,mx
!  			  xglobal=i+xborders(procxind(proc))-3
! 			  yglobal=j+0.5_psn+yborders(procyind(proc))-3
! 		      if(yglobal.ge.BottomEdge_Ypos.and.yglobal.le.TopEdge_Ypos) then
! 			      Bx(i,j,k)=Bx(i,j,k)+AmpB_prtb*cos(2*pi*(xglobal-LeftEdge_Xpos)/lx)*sin(2*pi*(yglobal-BottomEdge_Ypos)/ly)
! 		      end if
!
!
!  			  xglobal=i+0.5_psn+xborders(procxind(proc))-3
! 			  yglobal=j+yborders(procyind(proc))-3
! 			  if(yglobal.ge.BottomEdge_Ypos.and.yglobal.le.TopEdge_Ypos) then
! 			      By(i,j,k)=By(i,j,k)-AmpB_prtb*sin(2*pi*(xglobal-LeftEdge_Xpos)/lx)*cos(2*pi*(yglobal-BottomEdge_Ypos)/ly)
! 		      end if
! 		   end do
!        end do
!     end do
! end subroutine InitBfldPrtb

end module setup