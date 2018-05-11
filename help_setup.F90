module help_setup
     use parameters
     use vars
	 use communication
     use memory
	 use prtl_stats 
     real(psn), dimension(:), allocatable :: flvrqmTemp
     integer, dimension(:), allocatable   :: FlvrSaveFldDataTemp, FlvrTypeTemp,FlvrSplitPrtlTemp,CurrentTagIDTemp,TagCounterTemp
contains 
!------------------------------------------------------------
! subroutines to set values associated with type of particles in case of multiple species 
!------------------------------------------------------------      
     subroutine SetQbyM(ind,value) ! Update this soubroutine such that the Flvrs can be added in any order, no incremently 
          integer :: ind
          real(psn) :: value 
        if(ind.gt.Nflvr) then 
               allocate(flvrqmTemp(ind),FlvrSaveFldDataTemp(ind),FlvrTypeTemp(ind),FlvrSplitPrtlTemp(ind),CurrentTagIDTemp(ind),TagCounterTemp(ind))
               flvrqmTemp(1:Nflvr)=flvrqm(1:Nflvr)
               FlvrSaveFldDataTemp(1:Nflvr)=FlvrSaveFldData(1:Nflvr)     
               FlvrTypeTemp(1:Nflvr)=FlvrType(1:Nflvr)
               FlvrSplitPrtlTemp(1:Nflvr)=FlvrSplitPrtl(1:Nflvr)
               CurrentTagIDTemp(1:Nflvr)=CurrentTagID(1:Nflvr)
               TagCounterTemp(1:Nflvr)=TagCounterTemp(1:Nflvr)
			   
               deallocate(flvrqm,FlvrSaveFldData,FlvrType,FlvrSplitPrtl,CurrentTagID,TagCounter)
               call move_alloc(flvrqmTemp,flvrqm)
               call move_alloc(FlvrSaveFldDataTemp,FlvrSaveFldData)
               call move_alloc(FlvrTypeTemp,FlvrType)
               call move_alloc(FlvrSplitPrtlTemp,FlvrSplitPrtl)
               call move_alloc(CurrentTagIDTemp,CurrentTagID)
               call move_alloc(TagCounterTemp,TagCounter)
               Nflvr=Nflvr+1
          end if
               flvrqm(ind)=value
     end subroutine SetQbyM
     subroutine SetFlvrSaveFldData(ind,value)
          integer :: ind
          integer :: value 
               FlvrSaveFldData(ind)=value
     end subroutine SetFlvrSaveFldData
     subroutine SetFlvrType(ind,value)
          integer :: ind
          integer :: value 
          FlvrType(ind)=value
     end subroutine SetFlvrType
	 
	 subroutine InitNewFlvr(FlvID,QbyM,Type,SaveFld,Split)
		 integer :: FlvID,Type,SaveFld,Split
		 real(psn) :: QbyM
		 call SetQbyM(FlvID,QbyM)
		 FlvrType(FlvID)=Type
		 FlvrSaveFldData(FlvID)=SaveFld
		 FlvrSplitPrtl(FlvID)=Split
		 CurrentTagID(FlvID)=1+proc*NtagProcLen 
		 TagCounter(FlvID)=1
	 end subroutine InitNewFlvr
	 
!------------------------------------------------------------
!subroutine to genrate unique tag for particles, warning:: tagging start from 1 if simulation restarted
!------------------------------------------------------------
subroutine GetTagID(tag,FlvID,tag_ratio)
	implicit none
	integer :: tag,FlvID,tag_ratio_this,offset_this
	integer, optional :: tag_ratio
	
    if(present(tag_ratio)) then 
		tag_ratio_this=tag_ratio
	else 
		tag_ratio_this=psave_ratio
	end if
	
	if(TagCounter(FlvID).gt.tag_ratio_this) TagCounter(FlvID)=1
	if(TagCounter(FlvID).eq.1) then  
		tag=CurrentTagID(FlvID)
	    if(mod(CurrentTagID(FlvID),NtagProcLen).eq.0) CurrentTagID(FlvID)=CurrentTagID(FlvID)+NtagProcLen*(nproc-1)
        CurrentTagID(FlvID)=CurrentTagID(FlvID)+1
    else
		tag=0
	end if	
	TagCounter(FlvID)=TagCounter(FlvID)+1
end subroutine GetTagID


!The following subroutine used to insert new particles in the simulation which are uniformly distributed over the entire box
! subroutine InsertMultipleParticles_ThermalAndHomogeneous(Nprtl,Temp,TagRatio,FlvID,Type,QbyM,SaveFldData)
!      integer,   intent(IN):: Nprtl,TagRatio,FlvID,Type,SaveFldData
!      real(psn), intent(IN):: Temp,QbyM
!      real(xpsn)           :: xlocal
!      real(ypsn)           :: ylocal
!      real(zpsn)           :: zlocal
!      real(psn)            :: ugamma,vgamma,wgamma
!      integer              :: Inserted
!      integer              :: InsertAtInd
!      integer              :: ptagID
!      save InsertAtInd
!      data InsertAtInd /1/
!      Inserted=0
!      CurrentTagID(FlvID)=1+proc*NtagProcLen
!      do while(Inserted.lt.Nprtl)
!
!           if(InsertAtInd.gt.prtl_arr_size) then !Make sure that prtl Arr is large enough
!                allocate(ptemp(prtl_arr_size+2*Nprtl))
!                ptemp(1:prtl_arr_size)=p(1:prtl_arr_size)
!                call move_alloc(ptemp,p)
!                p(prtl_arr_size+1:prtl_arr_size+2*Nprtl)%q=0
!                prtl_arr_size=prtl_arr_size+2*Nprtl
!           end if
!
!           if(p(InsertAtInd)%q.ne.0) then
!                InsertAtInd=InsertAtInd+1 ! the slot is not empty, try the next one
!           else
!                 Inserted=Inserted+1
!                 if(mod(Inserted,TagRatio).eq.0) then
!                     ptagID=CurrentTagID(FlvID)
!  				    if(mod(CurrentTagID(FlvID),NtagProcLen).eq.0) CurrentTagID(FlvID)=CurrentTagID(FlvID)+NtagProcLen*(nproc-1)
!  			        CurrentTagID(FlvID)=CurrentTagID(FlvID)+1
!                 else
!                      ptagID=0
!                 end if
!                call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)
!                call GetVelGamma_MaxwellJuttner(0.0_psn,Temp,ugamma,vgamma,wgamma)
!                call InsertParticleAt(InsertAtInd,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,ptagID,FlvID)
!           end if
!      end do
!      call SetQbyM(FlvID,QbyM)
!      call SetFlvrSaveFldData(FlvID,SaveFldData)
!      call SetFlvrType(FlvID,Type)
!      FlvrSplitPrtl(FlvID)=0
! end subroutine InsertMultipleParticles_ThermalAndHomogeneous


!------------------------------------------------------------------------------------------------
! The following subroutine return postion of a particle assuming that all particles are distributed uniformly in the simulation box 
!------------------------------------------------------------------------------------------------

subroutine GetRandomPositionHomogeneous(x,y,z)
     real(xpsn), intent(inout):: x
     real(ypsn), intent(inout):: y
     real(zpsn), intent(inout):: z
     real(psn) :: r1,r2,r3
     
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      
      x=xmin+(xmax-xmin)*r1
      y=ymin+(ymax-ymin)*r2
      z=zmin+(zmax-zmin)*r3
     
end subroutine GetRandomPositionHomogeneous 

subroutine GetRandomPositionHomogeneousInSubDomain(x,y,z,xi,xf,yi,yf,zi,zf)
     real(xpsn), intent(inout):: x
     real(ypsn), intent(inout):: y
     real(zpsn), intent(inout):: z
	 real(psn), intent(in) ::xi,xf,yi,yf,zi,zf
     real(psn) :: r1,r2,r3
     
      call random_number(r1)
      call random_number(r2)
      call random_number(r3)
      
      x=xi+(xf-xi)*r1
      y=yi+(yf-yi)*r2
      z=zi+(zf-zi)*r3
     
end subroutine GetRandomPositionHomogeneousInSubDomain 

! subroutine ReInitPrtlArr(new_size)
! 	integer :: new_size
!     prtl_arr_size=new_size
!     if(allocated(p)) deallocate(p)
! 	allocate(p(prtl_arr_size)) !allocate memory
!     p(1:prtl_arr_size)%q=0 !declare all slots empty at the begining
! end subroutine ReInitPrtlArr

!--------------------------------------------------------------------------------------------
!
!        Particle Splitting Algorithms 
!
!--------------------------------------------------------------------------------------------
! subroutine SplitHighEnergyTailSpeedSpec
!      integer, dimension(Nflvr) :: SplitVminInd,VmaxInd
!      integer :: i,n,bin_ind
!      real(psn) ::speed
!      !First determine the velocity range of the particles which are to be split
!      do i=1,Nflvr
!           call GetVelSpecIndRangeToSplit(i,SplitVminInd(i),VmaxInd(i))
!           SplitVminInd=SplitVminInd+(VmaxInd-SplitVminInd)/2
!      end do
!
!      do n=1,prtl_arr_size
!           if(abs(qp(n)).le.SplitQ_MIN) cycle
!           if(FlvrSplitPrtl(flvp(n)).eq.0) cycle
!           speed=sqrt(1.0_psn - 1.0_psn/(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n)))
!           bin_ind=floor(speed/Speed_spec_binwidth) +1
!           if(bin_ind.gt.SplitVminInd(flvp(n))) then
!                qp(n)=qp(n)*0.5_psn
!                if(pfree_ind.eq.0) call ReshapePrtlArr(int(1.2*np)+100) ! increase the length of prtl array if all slots are exhausted
!                free_slot=pfree(pfree_ind)
!                p(free_slot)=p(n)
!                if(tagp(n).gt.0) then !Now tag the new split particle for tracking
!                     p(free_slot)%a=CurrentTagID(p(free_slot)%flv)
!  				    if(mod(CurrentTagID(p(free_slot)%flv),NtagProcLen).eq.0) CurrentTagID(p(free_slot)%flv)=CurrentTagID(p(free_slot)%flv)+NtagProcLen*(nproc-1)
!                     CurrentTagID(p(free_slot)%flv)=CurrentTagID(p(free_slot)%flv)+1
!                end if
!                call ScatterParticleDoublet(n,free_slot,0.01_psn)
!                pfree_ind=pfree_ind-1
!                np=np+1
!           end if
!      end do
! end subroutine SplitHighEnergyTailSpeedSpec
!
! subroutine GetVelSpecIndRangeToSplit(FlvID,ind1,ind2)
!      integer :: FlvID,ind1,ind2,i
!      real(dbpsn):: MaxWeight
!      MaxWeight=0.0_dbpsn
!      do i=1,Speed_spec_binlen
!           if(spec_speed(FlvID,i).gt.MaxWeight) then
!                MaxWeight=spec_speed(FlvID,i)
!                ind1=i
!           end if
!           if(spec_speed(FlvID,i).gt.0) ind2=i
!      end do
! end subroutine GetVelSpecIndRangeToSplit
!
! subroutine ScatterParticleDoublet(ind1,ind2,delp_p)
!      implicit none
!      integer :: ind1,ind2
!      real(psn) :: delp_p
!      real(psn)::gbeta
!      real(psn)::Dtheta,CosTheta,SinTheta
!      real(psn)::Phi,CosPhi,SinPhi
!      real(psn)::CosDtheta,SinDtheta,PhiRandom,r1
!      real(psn)::u1,v1,w1,u2,v2,w2
!
!      gbeta=sqrt(p(ind1)%u**2+p(ind1)%v**2+p(ind1)%w**2) !magnitude of gamma*beta
!      !Theta and Phi determine direction of 3-velocity
!      CosTheta=p(ind1)%w/gbeta
!      SinTheta=sqrt(1-CosTheta*CosTheta)
!      Phi=atan2(p(ind1)%v,p(ind1)%u)
!      CosPhi=cos(Phi)
!      SinPhi=sin(Phi)
!
!      Dtheta=atan(delp_p)
!      SinDtheta=sin(Dtheta)
!      CosDtheta=cos(Dtheta)
!      call random_number(r1)
!      PhiRandom=2*pi*r1
!
!      w1=gbeta*CosDtheta
!      u1=gbeta*SinDtheta*cos(PhiRandom)
!      v1=gbeta*SinDtheta*sin(PhiRandom)
!
!      w2=gbeta*CosDtheta
!      u2=gbeta*SinDtheta*cos(PhiRandom+pi)
!      v2=gbeta*SinDtheta*sin(PhiRandom+pi)
!
!      u(ind1)=(u1*CosTheta +w1*SinTheta)*CosPhi-v1*SinPhi
!      v(ind1)=(u1*CosTheta +w1*SinTheta)*SinPhi+v1*CosPhi
!      w(ind1)=-u1*SinTheta +w1*CosTheta
!
!      u(ind2)=(u2*CosTheta +w2*SinTheta)*CosPhi-v2*SinPhi
!      v(ind2)=(u2*CosTheta +w2*SinTheta)*SinPhi+v2*CosPhi
!      w(ind2)=-u2*SinTheta +w2*CosTheta
!
!      !call RotateVector(u1,v1,w1,p(ind1)%u,p(ind1)%v,p(ind1)%w,CosTheta,Phi)
!      !call RotateVector(u2,v2,w2,p(ind2)%u,p(ind2)%v,p(ind2)%w,CosTheta,Phi)
!
! end subroutine ScatterParticleDoublet
subroutine SmallAngleScattering(u,v,w,ThetaMean)
	real(psn) :: u,v,w,ThetaMean
	real(psn) :: theta,phi
	real(psn) :: u1,v1,w1,u2,v2,w2,u3,v3,w3,p1,p2,p3
	real(psn) :: mag
	real(psn) :: r1,r2 
	
	call random_number(r2)
	theta=ThetaMean
	phi=2*pi*r2 
	
	mag=sqrt(u**2+v**2+w**2)
    p1=mag*cos(theta)
	p2=mag*sin(theta)*cos(phi) !New vector  
	p3=mag*sin(theta)*sin(phi)

	
	mag=sqrt(u**2+v**2+w**2)
	if(mag.eq.0) return
	u1=u/mag
	v1=v/mag
	w1=w/mag
	
    u2=-v1 
    v2=u1
    w2=0.0_psn
	mag=sqrt(u2**2+v2**2+w2**2)
    if(mag.ne.0) then 
	    u2=u2/mag 
	    v2=v2/mag
	else 
		u2=1.0_psn
		v2=0.0_psn
	end if 

    u3=v2*w1-w2*v1
    v3=w2*u1-u2*w1
    w3=u2*v1-v2*u1
	mag=sqrt(u3**2+v3**2+w3**2)
	u3=u3/mag
	v3=v3/mag
	w3=w3/mag

	u= u1*p1 + u2*p2 + u3*p3
	v= v1*p1 + v2*p2 + v3*p3
	w= w1*p1 + w2*p2 + w3*p3		
end subroutine SmallAngleScattering



! subroutine RotateVector(u,v,w,CosTheta,Phi)
!      real(psn) :: ui,vi,wi,uf,vf,wf
!      real(psn) :: CosTheta,Phi
!      real(psn) :: SinTheta,CosPhi,SinPhi
!
!      SinTheta=sqrt(1-CosTheta*CosTheta)
!      CosPhi=cos(Phi)
!      SinPhi=sin(Phi)
!
!      !first rotation about y-axis
!      vx2= vx1*vCosTheta +wi*vSinTheta
!      vy2= vy1
!      vz2=-vx1*vSinTheta +wi*vCosTheta
!
!      !now rotation about z-axis
!      uf=(vx1*vCosTheta +wi*vSinTheta)*CosPhi-vy1*SinPhi
!      vf=(vx1*vCosTheta +wi*vSinTheta)*SinPhi+vy1*CosPhi
!      wf=-ui*SinTheta +wi*CosTheta
! end subroutine RotateVector

subroutine ReflectPrtlMovingFrameX(Drift,ugamma,vgamma,wgamma)
	real(psn) :: Drift,DriftGamma,DriftBeta
	real(psn) :: ugamma,vgamma,wgamma 
	real(psn) :: gamma
    
	if(abs(Drift).lt.1.0) then 
         DriftBeta=Drift
         DriftGamma=1.0_psn/sqrt((1.0_psn-Drift)*(1.0_psn+Drift))          
    else 
         DriftGamma=abs(Drift)
         DriftBeta=sqrt((Drift-1.0_psn)*(Drift+1.0_psn))/Drift          
    end if
	
	gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
	DriftBeta=-DriftBeta
	ugamma=DriftGamma*ugamma + DriftGamma*DriftBeta*gamma!in the drifting frame 
	
	ugamma=-ugamma
	gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
	DriftBeta=-DriftBeta
	ugamma=DriftGamma*ugamma + DriftGamma*DriftBeta*gamma! back to the simualtion frame	
end subroutine ReflectPrtlMovingFrameX

subroutine HardScattering(ugamma,vgamma,wgamma)
	real(psn) :: ugamma,vgamma,wgamma
	real(psn) :: CosTheta,SinTheta,phi,mom
    real(psn) :: r1,r2
    call random_number(r1)
    call random_number(r2)
	CosTheta=2*r1-1
	SinTheta=sin(acos(CosTheta))
	phi=2*pi*r2
	mom=sqrt(ugamma**2+vgamma**2+wgamma**2)
    ugamma=mom*CosTheta 
	vgamma=mom*SinTheta*cos(phi)
	wgamma=mom*SinTheta*sin(phi)	
end subroutine HardScattering 

subroutine HardScatteringMovingFrame(ugamma,vgamma,wgamma,vx,vy,vz)
	real(psn) :: ugamma,vgamma,wgamma
	real(psn) :: vx,vy,vz
	
end subroutine HardScatteringMovingFrame 


! subroutine SplitPrtlInSubDomain_Rel(xi,xf,yi,yf,zi,zf)
! 	integer :: xi,xf,yi,yf,zi,zf
!     integer   :: x1,x2,y1,y2,z1,z2
! 	real(psn), dimension(:,:), allocatable :: NormalisedSpec
!     real(dbpsn) :: spec_peak_value
! 	integer, dimension(Nflvr) :: SpecPeakInd
! 	integer :: i,n,bin_ind
! 	real(psn) :: gamma,r1
!
! 	call CalcGmaxGminLocalInSubDomain(xi,xf,yi,yf,zi,zf)
!     call GetGmaxGminGlobal
!     call CreateGammaSpecBin
!     call CalcGammaSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
!     call AllReduceGammaSpectrum
!     call DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
!
!     if(.not.allocated(NormalisedSpec)) allocate(NormalisedSpec(Nflvr,0))
! 	if(size(NormalisedSpec,2).lt.size(spec_gamma,2)) then
! 		deallocate(NormalisedSpec)
! 		allocate(NormalisedSpec(Nflvr,size(spec_gamma,2)))
! 	end if
!     !Locate the peak of the gamma spectrum
!     SpecPeakInd=1
! 	do n=1,Nflvr
! 		spec_peak_value=0
! 	    do i=1,Gamma_spec_binlen
! 			if(spec_gamma(n,i).gt.spec_peak_value) then
! 				spec_peak_value=spec_gamma(n,i)
! 				SpecPeakInd(n)=i
! 			end if
! 		end do
! 	end do
! 	!Spectrum is normalised such that the value at peak=0 and at max_gamma=1
! 	NormalisedSpec=0*spec_gamma ! it determines the splitting prob. of a particle in the tail
! 	do n=1,Nflvr
! 	    do i=1,Gamma_spec_binlen
! 			if(i.gt.SpecPeakInd(n)) then
! 				NormalisedSpec(n,i)=(spec_gamma(n,SpecPeakInd(n))-spec_gamma(n,i))/spec_gamma(n,SpecPeakInd(n))
! 			end if
! 		end do
! 	end do
!
! 	do n=1,prtl_arr_size
! 	     if(abs(qp(n)).le.SplitQ_MIN) cycle
!   	     if((xp(n).lt.x1.or.xp(n).gt.x2).or.(yp(n).lt.y1.or.yp(n).gt.y2).or.(zp(n).lt.z1.or.zp(n).gt.z2)) cycle
!          gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
!          bin_ind=floor((log(gamma)-log(gmin_allflv))/Gamma_spec_binwidth)+1
!          if(bin_ind.gt.SpecPeakInd(flvp(n))) then
! 	          call random_number(r1)
! 			  if(r1.gt.NormalisedSpec(flvp(n),bin_ind)) then
! 				  qp(n)=qp(n)*0.5_psn
! 	              if(pfree_ind.eq.0) call ReshapePrtlArr(int(1.2*np)+100) ! increase the length of prtl array if all slots are exhausted
! 	              free_slot=pfree(pfree_ind)
! 	              p(free_slot)=p(n)
! 	              call ScatterParticleDoublet(n,free_slot,0.01_psn)
! 	              pfree_ind=pfree_ind-1
! 	              np=np+1
! 			  end if
! 		  end if
! 	end do
!
! end subroutine SplitPrtlInSubDomain_Rel
!
! subroutine SplitPrtlInSubDomain_NonRel(xi,xf,yi,yf,zi,zf)
! 	integer :: xi,xf,yi,yf,zi,zf
!     integer   :: x1,x2,y1,y2,z1,z2
! 	real(psn), dimension(:,:), allocatable :: NormalisedSpec
!     real(dbpsn) :: spec_peak_value
! 	integer, dimension(Nflvr) :: SpecPeakInd
! 	integer :: i,n,bin_ind
! 	real(psn) :: speed,r1
!
!     call CalcSpeedSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
!     call AllReduceSpeedSpectrum
!     call DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
!
!     if(.not.allocated(NormalisedSpec)) allocate(NormalisedSpec(Nflvr,Speed_spec_binlen))
!     !Locate the peak of the gamma spectrum
!     SpecPeakInd=1
! 	do n=1,Nflvr
! 		spec_peak_value=0
! 	    do i=1,Speed_spec_binlen
! 			if(spec_speed(n,i).gt.spec_peak_value) then
! 				spec_peak_value=spec_speed(n,i)
! 				SpecPeakInd(n)=i
! 			end if
! 		end do
! 	end do
! 	!Spectrum is normalised such that the value at peak=0 and at max_gamma=1
! 	NormalisedSpec=0*spec_speed ! it determines the splitting prob. of a particle in the tail
! 	do n=1,Nflvr
! 	    do i=1,Speed_spec_binlen
! 			NormalisedSpec(n,i)=(spec_speed(n,SpecPeakInd(n))-spec_speed(n,i))/spec_speed(n,SpecPeakInd(n))
! 		end do
! 	end do
!
! 	do n=1,prtl_arr_size
! 	     if(abs(qp(n)).le.SplitQ_MIN) cycle
!   	     if((xp(n).lt.x1.or.xp(n).gt.x2).or.(yp(n).lt.y1.or.yp(n).gt.y2).or.(zp(n).lt.z1.or.zp(n).gt.z2)) cycle
!          speed=sqrt(1.0_psn - 1.0_psn/(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n)))
!          bin_ind=floor(speed/Speed_spec_binwidth) +1
!          if(bin_ind.gt.SpecPeakInd(flvp(n))) then
! 	          call random_number(r1)
! 			  if(r1.gt.NormalisedSpec(flvp(n),bin_ind)) then
! 				  qp(n)=qp(n)*0.5_psn
! 	              if(pfree_ind.eq.0) call ReshapePrtlArr(int(1.2*np)+100) ! increase the length of prtl array if all slots are exhausted
! 	              free_slot=pfree(pfree_ind)
! 	              p(free_slot)=p(n)
! 	              call ScatterParticleDoublet(n,free_slot,1.0_psn)
! 	              pfree_ind=pfree_ind-1
! 	              np=np+1
! 			  end if
! 		  end if
! 	end do
!
! end subroutine SplitPrtlInSubDomain_NonRel



!========================================================================================================
!
!      Momentum distribution functions 
!
!========================================================================================================


!--------------------------------------------------------------------------------------------------------
!Distribution of Gamma of relativistic particles with drift
!a) Drift is along the x-axis (if < 1 it is beta, otherwise it is Lorentz factor Gamma)
!b) Here Temperature is in units of mc^2/k, i.e. Temp = kT/mc^2 , where T is the actual physical Temperature
!--------------------------------------------------------------------------------------------------------

subroutine GetVelGamma_MaxwellJuttner(Drift,Temp,ugamma,vgamma,wgamma,direction)     
     real(psn), intent(in)  :: Drift,Temp
     real(psn), dimension(3), intent(in), optional :: direction
     real(psn), intent(out) :: ugamma,vgamma,wgamma !these are gamma*beta
     real(psn)              :: Temp0=-1.0 ! values in the last call  
     integer, parameter :: GammaTableSize=10000 
     real(psn), save, dimension(GammaTableSize):: Gamma_Table,PDF_Table    
     real(psn) :: r1,gamma,beta
     integer :: index
	 
	 save Temp0

     if(Temp0.ne.Temp) then 
          call InitMaxwellJuttnerTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
          Temp0=Temp
     end if
          
     call random_number(r1)     
     call BinarySearch(GammaTableSize,PDF_Table,r1,index)   
	 if(PDF_Table(index+1).eq.PDF_Table(index)) then   
		 gamma=Gamma_Table(index)
	 else 
         gamma=Gamma_Table(index)+ ( (Gamma_Table(index+1)-Gamma_Table(index))/(PDF_Table(index+1)-PDF_Table(index)) )*(r1-PDF_Table(index))
     end if 
	 beta=sqrt((gamma-1)*(gamma+1))/gamma
     
	 call AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)    
end subroutine GetVelGamma_MaxwellJuttner

!-------------------------------------------------------------------------------------------
! Following subroutine requires PDF table to be suppplied by the user, useful in case of shock problem  
! Currently drift is only in the x-direction 
!-------------------------------------------------------------------------------------------
subroutine GetVelGamma_MaxwellJuttner_StaticPDF(Drift,ugamma,vgamma,wgamma,GammaTableSizeStatic,Gamma_TableStatic,PDF_TableStatic,direction)     
     real(psn), intent(in)  :: Drift
     real(psn), dimension(3), intent(in), optional :: direction
     real(psn), intent(out) :: ugamma,vgamma,wgamma !these are gamma*beta
     integer                :: GammaTableSizeStatic
     real(psn), dimension(GammaTableSizeStatic):: Gamma_TableStatic,PDF_TableStatic   
     real(psn) :: r1,gamma,beta
     integer :: index
 
     call random_number(r1)     
     call BinarySearch(GammaTableSizeStatic,PDF_TableStatic,r1,index)  
	 if(PDF_TableStatic(index+1).eq.PDF_TableStatic(index)) then 
		  gamma=Gamma_TableStatic(index)
	 else 	    
          gamma=Gamma_TableStatic(index)+ ( (Gamma_TableStatic(index+1)-Gamma_TableStatic(index))&
	                                     /(PDF_TableStatic(index+1)-PDF_TableStatic(index)) )*(r1-PDF_TableStatic(index))
	 end if									   
	 beta=sqrt((gamma-1.0_psn)*(gamma+1.0_psn))/gamma
	 call AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)
end subroutine GetVelGamma_MaxwellJuttner_StaticPDF

!--------------------------------------------------------------------------------------------------------
!Distribution of v/c of non-rel particles with drift
! GammaTable and Gamma_PDF_Table in this case should be understood to store velocity like in VelTable and Vel_PDF_Table
!--------------------------------------------------------------------------------------------------------
subroutine GetVelGamma_MaxwellNonRel(Drift,Temp,ugamma,vgamma,wgamma,direction)     
     real(psn), intent(in)  :: Drift,Temp
     real(psn), dimension(3), intent(in), optional :: direction
     real(psn), intent(out) :: ugamma,vgamma,wgamma !these are gamma*beta
     real(psn)              :: Temp0=-1.0 ! values in the last call  
     integer, parameter :: GammaTableSize=10000 
     real(psn), save, dimension(GammaTableSize):: Gamma_Table,PDF_Table    
     real(psn) :: r1,gamma,beta
     integer :: index

     if(Temp0.ne.Temp) then 
          call InitMaxwellNonRelTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
          Temp0=Temp
     end if
          
     call random_number(r1)     
     call BinarySearch(GammaTableSize,PDF_Table,r1,index)    
	 if(PDF_Table(index+1).eq.PDF_Table(index)) then 
		 beta=Gamma_Table(index)
	 else 
	     beta=Gamma_Table(index)+ ( (Gamma_Table(index+1)-Gamma_Table(index))/(PDF_Table(index+1)-PDF_Table(index)) )*(r1-PDF_Table(index))
	 end if 
	 if(beta.ge.1.0) beta=0.0 !in some very rare cases beta can be equal to 1
	 
     gamma=1.0_psn/sqrt((1.0_psn+beta)*(1.0_psn-beta))
	 call AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)    
end subroutine GetVelGamma_MaxwellNonRel

subroutine GetVelGamma_MaxwellNonRel_StaticPDF(Drift,ugamma,vgamma,wgamma,GammaTableSizeStatic,Gamma_TableStatic,PDF_TableStatic,direction)     
     real(psn), intent(in)  :: Drift
     real(psn), dimension(3), intent(in), optional :: direction
     real(psn), intent(out) :: ugamma,vgamma,wgamma !these are gamma*beta
     integer                :: GammaTableSizeStatic
     real(psn), dimension(GammaTableSizeStatic):: Gamma_TableStatic,PDF_TableStatic   
     real(psn) :: r1,gamma,beta
     integer :: index
     call random_number(r1)     
     call BinarySearch(GammaTableSizeStatic,PDF_TableStatic,r1,index)  
	 if(PDF_TableStatic(index+1).eq.PDF_TableStatic(index)) then 
		 beta=Gamma_TableStatic(index)
	 else    
         beta=Gamma_TableStatic(index)+ ( (Gamma_TableStatic(index+1)-Gamma_TableStatic(index))&
	                                   /(PDF_TableStatic(index+1)-PDF_TableStatic(index)) )*(r1-PDF_TableStatic(index))
	 end if			

     if(beta.ge.1.0) beta=0.0 !in some very rare cases beta can be equal to 1
	       									   
	 gamma=1.0_psn/sqrt((1.0_psn+beta)*(1.0_psn-beta))
	 !print*,'gamma',gamma,beta,Gamma_TableStatic(index+1),Gamma_TableStatic(index),index,r1,PDF_TableStatic(index+1),PDF_TableStatic(index)
	 call AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)
end subroutine GetVelGamma_MaxwellNonRel_StaticPDF



!the following subroutine adds drift speed to the thermal speed 
subroutine AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)
    real(psn), intent(in)  :: Drift
	real(psn),intent(inout):: ugamma,vgamma,wgamma
	real(psn),intent(in)   :: gamma,beta
    real(psn), dimension(3), intent(in), optional :: direction
    real(psn) :: DriftBeta,DriftGamma
    real(psn) :: phi,CosTheta,SinTheta
    real(psn) :: r2,r3
    !variables used in the case when boost is not along the x-axis 
    real(psn), dimension(3) :: dirn
    real(psn) :: dirn_mag,projn
	
    if(abs(Drift).lt.1) then 
         DriftBeta=Drift
         DriftGamma=1/sqrt((1-Drift)*(1+Drift))          
    else 
         DriftGamma=abs(Drift)
         DriftBeta=sqrt((Drift-1)*(Drift+1))/Drift          
    end if
	
    call random_number(r2)
    call random_number(r3)
    CosTheta=2*r2-1 
    SinTheta=sqrt(1-CosTheta**2)
    phi=2*pi*r3
    
    ugamma=gamma*beta*CosTheta
    vgamma=gamma*beta*SinTheta*cos(phi)
    wgamma=gamma*beta*SinTheta*sin(phi)
	
    if(.not.present(direction)) then 
        ugamma=(ugamma+DriftBeta*gamma)*DriftGamma !Lorentz boost 		
    else		
         dirn=direction  ! store the values in local memory 
         dirn_mag=sqrt(dirn(1)**2+dirn(2)**2+dirn(3)**2)
   
         if(dirn_mag.eq.0) then 
              dirn=0
         else 
              dirn=dirn/dirn_mag
         end if
               
         projn=ugamma*dirn(1)+vgamma*dirn(2)+wgamma*dirn(3)
         
         ugamma=ugamma+(DriftGamma-1)*projn*dirn(1) + DriftGamma*DriftBeta*gamma*dirn(1)
         vgamma=vgamma+(DriftGamma-1)*projn*dirn(2) + DriftGamma*DriftBeta*gamma*dirn(2)
         wgamma=wgamma+(DriftGamma-1)*projn*dirn(3) + DriftGamma*DriftBeta*gamma*dirn(3)
    end if
end subroutine AddDriftToThermalVel

subroutine InitMaxwellJuttnerTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
     integer :: GammaTableSize
     real(psn), intent(inout), dimension(GammaTableSize) :: Gamma_Table,PDF_Table
     real(psn), intent(in) :: Temp
     real(psn) ::TempThis
     real(psn) :: GammaMax,GammaMin,dGamma
     real(psn) :: PDF_sum
     integer :: i
     
     !Create a table of Log(Gamma) in a specified range, linearly spaced in log scale 
     TempThis=max(Temp,1.0e-8)
     GammaMax=10.0_psn*TempThis+1.0_psn
     GammaMin=1.0_psn
     dGamma=max(((GammaMax-GammaMin)/(GammaTableSize-1)),0.00001_psn)

     do i=1,GammaTableSize
          Gamma_Table(i)=GammaMin+dGamma*(i-1)
     end do 
     
     PDF_sum=0
     PDF_Table(1)=0
     do i=2,GammaTableSize
#ifdef twoD 
       PDF_sum=PDF_sum+Gamma_Table(i)*exp(-(Gamma_Table(i)-1.0_psn)/TempThis)
#else           
        PDF_sum=PDF_sum+Gamma_Table(i)*sqrt(Gamma_Table(i)**2-1.0_psn)*exp(-(Gamma_Table(i)-1)/TempThis) 
#endif
       PDF_Table(i)=PDF_sum
     end do 
     PDF_Table=PDF_Table/PDF_sum !Normalisation     
end subroutine InitMaxwellJuttnerTable

subroutine InitMaxwellNonRelTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
	implicit none 
	 integer :: GammaTableSize
     real(psn), intent(inout), dimension(GammaTableSize) :: Gamma_Table,PDF_Table
     real(psn), intent(in) :: Temp
     real(psn) ::TempThis
     real(psn) :: VelMax,VelMin,dVel
     real(psn) :: PDF_sum
     integer :: i
     
     !Create a table of Log(Gamma) in a specified range, linearly spaced in log scale 
     TempThis=Temp!max(Temp,1.0e-15)
     VelMax=min(1.0,5.0_psn*sqrt(TempThis))
     VelMin=0.0_psn
     dVel=(VelMax-VelMin)/(GammaTableSize-1)
 
     do i=1,GammaTableSize
          Gamma_Table(i)=VelMin+dVel*(i-1)
     end do  
	    
     PDF_sum=0
     PDF_Table(1)=0
     do i=2,GammaTableSize
       PDF_sum=PDF_sum+(Gamma_Table(i)**2)*exp(-(Gamma_Table(i)**2)/(2*TempThis))
	   !if(proc.eq.0) print*,'GammaTable',Gamma_Table(i),'PDF_sum',PDF_sum,'Temp',TempThis,'exp',exp(-(Gamma_Table(i)**2)/(2*TempThis)),'dv',dVel
       PDF_Table(i)=PDF_sum
     end do
	 if((PDF_sum.eq.0).or.(PDF_sum.ne.PDF_sum)) then !to take care of zero temprature case 
		 PDF_Table=1.0_psn
		 PDF_Table(1)=0.0_psn
		 PDF_sum=1.0_psn
	 end if 
     PDF_Table=PDF_Table/PDF_sum !Normalisation
end subroutine InitMaxwellNonRelTable

subroutine Init_PUI_ShellDist_Table(GammaTableSize,Gamma_Table,PDF_Table,WindSpeed)
	 implicit none 
	 integer :: GammaTableSize
     real(psn), intent(inout), dimension(GammaTableSize) :: Gamma_Table,PDF_Table
     real(psn), intent(in) :: WindSpeed ! Solar wind speed, same as the upstream drift speed 
     real(psn) :: VelMax,VelMin,dVel,Vratio
     real(psn) :: PDF_sum
     integer :: i
     
     !Create a table of Log(Gamma) in a specified range, linearly spaced in log scale 
     VelMax=min(1.0_psn,5.0_psn*WindSpeed)
     VelMin=0.0_psn
     dVel=(VelMax-VelMin)/(GammaTableSize-1)

     do i=1,GammaTableSize
          Gamma_Table(i)=VelMin+dVel*(i-1)
     end do  
	    
     PDF_sum=0
     PDF_Table(1)=0
     do i=2,GammaTableSize
	   Vratio= (Gamma_Table(i)/WindSpeed)**(2.0_psn -(3.0_psn/2.0_psn))
       if(Gamma_Table(i).gt.WindSpeed) Vratio=0.0_psn
       if((Gamma_Table(i-1).le.WindSpeed).and.(Gamma_Table(i).ge.WindSpeed)) then
		   Vratio=1.0_psn
	   else
		   Vratio=0.0_psn
	   end if
	   
	   !PDF_sum=PDF_sum+Vratio*exp(-(6.0_psn/90.0_psn)*Vratio)
       PDF_sum=PDF_sum+Vratio!*exp(-(6.0_psn/90.0_psn)*Vratio)
       PDF_Table(i)=PDF_sum
     end do
     PDF_Table=PDF_Table/PDF_sum !Normalisation
end subroutine Init_PUI_ShellDist_Table

subroutine Init_KappaDist_Table(GammaTableSize,Gamma_Table,PDF_Table,Temp,KappaIndex)
	implicit none 
	 integer :: GammaTableSize
     real(psn), intent(inout), dimension(GammaTableSize) :: Gamma_Table,PDF_Table
     real(psn), intent(in) :: Temp ! Solar wind speed, same as the upstream drift speed 
	 integer               :: KappaIndex
     real(psn) :: VelMax,VelMin,dVel,Vratio
     real(psn) :: PDF_sum
     integer :: i
     
     !Create a table of Log(Gamma) in a specified range, linearly spaced in log scale 
     VelMax=min(1.0_psn,20.0_psn*sqrt(Temp))
     VelMin=0.0_psn
     dVel=(VelMax-VelMin)/(GammaTableSize-1)

     do i=1,GammaTableSize
          Gamma_Table(i)=VelMin+dVel*(i-1)
     end do  
	    
     PDF_sum=0
     PDF_Table(1)=0
     do i=2,GammaTableSize
	   Vratio= (1+(Gamma_Table(i)**2)/(2*KappaIndex*Temp))**(-KappaIndex-1)
	   !PDF_sum=PDF_sum+Vratio*exp(-(6.0_psn/90.0_psn)*Vratio)
       PDF_sum=PDF_sum+Vratio
	   PDF_Table(i)=PDF_sum
     end do
     PDF_Table=PDF_Table/PDF_sum !Normalisation
end subroutine Init_KappaDist_Table

subroutine BinarySearch(N,Arr,val,index)
     integer :: N
     real(psn),intent(in), dimension(N):: Arr
     real(psn),intent(in) :: val
     integer :: index,imin,imax,imid,imid_prev
     logical :: search
     
     imin=1
     imax=N
	 imid=(imin+imax)/2
     search=.true. 
	 do while(search)
		  imid_prev=imid 
          if(Arr(imid).ge.val) imax=imid
          if(Arr(imid).le.val) imin=imid
		  imid=(imin+imax)/2
          if(imin+1.ge.imax) search=.false.
		  if(imid.eq.imid_prev) search=.false. !to ensure exit from the loop if Arr(imin)=Arr(imax) and imin+1<imax 		  
     end do 
     index=imin	 
end subroutine BinarySearch 


!--------------------------------------------------------------------------------------------
!                    Poisson distribution 
!--------------------------------------------------------------------------------------------
subroutine GetIntPoissonDist(TableSize,Table,num)
	integer :: TableSize
	real(psn), dimension(TableSize) ::Table
	integer    :: num
	real(psn)  :: r1
	integer    :: index
    call random_number(r1)   
    call BinarySearch(TableSize,Table,r1,index)     
	num=index-1
	!print*,'r1',r1,'index',index,'Table',Table(index),Table(index+1)
end subroutine GetIntPoissonDist
subroutine InitPoissonDist(TableSize,Table,Mean)
	integer :: TableSize
	real(psn), dimension(TableSize) ::Table
	real(psn) :: Mean,Table_sum
	real(dbpsn) :: fact
	integer   :: i
	
	Table(1)=0.0_psn
	Table_sum=Table(1)
	
	if(Mean.eq.0) then ! Mean=0 is a special case 
		Table=1.0
		Table(1)=0.0
		return
	end if 
	
	do i=2,TableSize
		call lnFacN(i-2,fact)
		Table_sum = Table_sum + exp(-Mean+(i-2)*log(Mean)-fact)!note that index=i+1 holds probablity of number i  
	    Table(i) = Table_sum
	end do 
	!Table=Table/Table_sum !Normalization
! 	do i=1,100
! 		print*,'i',i,Table(i),Mean
! 	end do
end subroutine InitPoissonDist

subroutine lnFacN(num,res) ! see the wikipedia page for Stirling series
	integer :: num,i
	real(dbpsn) :: res,n,n1
	if(num.lt.5) then 
		n=1
		do i=1,num
			n=n*i
		end do  
		res=log(n)
	else 
		n=num!store the value in double precision variable 
		res=n*log(n)-n+0.5d0*log(2.0d0*3.14159265358979323846264d0*n)
		n1=1.0d0/n
		res=res+(1.0d0/12.0d0)*n1 ! first term 
		n1=n1/(n*n)
		res=res-(1.0d0/360.0d0)*n1 ! second term
		n1=n1/(n*n)
		res=res+(1.0d0/1260.0d0)*n1 ! Third term
		n1=n1/(n*n)
		res=res-(1.0d0/1680.0d0)*n1 ! Fourth term
    end if 
end subroutine lnFacN

!-------End of Poisson distribution subroutines --------------------------------------------

!--------------------------------------------------------------------------------------------
!
!        Some Auxiliary Subroutines and Functions 
!
!--------------------------------------------------------------------------------------------
 
!The Following subroutine returns a unit vector perpendicular to input vectors A and B
subroutine GetUnitVectorAlongAXB(Ax,Ay,Az,Bx,By,Bz,Ox,Oy,Oz)
     real(psn), intent(IN)    :: Ax,Ay,Az
     real(psn), intent(IN)    :: Bx,By,Bz
     real(psn), intent(INOUT) :: Ox,Oy,Oz
     real(psn)                :: Omag
     Ox=Ay*Bz-Az*By
     Oy=Az*Bx-Ax*Bz
     Oz=Ax*By-Ay*Bx
     Omag=sqrt(Ox*Ox+Oy*Oy+Oz*Oz)
     if(Omag.eq.0) then 
          Ox=0.0_psn
          Oy=0.0_psn
          Oz=0.0_psn
     else 
          Ox=Ox/Omag
          Oy=Oy/Omag
          Oz=Oz/Omag
     end if
end subroutine GetUnitVectorAlongAXB 

!The following subroutine returns global cordinate (in double precision) of a point on this proc. 
subroutine GetGlobalPosition_DP(xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
     real        :: xlocal,ylocal,zlocal
     real(dbpsn) :: xglobal,yglobal,zglobal
     
    xglobal=xlocal-3+xborders(procxind(proc))
    yglobal=ylocal-3+yborders(procyind(proc))
#ifdef twoD     
     zglobal=zlocal
#else
    zglobal=zlocal-3+zborders(proczind(proc)) 
#endif 

end subroutine GetGlobalPosition_DP


!--------------------------------------------------------------------------------------------
!
!        Some Auxiliary Subroutines useful for Shock problem 
!
!--------------------------------------------------------------------------------------------

!the following subroutines tries to find the current location of the shock by following total magnetic energy vs. x 
subroutine FindShockXpos_BmagPeak(Xshock)
	implicit none 
	integer :: Xshock
    integer :: i,j,k,i1,j1,k1
	integer :: xpeak
	real(psn) :: mag_peak
    real(psn), dimension(mx) :: mag1D_slice,mag1D_this_proc
	real(psn), dimension(0:nproc-1) :: mag_peak_all_proc
	integer, dimension(0:nproc-1)   :: xpeak_all_proc
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat 
	integer :: mpi_err
	!first compute the local peak and its location 	
	
	mag1D_this_proc=0
    do i=3,mx-3    
#ifndef twoD 
          do k=3,mz-3
#else
          do k=1,1
#endif
                do j=3,my-3
                    mag1D_this_proc(i)=mag1D_this_proc(i)+(Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)
                end do
          end do
	 end do	
	 
	 !now communicate with all relevant proc to get the global sum in Y-Z plane
	 i1=procxind(proc)
	 j1=procyind(proc)
	 k1=proczind(proc)
	 mag1D_slice=0.0_psn
	 if(j1.eq.0.and.k1.eq.0) then !first reduce all local 1D array one one proc.
		 mag1D_slice=mag1D_this_proc
#ifdef twoD
         do k=0,0
#else          		 
		 do k=0,nSubDomainsZ-1
#endif
			 do j=0,nSubDomainsY-1
				 if(j.eq.0.and.k.eq.0) cycle
				     call MPI_RECV(mag1D_this_proc,mx,mpi_psn,proc_grid(i1,j,k),0,MPI_COMM_WORLD,mpi_stat,mpi_err)
					 mag1D_slice=mag1D_slice+mag1D_this_proc ! the local array is used as recv buffer
			 end do
		 end do
	 else
		 call MPI_SEND(mag1D_this_proc,mx,mpi_psn,proc_grid(i1,0,0),0,MPI_COMM_WORLD,mpi_err)
	 end if
	 
	 
	 if(j1.eq.0.and.k1.eq.0) then !now send the reduced array to all proc. working on this YZ slice
#ifdef twoD 		 
		 do k=0,0
#else
         do k=0,nSubDomainsZ-1
#endif 
			 do j=0,nSubDomainsY-1
				 if(j.eq.0.and.k.eq.0) cycle
				     call MPI_SEND(mag1D_slice,mx,mpi_psn,proc_grid(i1,j,k),0,MPI_COMM_WORLD,mpi_err)
			 end do
		 end do
	 else
		 call MPI_RECV(mag1D_slice,mx,mpi_psn,proc_grid(i1,0,0),0,MPI_COMM_WORLD,mpi_stat,mpi_err)
	 end if


	 !Compute the local peak and the index corresponding to the peak
     mag_peak=0.0_psn
	 xpeak=0
	 do i=3,mx-3
	     if(mag1D_slice(i).gt.mag_peak) then
		     mag_peak=mag1D_slice(i)
		     xpeak=i+xborders(procxind(proc))-3
	      end if
	 end do

	 !Now communicate the local location of the magentic peak and find the global peak
     call MPI_ALLGATHER(mag_peak,1,mpi_psn,mag_peak_all_proc,1,mpi_psn,MPI_COMM_WORLD,mpi_err)
     call MPI_ALLGATHER(xpeak,1,MPI_INTEGER,xpeak_all_proc,1,MPI_INTEGER,MPI_COMM_WORLD,mpi_err)

	 mag_peak=0.0_psn
	 do i=0,nproc-1 !compute the global peak
		 if(mag_peak_all_proc(i).gt.mag_peak) then
			 mag_peak=mag_peak_all_proc(i)
			 Xshock=xpeak_all_proc(i)
		 end if
	 end do
end subroutine FindShockXPos_BmagPeak 


end module help_setup