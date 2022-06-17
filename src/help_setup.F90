module help_setup
     use parameters
     use vars
	 use communication
     use memory
	 use prob_dist
	 use prtl_stats 
	 use savedata
	 use fields
	 use prtl_tag
     implicit none 

contains 
!---------------------------------------------------------------------------------
!  subroutines to initialise particles 
!---------------------------------------------------------------------------------		
   
	subroutine InitPrtl(Flvr1,Flvr2,Density,Temperature1,Temperature2,DriftVelocity1,DriftVelocity2,SpeedDist1,SpeedDist2,Vmax1,Vmax2)
		integer, optional :: Flvr1
		integer, optional :: Flvr2
		procedure(scalar_global) :: Density
		procedure(scalar_global), optional :: Temperature1
		procedure(scalar_global), optional :: Temperature2
		procedure(func1D),        optional :: SpeedDist1
		procedure(func1D),        optional :: SpeedDist2
		real(psn), optional :: Vmax1 !the maximum particle speed in the plasma frame for Flvr1, default is c 
		real(psn), optional :: Vmax2 ! max. speed for Flvr2; Vmax1 and Vmax2 are only used when speed dist is provided
		procedure(vector_global), optional :: DriftVelocity1  
		procedure(vector_global), optional :: DriftVelocity2  
		!!----------- End of input variables ------------!!		
		integer :: i1, i2
		integer, dimension(:) , allocatable :: pid
		
		call DefinePhaseSpaceProperty(Flvr1,Density,Temperature1,SpeedDist1,Vmax1,DriftVelocity1,i1)
		allocate(pid(1))
		pid = (/ i1 /)
		
		if(present(Flvr2)) then
			call DefinePhaseSpaceProperty(Flvr2,Density,Temperature2,SpeedDist2,Vmax2,DriftVelocity1,i2)
			deallocate(pid)
			allocate(pid(2))
			pid = (/ i1, i2 /)
		end if 
		
		call InitPrtlFromPSP(pid, Density)
			
	end subroutine InitPrtl
	
	subroutine InitPrtlFromPSP(pid, Density)
		integer, dimension(:) :: pid
		procedure(scalar_global) :: Density
		
		real(dbpsn) :: r1,r2,r3,rnd_acpt
		real(dbpsn) :: xglobal, yglobal, zglobal
		real(psn) :: xlocal, ylocal, zlocal
		integer   :: n, est_np
		
		!make sure that the particle array is large enough 
		est_np = size(pid)*EstimatePrtlCount(Density, Nelc)
		if(used_prtl_arr_size + est_np .gt. prtl_arr_size) call ReshapePrtlArr( new_size=int(prtl_arr_size + est_np) )
		
		do n=1,Nelc
			call random_number(r1)
		    call random_number(r2)
			call random_number(r3)
			call random_number(rnd_acpt)
		
			xglobal= xborders(procxind(proc)) + r1*(xborders(procxind(proc)+1)-xborders(procxind(proc)))
			yglobal= yborders(procyind(proc)) + r2*(yborders(procyind(proc)+1)-yborders(procyind(proc)))
			zglobal= zborders(proczind(proc)) + r3*(zborders(proczind(proc)+1)-zborders(proczind(proc)))
		
			if(rnd_acpt.le.Density(xglobal,yglobal,zglobal)) then
				
			    xlocal= xglobal - xborders(procxind(proc)) + xmin
			    ylocal= yglobal - yborders(procyind(proc)) + ymin
			    zlocal= zglobal - zborders(proczind(proc)) + zmin
				call InsertPrtl_PSP(pid,xglobal,yglobal,zglobal, xlocal,ylocal,zlocal ) 
		 
			end if 
		end do 
		
		
	end subroutine InitPrtlFromPSP
	
!---------------------------------------------------------------------------------
!  Define a new Phase Space Property
!  pid (INTEGER) = handle/pointer to this definition
!---------------------------------------------------------------------------------			
	subroutine DefinePhaseSpaceProperty(Flvr,Density,Temperature,SpeedDist,Vmax,DriftVelocity,pid)
		integer :: Flvr
		procedure(scalar_global), optional :: Density
		procedure(scalar_global), optional :: Temperature
		procedure(func1D),        optional :: SpeedDist
		procedure(vector_global), optional :: DriftVelocity  
		real(psn), optional :: Vmax !the maximum particle speed in the plasma frame, default is c 
		integer :: pid
		
		used_PSP_list_size = used_PSP_list_size +1
		pid = used_PSP_list_size
		
		call SetPhaseSpaceProperty(PSP_list(pid),Flvr,Density,Temperature,SpeedDist,Vmax,DriftVelocity)
		
	end subroutine DefinePhaseSpaceProperty 
!---------------------------------------------------------------------------------
!  Set Load Balancing Type. The loadbalancing scheme is applied (loadbalance.F90) based on the type
!---------------------------------------------------------------------------------	
    subroutine SetLoadBalancing(Type)
		character (len=*) :: Type
		if(Type.eq.'shock') load_balancing_type = 1 ! 1 = shock
	end subroutine SetLoadBalancing

!---------------------------------------------------------------------------------
!  Set an external current source ::  pass a subroutine defined in the setup file
!---------------------------------------------------------------------------------	
	subroutine SetExternalCurrent(J)
		procedure(vector_global) :: J

		ext_current_present = .true.
		J_ext => J
	end subroutine SetExternalCurrent		
	
!---------------------------------------------------------------------------------
!  subroutines to initialise Electric and Magnetic fields
!---------------------------------------------------------------------------------		
    subroutine InitElectricField(Efld)
		external :: Efld
		integer :: i,j,k
		real(dbpsn) :: x,y,z
		real(psn) :: e_x,e_y,e_z

#ifdef twoD
        do k=1,1
#else 	
	    do k=1,mz
#endif 
        do j=1,my 
		    do i=1,mx
				x= i -3.0_dbpsn + xborders(procxind(proc))
				y= j -3.0_dbpsn + yborders(procyind(proc))
#ifdef twoD    
                z=0.0_dbpsn
#else				 				
				z= k - 3.0_dbpsn + zborders(proczind(proc))
#endif 				
			    call Efld(x + 0.5_dbpsn,y,z,e_x,e_y,e_z)
				Ex(i,j,k)=e_x
			    call Efld(x,y + 0.5_dbpsn,z,e_x,e_y,e_z)
				Ey(i,j,k)=e_y
#ifdef twoD				
				call Efld(x,y,0.0_dbpsn,e_x,e_y,e_z)
#else				
				call Efld(x,y,z + 0.5_dbpsn,e_x,e_y,e_z)
#endif				
			
				Ez(i,j,k)=e_z
		    end do 
        end do  
        end do 
	end subroutine InitElectricField
	
    subroutine InitMotionalElectricField(Bfld,Vfld)
		external :: Bfld, Vfld
		integer :: i,j,k
		real(dbpsn) :: x,y,z
		real(psn) :: b_x,b_y,b_z,vx,vy,vz

#ifdef twoD
        do k=1,1
#else 	
	    do k=1,mz
#endif 
        do j=1,my 
		    do i=1,mx
				x= i -3.0_dbpsn + xborders(procxind(proc))
				y= j -3.0_dbpsn + yborders(procyind(proc))
#ifdef twoD    
                z=0.0_dbpsn
#else				 				
				z= k -3.0_dbpsn + zborders(proczind(proc))
#endif 				
			    call Bfld(x + 0.5_dbpsn,y,z,b_x,b_y,b_z)
				call Vfld(x + 0.5_dbpsn,y,z,vx,vy,vz)
				Ex(i,j,k) = vz * b_y - vy * b_z 
				
			    call Bfld(x,y + 0.5_dbpsn,z,b_x,b_y,b_z)
				call Vfld(x,y + 0.5_dbpsn,z,vx,vy,vz)
				Ey(i,j,k) = vx * b_z - vz * b_x
				
#ifdef twoD				
				call Bfld(x,y,0.0_dbpsn,b_x,b_y,b_z)
				call Vfld(x,y,0.0_dbpsn,vx,vy,vz)
#else				
				call Bfld(x,y,z + 0.5_dbpsn,b_x,b_y,b_z)
				call Vfld(x,y,z + 0.5_dbpsn,vx,vy,vz)
#endif				
			
				Ez(i,j,k) = vy * b_x - vx * b_y
		    end do 
        end do  
        end do 
	end subroutine InitMotionalElectricField
	
    subroutine InitMagneticField(Bfld)
		external :: Bfld
		integer :: i,j,k
		real(dbpsn) :: x,y,z,zp_half
		real(psn) :: b_x,b_y,b_z

#ifdef twoD
        do k=1,1
#else 	
	    do k=1,mz
#endif 
        do j=1,my 
		    do i=1,mx
				x=i -3.0_dbpsn + xborders(procxind(proc))
				y=j -3.0_dbpsn + yborders(procyind(proc))
#ifdef twoD    
                z=0.0_dbpsn
				zp_half=0.0_dbpsn
#else				 				
				z= k - 3.0_dbpsn + zborders(proczind(proc))
				zp_half=z+0.5_dbpsn
#endif 				
				call Bfld(x,y+0.5_dbpsn,zp_half, b_x,b_y,b_z)
				Bx(i,j,k)=b_x
			    call Bfld(x+0.5_dbpsn,y,zp_half,b_x,b_y,b_z)
				By(i,j,k)=b_y				
				call Bfld(x+0.5_dbpsn,y+0.5_dbpsn,z,b_x,b_y,b_z)	
				Bz(i,j,k)=b_z
		    end do 
        end do  
        end do 
	end subroutine InitMagneticField
	
!------------------------------------------------------------
! Define a new particle species (or "Flavor")
!------------------------------------------------------------       
	 subroutine DefineNewFlvr(FlvID,Charge,QbyM,Type,SaveFld,SaveRatio)
		 integer :: FlvID
		 integer, optional :: Type,SaveFld,SaveRatio
		 real(psn) :: QbyM,Charge
		 
		 call SetQbyM(FlvID,QbyM)
		 FlvrCharge(FlvID)=Charge
		 if(present(Type)) then 
			 FlvrType(FlvID)=Type
		 else 
			 FlvrType(FlvID)=0
		 endif 
		 if(present(SaveFld)) then 
			 FlvrSaveFldData(FlvID)=SaveFld
		 else 
			 FlvrSaveFldData(FlvID)=0
		 end if 
		 if(present(SaveRatio)) then 
			 FlvrSaveRatio(FlvID)=SaveRatio
		 else 
			 FlvrSaveRatio(FlvID)=psave_ratio
		 end if 
		 CurrentTagID(FlvID)=0!1+proc*NtagProcLen 
		 call ReallocatePrtlTagArr
	 end subroutine DefineNewFlvr
	 
!------------------------------------------------------------
!Set Guide Field
!------------------------------------------------------------	
	subroutine SetBackgroundMagneticField(Bx,By,Bz)
		real(psn) :: Bx,By,Bz
		Bx_ext0=Bx
		By_ext0=By
		Bz_ext0=Bz
	end subroutine SetBackgroundMagneticField 


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


! !------------------------------------------------------------------------------------------------
! ! The following subroutine return postion of a new particle assuming that all particles are distributed uniformly in the simulation box
! !------------------------------------------------------------------------------------------------
!
!
! subroutine SmallAngleScattering(u,v,w,ThetaMean)
! 	real(psn) :: u,v,w,ThetaMean
! 	real(psn) :: theta,phi
! 	real(psn) :: u1,v1,w1,u2,v2,w2,u3,v3,w3,p1,p2,p3
! 	real(psn) :: mag
! 	real(psn) :: r1,r2
!
! 	call random_number(r2)
! 	theta=ThetaMean
! 	phi=2*pi*r2
!
! 	mag=sqrt(u**2+v**2+w**2)
!     p1=mag*cos(theta)
! 	p2=mag*sin(theta)*cos(phi) !New vector
! 	p3=mag*sin(theta)*sin(phi)
!
!
! 	mag=sqrt(u**2+v**2+w**2)
! 	if(mag.eq.0) return
! 	u1=u/mag
! 	v1=v/mag
! 	w1=w/mag
!
!     u2=-v1
!     v2=u1
!     w2=0.0_psn
! 	mag=sqrt(u2**2+v2**2+w2**2)
!     if(mag.ne.0) then
! 	    u2=u2/mag
! 	    v2=v2/mag
! 	else
! 		u2=1.0_psn
! 		v2=0.0_psn
! 	end if
!
!     u3=v2*w1-w2*v1
!     v3=w2*u1-u2*w1
!     w3=u2*v1-v2*u1
! 	mag=sqrt(u3**2+v3**2+w3**2)
! 	u3=u3/mag
! 	v3=v3/mag
! 	w3=w3/mag
!
! 	u= u1*p1 + u2*p2 + u3*p3
! 	v= v1*p1 + v2*p2 + v3*p3
! 	w= w1*p1 + w2*p2 + w3*p3
! end subroutine SmallAngleScattering
!
!

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





!========================================================================================================
!
!     Subroutines to generate mommetum distribution functions (outdated)
!
!========================================================================================================


!--------------------------------------------------------------------------------------------------------
!Distribution of Gamma of relativistic particles with drift
!a) Drift is along the x-axis (if < 1 it is beta, otherwise it is Lorentz factor Gamma)
!b) Here Temperature is in units of mc^2/k, i.e. Temp = kT/mc^2 , where T is the actual physical Temperature
!--------------------------------------------------------------------------------------------------------
!
! subroutine GetVelGamma_MaxwellJuttner(Drift,Temp,ugamma,vgamma,wgamma,direction)
!      real(psn), intent(in)  :: Drift,Temp
!      real(psn), dimension(3), intent(in), optional :: direction
!      real(psn), intent(out) :: ugamma,vgamma,wgamma !these are gamma*beta
!      real(psn)              :: Temp0=-1.0 ! values in the last call
!      integer, parameter :: GammaTableSize=10000
!      real(psn), save, dimension(GammaTableSize):: Gamma_Table,PDF_Table
!      real(psn) :: r1,gamma,beta
!      integer :: index
!
! 	 save Temp0
!
!      if(Temp0.ne.Temp) then
!           call InitMaxwellJuttnerTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
!           Temp0=Temp
!      end if
!
!      call random_number(r1)
!      call BinarySearch(GammaTableSize,PDF_Table,r1,index)
! 	 if(PDF_Table(index+1).eq.PDF_Table(index)) then
! 		 gamma=Gamma_Table(index)
! 	 else
!          gamma=Gamma_Table(index)+ ( (Gamma_Table(index+1)-Gamma_Table(index))/(PDF_Table(index+1)-PDF_Table(index)) )*(r1-PDF_Table(index))
!      end if
! 	 beta=sqrt((gamma-1)*(gamma+1))/gamma
!
! 	 call AddDriftToThermalVel(Drift,ugamma,vgamma,wgamma,gamma,beta,direction)
! end subroutine GetVelGamma_MaxwellJuttner
!
!
!
! !--------------------------------------------------------------------------------------------------------
! !Distribution of v/c of non-rel particles with drift
! ! GammaTable and Gamma_PDF_Table in this case should be understood to store velocity like in VelTable and Vel_PDF_Table
! !--------------------------------------------------------------------------------------------------------
!
!
!
!
!
! subroutine InitMaxwellJuttnerTable(GammaTableSize,Gamma_Table,PDF_Table,Temp)
!      integer :: GammaTableSize
!      real(psn), intent(inout), dimension(GammaTableSize) :: Gamma_Table,PDF_Table
!      real(psn), intent(in) :: Temp
!      real(psn) ::TempThis
!      real(psn) :: GammaMax,GammaMin,dGamma
!      real(psn) :: PDF_sum
!      integer :: i
!
!      !Create a table of Log(Gamma) in a specified range, linearly spaced in log scale
!      TempThis=max(Temp,0.00000001_psn)
!      GammaMax=10.0_psn*TempThis+1.0_psn
!      GammaMin=1.0_psn
!      dGamma=max(((GammaMax-GammaMin)/(GammaTableSize-1)),0.00001_psn)
!
!      do i=1,GammaTableSize
!           Gamma_Table(i)=GammaMin+dGamma*(i-1)
!      end do
!
!      PDF_sum=0
!      PDF_Table(1)=0
!      do i=2,GammaTableSize
! #ifdef twoD
!        PDF_sum=PDF_sum+Gamma_Table(i)*exp(-(Gamma_Table(i)-1.0_psn)/TempThis)
! #else
!         PDF_sum=PDF_sum+Gamma_Table(i)*sqrt(Gamma_Table(i)**2-1.0_psn)*exp(-(Gamma_Table(i)-1)/TempThis)
! #endif
!        PDF_Table(i)=PDF_sum
!      end do
!      PDF_Table=PDF_Table/PDF_sum !Normalisation
! end subroutine InitMaxwellJuttnerTable
!
!


! Older subroutines : Not in Use

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
!           if(FlvrSpare(flvp(n)).eq.0) cycle
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


end module help_setup