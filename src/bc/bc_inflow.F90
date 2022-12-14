module bc_inflow 
	use parameters
	use vars
	use deposit
	use mem_prtl
	use init_prtl
#ifdef gpu
    use fields_gpu
	use particles_gpu
#endif	
	implicit none 
	
	integer :: nInjPrtl=0
	
	type InjPrtlType
		integer :: side=0 ! 1: x,left; 2: x,right; 3: y,bottom; 4: y,top  
		procedure(scalar_global), pointer, nopass :: Den =>null()
		procedure(vector_global), pointer, nopass :: Drift =>null()		
		integer, dimension(:), allocatable :: pid
	    character (len=4) :: prtl_placement	 
		real(psn) :: max_displacement ! max. dispalcement is used to add some noise in psitioning the particles
	end type InjPrtlType
	
	type(InjPrtlType), dimension(:), allocatable :: InjPrtl
	
contains 
!------------------------------------------	
! Particles
!------------------------------------------	
	
    subroutine InflowBC_Prtl
		integer :: n 
		real(psn) :: pos
        
		if(bc_face(2)%type_prtl.eq.'iflw') then 
			do n=1,nInjPrtl
				if(InjPrtl(n)%prtl_placement.eq.'rndm') then 
					call RefPrtlFluidFrameRight(bc_face(2)%flw)
					exit
				end if
			end do
		end if

				
		do n=1,nInjPrtl	
			if(InjPrtl(n)%prtl_placement.eq.'rndm') call InjectNewPrtl(InjPrtl(n))
			if(InjPrtl(n)%prtl_placement.eq.'even') call InjectNewPrtl_UniformSeparation(InjPrtl(n))
		end do
#ifdef gpu
        call SendPrtlToGPU(Nchunk_prtl_gpu)
#endif

	end subroutine InflowBC_Prtl



!---------------------------------------------------------------------------------------------------
!
!       Reflect the backward-moving particles in the fluid frame
!
!---------------------------------------------------------------------------------------------------

	subroutine RefPrtlFluidFrameRight(flw)
		type(FlowFldType) :: flw
		real(psn) :: x1,vx,vy,vz,x0,y0,z0
		real(dbpsn) :: ox,oy,oz
		real(psn) :: xinj_local
		integer :: n, i1

		xinj_local = bc_face(2)%pos_prtl - xborders(procxind)+3
		if( xinj_local .gt. mx-1) return

#ifdef gpu 
		call LoadGPUOutlierRight( real(xinj_local,psn) - flw%vmax * c) ! bring particles near the right edge from GPU
        call RecvFullDomainCurrFromGPU
#endif		

		call GetOxyz(ox,oy,oz)
		

		do n=1,used_prtl_arr_size
			call flw%Drift(xp(n)+ox,yp(n)+oy,zp(n)+oz,vx,vy,vz) 
			x1= xinj_local + vx*c ! instantaneous location of the drifted boundary layer
			if((qp(n).ne.0).and.(xp(n).gt.x1)) then
				call ReflectPrtlMovingFrameX(vx,up(n),vp(n),wp(n))
				x0=xp(n)
				y0=yp(n)
				z0=zp(n)
				xp(n)=2.0_psn*x1-x0
				call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
			end if
		end do
		
#ifdef gpu
        call SendFullDomainCurrToGPU
#endif		
	end subroutine RefPrtlFluidFrameRight
	
	
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
	
	
	subroutine GetOxyz(ox,oy,oz)
		real(dbpsn) ::  ox,oy,oz ! x,y, and z poistion (global) of the lower left corner of this subdomain
		ox=xborders(procxind)-3.0_psn
		oy=yborders(procyind)-3.0_psn
		oz=zborders(proczind)-3.0_psn
	end subroutine GetOxyz

	!------------------------------------------------------------
	!
	!    Inject new particles in the emptied out domain, random position
	!
	!------------------------------------------------------------
	
	subroutine InjectNewPrtl(inj)
		type(InjPrtlType) :: inj
		real(psn) :: x1,x2,y1,y2,z1,z2,q1,q2
		real(dbpsn) :: xglobal,yglobal,zglobal
		integer :: nprtl_new,n,tag, side
		real(psn) :: rnd_acpt,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
		real(psn) :: vx,vy,vz
		real(psn) :: Temp
		logical   :: acpt

		side = inj%side
		x1 = xmin; x2 = xmax; y1 = ymin; y2 = ymax; z1 = zmin; z2 = zmax 
		call SetInjDomain(inj,nprtl_new,x1,x2,y1,y2,z1,z2)
		
		if(used_prtl_arr_size+size(inj%pid)*nprtl_new.gt.prtl_arr_size) call ReshapePrtlArr(used_prtl_arr_size+size(inj%pid)*nprtl_new+1000000)

		do n=1,nprtl_new
			 
		    call GetNewPrtlPos(x1,x2,y1,y2,z1,z2,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
			call inj%Drift(xglobal,yglobal,zglobal,vx,vy,vz)
			
			acpt = .false. 
			
			call random_number(rnd_acpt)
			if(rnd_acpt.le.inj%Den(xglobal,yglobal,zglobal)) acpt = .true. ! comply with the local density
			
			!fill only the emptied out domain
			select case (side)
				case(1)
					if(xglobal.gt.bc_face(1)%pos_prtl+vx*c) acpt = .false.
				case(2)
					if(xglobal.lt.bc_face(2)%pos_prtl+vx*c) acpt = .false.
				case(3)
					if(yglobal.gt.bc_face(3)%pos_prtl+vy*c) acpt = .false.
				case(4)
					if(yglobal.lt.bc_face(4)%pos_prtl+vy*c) acpt = .false.
			end select
			
			
			if(acpt) then	

				call InsertPrtl_PSP(inj%pid,xglobal,yglobal,zglobal, xlocal,ylocal,zlocal, 1.0_psn ) 
								
			end if 
		end do 
	end subroutine InjectNewPrtl
	
	subroutine SetInjDomain(inj,nprtl_new,x1,x2,y1,y2,z1,z2)
		type(InjPrtlType) :: inj
		real(psn) :: x1,x2,y1,y2,z1,z2
		real(dbpsn) :: pos
		real(psn) :: vmax, shift
		integer   :: nprtl_new
		integer   :: side
		
		side = inj%side
		pos = bc_face(side)%pos_prtl
		vmax = bc_face(side)%flw%vmax
		shift = bc_face(side)%speed*c
		
	   	if(side.eq.1) then 
		     x1=max(pos-shift,real(xborders(procxind),dbpsn)) -xborders(procxind)+3
	   	     x2=min(pos+vmax*c,real(xborders(procxind+1),dbpsn)) -xborders(procxind)+3
		end if
		if(side.eq.2) then 
		     x1=max(pos-vmax*c,real(xborders(procxind),dbpsn)) -xborders(procxind)+3
	   	     x2=min(pos+shift,real(xborders(procxind+1),dbpsn)) -xborders(procxind)+3	  
		end if
	   	if(side.eq.3) then 
		     y1=max(pos-shift,real(yborders(procyind),dbpsn)) -yborders(procyind)+3
	   	     y2=min(pos+vmax*c,real(yborders(procyind+1),dbpsn)) -yborders(procyind)+3
		end if
		if(side.eq.4) then 
		     y1=max(pos-vmax*c,real(yborders(procyind),dbpsn)) -yborders(procyind)+3
	   	     y2=min(pos+shift,real(yborders(procyind+1),dbpsn)) -yborders(procyind)+3	 
		end if
	   	if(side.eq.5) then 
		     z1=max(pos-shift,real(zborders(proczind),dbpsn)) -zborders(proczind)+3
	   	     z2=min(pos+vmax*c,real(zborders(proczind+1),dbpsn)) -zborders(proczind)+3
		end if
		if(side.eq.6) then 
		     z1=max(pos-vmax*c,real(zborders(proczind),dbpsn)) -zborders(proczind)+3
	   	     z2=min(pos+shift,real(zborders(proczind+1),dbpsn)) -zborders(proczind)+3	 
		end if

		nprtl_new= epc*max(0.0_psn,x2-x1)*max(0.0_psn,y2-y1)*max(0.0_psn,z2-z1)
		if(nprtl_new.eq.0) return 
		if(nprtl_new.lt.1000) call GetIntPoissonDist(real(real(nprtl_new),psn),nprtl_new)
		
	end subroutine SetInjDomain
	
	subroutine GetNewPrtlPos(x1,x2,y1,y2,z1,z2,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
		real(psn)   :: x1,x2,y1,y2,z1,z2
		real(psn)   :: xlocal,ylocal,zlocal
		real(dbpsn) :: xglobal,yglobal,zglobal
		real(psn)   :: r1,r2,r3
	 		call random_number(r1)
	 	    call random_number(r2)
	 		call random_number(r3)
	 		xlocal=x1+r1*(x2-x1)
			ylocal=y1+r2*(y2-y1)
			zlocal=z1+r3*(z2-z1)
			
			xglobal=xlocal+(xborders(procxind)-3)
			yglobal=ylocal+(yborders(procyind)-3)
#ifndef twoD			
			zglobal=zlocal+(zborders(proczind)-3)
#else
            zglobal=zlocal
#endif 	

	end subroutine GetNewPrtlPos
	
	
	!------------------------------------------------------------
	!
	!    Inject new particles in the empty domain, uniform grid-like placement
	!    Note: this scheme is currenly not designed to work with a drifting flow
	!    ideal : if vmax is not zero and also pass drift function to accept and reject within  "InitPrtlFromPSP_UniformSeparation"
	!------------------------------------------------------------
	subroutine InjectNewPrtl_UniformSeparation(inj)
		type(InjPrtlType) :: inj
		real(dbpsn) :: x1,x2,y1,y2,z1,z2
		real(dbpsn) :: pos
		real(psn) :: vmax, shift
		integer :: side
		
		side = inj%side
		x1 = xborders(procxind); x2 = xborders(procxind+1); y1 = yborders(procyind); y2 = yborders(procyind+1); z1 = zborders(proczind); z2 = zborders(proczind+1) 
		pos = bc_face(side)%pos_prtl
		vmax = 0!bc_face(side)%flw%vmax 
		shift = bc_face(side)%speed*c
		
	   	if(side.eq.1) then 
		     x1=max(pos-shift,real(xborders(procxind),dbpsn)) 
	   	     x2=min(pos+vmax*c,real(xborders(procxind+1),dbpsn))
		end if
		if(side.eq.2) then 
		     x1=max(pos-vmax*c,real(xborders(procxind),dbpsn)) 
	   	     x2=min(pos+shift,real(xborders(procxind+1),dbpsn)) 
		end if
	   	if(side.eq.3) then 
		     y1=max(pos-shift,real(yborders(procyind),dbpsn)) 
	   	     y2=min(pos+vmax*c,real(yborders(procyind+1),dbpsn))
		end if
		if(side.eq.4) then 
		     y1=max(pos-vmax*c,real(yborders(procyind),dbpsn)) 
	   	     y2=min(pos+shift,real(yborders(procyind+1),dbpsn))  
		end if
	   	if(side.eq.5) then 
		     z1=max(pos-shift,real(zborders(proczind),dbpsn)) 
	   	     z2=min(pos+vmax*c,real(zborders(proczind+1),dbpsn))
		end if
		if(side.eq.6) then 
		     z1=max(pos-vmax*c,real(zborders(proczind),dbpsn)) 
	   	     z2=min(pos+shift,real(zborders(proczind+1),dbpsn))  
		end if
		

		call InitPrtlFromPSP_UniformSeparation(inj%pid, inj%Den, x1,x2,y1,y2,z1,z2, inj%max_displacement)
		
		
	end subroutine InjectNewPrtl_UniformSeparation

	
end module bc_inflow
