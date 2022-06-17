module injector 
	use parameters
	use vars
	use deposit
	use prob_dist
	use memory
	use initialise
	use communication
	use prtl_tag
#ifdef gpu
    use fields_gpu
	use initialise_gpu
	use particles_gpu
	use injector_gpu
#endif	
	implicit none 
	
	integer :: nInjPrtl=0
	
	type InflowFldType
		real(psn) :: Attenuate  
		procedure(vector_global), pointer, nopass  :: MagFld =>null()
		procedure(vector_global), pointer, nopass :: Drift =>null()
		real(psn) :: vmax	
	end type InflowFldType
	
	type InjPrtlType
		integer :: side=0 ! 1: x,left; 2: x,right; 3: y,bottom; 4: y,top  
		procedure(scalar_global), pointer, nopass :: Den =>null()
		procedure(vector_global), pointer, nopass :: Drift =>null()		
		integer, dimension(:), allocatable :: pid
	end type InjPrtlType
	
	type(InjPrtlType), dimension(:), allocatable :: InjPrtl
	type(InflowFldType), dimension(6) :: InflowFld ! 1: x,left; 2: x,right; 3: y,bottom; 4: y,top 
	
contains 
!------------------------------------------	
! Particles
!------------------------------------------	
	
    subroutine InflowBC_Prtl
		integer :: n 
		real(psn) :: pos
		!older method deleted older particles behind the injector and intialised a new population of particles 
! 		if(BC_Xmin_Prtl_Type.eq.'iflw') call ClearOldPrtlLeft(BC_Xmin_Prtl)
! 		if(BC_Xmax_Prtl_Type.eq.'iflw') call ClearOldPrtlRight(BC_Xmax_Prtl)
! 		if(BC_Ymin_Prtl_Type.eq.'iflw') call ClearOldPrtlBottom(BC_Ymin_Prtl)
! 		if(BC_Ymax_Prtl_Type.eq.'iflw') call ClearOldPrtlTop(BC_Ymax_Prtl)
        
		if(BC_Xmax_Prtl_Type.eq.'iflw') call RefPrtlFluidFrameRight(InflowFld(2))

				
		do n=1,nInjPrtl	
			call InjectNewPrtlPair(InjPrtl(n))
		end do
#ifdef gpu
        call SendPrtlToGPU(prtl_arr_size,xp,yp,zp,up,vp,wp,qp,tagp,flvp,var1p,buff_size_prtl_gpu,used_prtl_arr_size)
		qp(1:used_prtl_arr_size) = 0 
		used_prtl_arr_size = 0 
#endif
		
        call updateInflowBC_position_prtl ! update the injector position

	end subroutine InflowBC_Prtl

	
!     subroutine InflowBC_Prtl_Top(yinj)
! 		real(psn) :: yinj
! 		!call RefPrtlFluidFrameUpperY(yinj)
! 		call ClearOldPrtlTop(yinj)
! 		call InjectNewPrtlY(1,yinj)
! 	end subroutine InflowBC_Prtl_Top
!     subroutine InflowBC_Prtl_Bottom(yinj)
! 		real(psn) :: yinj
! 		!call RefPrtlFluidFrameLowerY(yinj)
! 		call ClearOldPrtlBottom(yinj)
! 		call InjectNewPrtlY(2,yinj)
! 	end subroutine InflowBC_Prtl_Bottom

!---------------------------------------------------------------------------------------------------
!
!       Reflect the backward-moving particles in the fluid frame
!
!---------------------------------------------------------------------------------------------------

	subroutine RefPrtlFluidFrameRight(flw)
		type(InflowFldType) :: flw
		real(psn) :: x1,vx,vy,vz,x0,y0,z0
		real(dbpsn) :: ox,oy,oz
		real(psn) :: xinj_local
		integer :: n, i1

		xinj_local = BC_Xmax_Prtl-xborders(procxind(proc))+3
		if( xinj_local .gt. mx-1) return

#ifdef gpu 
		call LoadGPUOutlierRight( real(xinj_local,psn) - InflowFld(2)%vmax * c) ! bring particles near the right edge from GPU
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
		ox=xborders(procxind(proc))-3.0_psn
		oy=yborders(procyind(proc))-3.0_psn
		oz=zborders(proczind(proc))-3.0_psn
	end subroutine GetOxyz

	!------------------------------------------------------------
	!
	!    Inject new particles in the emptied out domain
	!
	!------------------------------------------------------------
	
	subroutine InjectNewPrtlPair(inj)
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
					if(xglobal.gt.BC_Xmin_Prtl+vx*c) acpt = .false.
				case(2)
					if(xglobal.lt.BC_Xmax_Prtl+vx*c) acpt = .false.
				case(3)
					if(yglobal.gt.BC_Ymin_Prtl+vy*c) acpt = .false.
				case(4)
					if(yglobal.lt.BC_Ymax_Prtl+vy*c) acpt = .false.
			end select
			
			
			if(acpt) then	
				
				call InsertPrtl_PSP(inj%pid,xglobal,yglobal,zglobal, xlocal,ylocal,zlocal ) 
								
			end if 
		end do 
	end subroutine InjectNewPrtlPair
	
	subroutine SetInjDomain(inj,nprtl_new,x1,x2,y1,y2,z1,z2)
		type(InjPrtlType) :: inj
		real(psn) :: x1,x2,y1,y2,z1,z2
		real(dbpsn) :: pos
		real(psn) :: vmax, shift
		integer   :: nprtl_new
		integer   :: side
		
		side = inj%side
		pos = inflow_pos_prtl(side)
		vmax = InflowFld(side)%vmax
		shift = inflowBC_speed(side)*c
		
	   	if(side.eq.1) then 
		     x1=max(pos-shift,real(xborders(procxind(proc)),dbpsn)) -xborders(procxind(proc))+3
	   	     x2=min(pos+vmax*c,real(xborders(procxind(proc)+1),dbpsn)) -xborders(procxind(proc))+3
		end if
		if(side.eq.2) then 
		     x1=max(pos-vmax*c,real(xborders(procxind(proc)),dbpsn)) -xborders(procxind(proc))+3
	   	     x2=min(pos+shift,real(xborders(procxind(proc)+1),dbpsn)) -xborders(procxind(proc))+3	  
		end if
	   	if(side.eq.3) then 
		     y1=max(pos-shift,real(yborders(procyind(proc)),dbpsn)) -yborders(procyind(proc))+3
	   	     y2=min(pos+vmax*c,real(yborders(procyind(proc)+1),dbpsn)) -yborders(procyind(proc))+3
		end if
		if(side.eq.4) then 
		     y1=max(pos-vmax*c,real(yborders(procyind(proc)),dbpsn)) -yborders(procyind(proc))+3
	   	     y2=min(pos+shift,real(yborders(procyind(proc)+1),dbpsn)) -yborders(procyind(proc))+3	 
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
			
			xglobal=xlocal+(xborders(procxind(proc))-3)
			yglobal=ylocal+(yborders(procyind(proc))-3)
#ifndef twoD			
			zglobal=zlocal+(zborders(proczind(proc))-3)
#else
            zglobal=zlocal
#endif 	

	end subroutine GetNewPrtlPos
			
	!------------------------------------------	
	! EM Field Boundary conditions for the inflow boundaries
	!------------------------------------------	

	subroutine InflowBC_Fld(flw,side)
		type(InflowFldType) :: flw
		integer  :: side	

		!enlarge the domain size if needed
		if(inflowBC_speed(side).ne.0)  call EnlargeDomain(side)
        
		if(flw%Attenuate.gt.0) call AttenuateFld(flw,side) 
		 
		if(inflowBC_speed(side).ne.0) call updateInflowBC_position_fld(side) ! update the injector position, assuming no further actions are required
			
	end subroutine InflowBC_Fld
	
	!----------------------------------------------------------------------------------------------
	!  Update magentic field in the domain outside of boundaries 
	!----------------------------------------------------------------------------------------------
	subroutine InflowBC_SetMagFld(pos,side)
		type(InflowFldType) :: flw
		real(dbpsn):: pos
		integer  :: side
		real(psn):: b_x,b_y,b_z
		real(dbpsn) :: x0,y0,z0,xg,yg,zg
		integer :: i,j,k
		
		flw = InflowFld(side)
		
		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif	
	
		 do k=1,mz
	          do j=1,my
				  do i=1,mx
					  
					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=b_x
					  end if 

					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=b_y
					  end if 

					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=b_z
					  end if 
  
				  end do
			  end do  
	     end do
		
	end subroutine InflowBC_SetMagFld
	
	subroutine InflowBC_SetElcFld(pos,side)
		type(InflowFldType) :: flw
		real(dbpsn):: pos
		integer  :: side
		real(psn)::  b_x,b_y,b_z,vx,vy,vz
		real(dbpsn) :: xg,yg,zg, x0,y0,z0
		integer :: i,j,k
		
		flw = InflowFld(side)
		
		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif	
		
		 do k=1,mz
	          do j=1,my
				  do i=1,mx
		  
					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ex(i,j,k)= vz * b_y - vy * b_z 
					  end if

					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ey(i,j,k)= vx * b_z - vz * b_x
					  end if 

					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ez(i,j,k)= vy * b_x - vx * b_y
					  end if 
  
				  end do
			  end do  
	     end do
		
	end subroutine InflowBC_SetElcFld
	
	logical function exterior(xg,yg,zg,pos,side)
		real(dbpsn) :: xg,yg,zg,pos
		integer :: side
		exterior = .false.
		if(side.eq.1 .and. xg.le. pos) exterior = .true. 
        if(side.eq.2 .and. xg.ge. pos) exterior = .true. 
		if(side.eq.3 .and. yg.le. pos) exterior = .true. 
	    if(side.eq.4 .and. yg.ge. pos) exterior = .true. 
	end function exterior
	
	logical function outside_subdomain(x0,y0,z0,pos,side)
		real(dbpsn) :: pos, x0, y0, z0
		integer :: side
		outside_subdomain = .false.
		if(side.eq.1 .and. pos-x0 .lt.  1) outside_subdomain = .true. 
        if(side.eq.2 .and. pos-x0 .gt. mx) outside_subdomain = .true. 
		if(side.eq.3 .and. pos-y0 .lt.  1) outside_subdomain = .true. 
	    if(side.eq.4 .and. pos-y0 .gt. my) outside_subdomain = .true. 
    end function outside_subdomain
	
	
	!------------------------------------------------------------------------
	! Attenuation of the EM field componnets
	!------------------------------------------------------------------------
	subroutine AttenuateFld(flw,side)
		type(InflowFldType) :: flw
		integer  :: side
		real(psn):: b_x,b_y,b_z,vx,vy,vz, f0
		real(dbpsn) :: x0,y0,z0,xg,yg,zg
		real(dbpsn):: pos, dl
	    integer :: i,j,k
		
		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif	
        pos = inflow_pos_fld(side)
		dl = flw%Attenuate
		if( (side.eq.1 .or. side.eq.3) .and. outside_subdomain(x0,y0,z0,pos+dl,side)) return
		if( (side.eq.2 .or. side.eq.4) .and. outside_subdomain(x0,y0,z0,pos-dl,side)) return
		
#ifdef gpu
		call RecvFullDomainEMFldFromGPU
#endif		
		
		do k=1,mz
	          do j=1,my
	               do i=1,mx
					   
 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then 
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0 = vz * b_y - vy * b_z 
						  Ex(i,j,k)=(Ex(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,dl)+f0
 					  end if

 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0 = vx * b_z - vz * b_x
						  Ey(i,j,k)=(Ey(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,dl)+f0
 					  end if 

 					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then 
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0= vy * b_x - vx * b_y
						  Ez(i,j,k)=(Ez(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,dl)+f0
 					  end if 
			  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=(Bx(i,j,k)-b_x)*AttenFactor(xg,yg,zg,pos,side,dl)+b_x
 					  end if 

 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=(By(i,j,k)-b_y)*AttenFactor(xg,yg,zg,pos,side,dl)+b_y
 					  end if 

 					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=(Bz(i,j,k)-b_z)*AttenFactor(xg,yg,zg,pos,side,dl)+b_z
 					  end if 
		
	               end do
	          end do
	    end do
		
#ifdef gpu
	    call SendFullDomainEMFldtoGPU
#endif			
	
	end subroutine AttenuateFld
	
	function AttenFactor(x,y,z,pos,side,dl)
		integer   :: side
		real(dbpsn) :: x,y,z,pos,dl
		real(psn) :: AttenFactor
		
		select case (side)
			case(1)
				AttenFactor = (abs(x-pos)/dl)**3
			case(2)
	            AttenFactor = (abs(x-pos)/dl)**3
			case(3)
			    AttenFactor = (abs(y-pos)/dl)**3
			case(4)
			    AttenFactor = (abs(y-pos)/dl)**3
		end select
		AttenFactor = min(1.0_psn,AttenFactor)
		AttenFactor = max(0.0_psn,AttenFactor)
	end function AttenFactor
	
	logical function atnn_domain(xg,yg,zg,pos,side,dl)
		real(dbpsn) :: xg,yg,zg,pos,dl
		integer :: side
		atnn_domain = .false.
		select case (side)
			case(1)
				if(xg.ge.pos.and.xg.le.pos+dl) atnn_domain = .true.
			case(2)
				if(xg.le.pos.and.xg.ge.pos-dl) atnn_domain = .true.
			case(3)
			    if(yg.ge.pos.and.yg.le.pos+dl) atnn_domain = .true.
			case(4)
				if(yg.le.pos.and.yg.ge.pos-dl) atnn_domain = .true.
		end select
	end function atnn_domain
	
	real(dbpsn) function inflow_pos_prtl(side)
		integer :: side
		inflow_pos_prtl = 0 
		if(side.eq.1) inflow_pos_prtl = BC_Xmin_Prtl
		if(side.eq.2) inflow_pos_prtl = BC_Xmax_Prtl
		if(side.eq.3) inflow_pos_prtl = BC_Ymin_Prtl
		if(side.eq.4) inflow_pos_prtl = BC_Ymax_Prtl
    end function inflow_pos_prtl
	
	real(dbpsn) function inflow_pos_fld(side)
		integer :: side
		inflow_pos_fld = 0 
		if(side.eq.1) inflow_pos_fld = BC_Xmin_Fld
		if(side.eq.2) inflow_pos_fld = BC_Xmax_Fld
		if(side.eq.3) inflow_pos_fld = BC_Ymin_Fld
		if(side.eq.4) inflow_pos_fld = BC_Ymax_Fld
    end function inflow_pos_fld
	
	! update the injector position
	subroutine updateInflowBC_position_fld(side)
		integer :: side
		real(psn) :: shift
		
		if(inflowBC_speed(side).eq.0) return
		shift = inflowBC_speed(side)*c
		
		select case (side)
		    case(1)
				BC_Xmin_Fld = BC_Xmin_Fld + shift
		    case(2)
				BC_Xmax_Fld = BC_Xmax_Fld + shift
		    case(3)
				BC_Ymin_Fld = BC_Ymin_Fld + shift
		    case(4)
				BC_Ymax_Fld = BC_Ymax_Fld + shift
		end select
	end subroutine updateInflowBC_position_fld
	
	subroutine updateInflowBC_position_prtl
		
			BC_Xmin_Prtl = BC_Xmin_Prtl + inflowBC_speed(1)*c
			BC_Xmax_Prtl = BC_Xmax_Prtl + inflowBC_speed(2)*c
			BC_Ymin_Prtl = BC_Ymin_Prtl + inflowBC_speed(3)*c
			BC_Ymax_Prtl = BC_Ymax_Prtl + inflowBC_speed(4)*c
	
	end subroutine updateInflowBC_position_prtl
	
!--------------------------------------------------------------------------------------------
!
!        Adjust Position of the injector automatically
!
!--------------------------------------------------------------------------------------------	
	subroutine  AdjustInjectorPositionFromShock(Start, Freq, MinDist, MaxDist, MinInjSpeed, BpeakXmin)
		integer :: Start, Freq, MinDist, MaxDist, BpeakXmin
		real(psn) :: MinInjSpeed, x1, x2, dist, x3
		integer :: Xshock 
		
		if(modulo(t,Freq).ne.0) return
		if(t.lt.Start) return
		
		call FindShockXpos_BmagPeak(Xshock, BpeakXmin)
		
		dist = BC_Xmax_Prtl - Xshock
		x3 = (MaxDist - MinDist)/3.0
		x1 = MinDist + x3
		x2 = MaxDist - x3
		
		if(dist .lt. x1) inflowBC_speed(2) = inflowBC_speed(2) + (x1 - dist)/x3 
		if(dist .gt. x2) inflowBC_speed(2) = inflowBC_speed(2) - (dist - x2)/x3
		
		inflowBC_speed(2) = min(inflowBC_speed(2), 1.0)
		inflowBC_speed(2) = max(inflowBC_speed(2), MinInjSpeed) 
	end subroutine  AdjustInjectorPositionFromShock

	
	
!--------------------------------------------------------------------------------------------
!
!        Some Auxiliary Subroutines useful for the shock problem 
!
!--------------------------------------------------------------------------------------------

	!the following subroutines tries to find the current location of the shock by following total magnetic energy vs. x 
	subroutine FindShockXpos_BmagPeak(Xshock, BpeakXmin)
		implicit none 
		integer :: Xshock, BpeakXmin, BpeakXmin_local
	    integer :: i,j,k,i1,j1,k1
		integer :: xpeak
		real(psn) :: mag_peak
	    real(psn), dimension(mx) :: mag1D_slice,mag1D_this_proc
		real(psn), dimension(0:nproc-1) :: mag_peak_all_proc
		integer, dimension(0:nproc-1)   :: xpeak_all_proc
		type(MPI_Status) :: mpi_stat    
		integer :: mpi_err
		
		Xshock = 0 

#ifdef gpu
		Bx = Bx_gpu; By=By_gpu; Bz=Bz_gpu;
#endif		
		
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
		 
		 ! make sure that any x less than BpeakXmin are not included in finding the peak
		 BpeakXmin_local = BpeakXmin -xborders(procxind(proc)) +3
		 do i=3,mx-3
			 if(i.lt.BpeakXmin_local) mag1D_this_proc(i) = 0 
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
	
!-----------------------------------------------------------------------------------------
!
!      Domain size adjustments
!	
!-----------------------------------------------------------------------------------------		

subroutine EnlargeDomain(side)
	integer :: side ! only side = 2 is implemented
	
	if(side.eq.2) call EnlargeDomainRight
end subroutine EnlargeDomain

subroutine EnlargeDomainRight
	integer :: inc
	
	inc = 0 
	if(BC_Xmax_Fld+4.gt.xborders(nSubDomainsX)) then 
		inc = 4*int(2.0*(xborders(nSubDomainsX)-xborders(nSubDomainsX-1))/4.0) !mutiple of 4
		inc = min(512,max(128,inc))
	end if

	
	if(inc.gt.0 .and. procxind(proc).eq.nSubDomainsX-1) then 
#ifdef gpu
        call RecvFullDomainEMFldFromGPU
#endif		
		call Enlarge_EM_Fld_X(inc)
	
		mx = mx + inc
		
		!Initialise EM Fld in the newly created domain
		call InflowBC_SetMagFld(BC_Xmax_Fld,2)
		call InflowBC_SetElcFld(BC_Xmax_Fld,2)
	
		call InitPrtlBoundaries
		call ReshapeAuxFld
		call ReshapeShortMoverFldArr
#ifdef gpu
		call ResetGridGPU
#endif 		
	end if 
	
    xborders(nSubDomainsX)=xborders(nSubDomainsX)+inc
	fdataxf_box=xborders(nSubDomainsX)
	  
end subroutine EnlargeDomainRight	

subroutine Enlarge_EM_Fld_X(inc)
     integer :: inc
     call EnlargeFldArrX(Ex,inc,mx,my,mz)
     call EnlargeFldArrX(Ey,inc,mx,my,mz)
     call EnlargeFldArrX(Ez,inc,mx,my,mz)
     call EnlargeFldArrX(Bx,inc,mx,my,mz)
     call EnlargeFldArrX(By,inc,mx,my,mz)
     call EnlargeFldArrX(Bz,inc,mx,my,mz)
	 if(nMoverEMfilter.gt.0) then 
	     call EnlargeFldArrX(filteredEx,inc,mx,my,mz)
	     call EnlargeFldArrX(filteredEy,inc,mx,my,mz)
	     call EnlargeFldArrX(filteredEz,inc,mx,my,mz)
	 end if
end subroutine Enlarge_EM_Fld_X

subroutine EnlargeFldArrX(Fld,inc,sx_old,sy_old,sz_old)
	 integer     :: inc, sx_old,sy_old,sz_old
	 real(psn), dimension(:,:,:), allocatable :: Fld
	 real(psn), dimension(:,:,:),allocatable :: FldTemp
     allocate(FldTemp(sx_old+abs(inc),sy_old,sz_old))  !negative inc is interpreted as enlargment to the left side of the box
	 FldTemp=0.0_psn   
	 if(inc.gt.0) then 
		 FldTemp(1:sx_old,:,:)=Fld(1:sx_old,:,:)
	 else 
		 FldTemp(abs(inc)+1:sx_old+abs(inc),:,:)=Fld(1:sx_old,:,:)
	 end if 
	 deallocate(Fld)
     call move_alloc(FldTemp,Fld)
end subroutine EnlargeFldArrX
	
	
! older injection algorithm: particles are reflected in the fluid frame	

! 	subroutine InjectNewPrtlY(side,yinj)
! 		integer :: side !1=upper, 2=lower
! 		real(psn) :: yinj,y1,y2,dy,dely,yinj_local
! 		real(psn) :: xglobal,yglobal,zglobal
! 		integer :: nprtl_new,n,tag
! 		real(psn) :: r1,r2,r3,rnd_acpt,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
! 		real(psn) :: vx,vy,vz
! 		real(psn) :: Temp
!
! 	   	 !new particles will be injected between y1 and y2, but only in the region emptied out by the older drifting particles
! 	   	 yinj_local=yinj-yborders(procyind(proc))+3
! 		 if(side.eq.1) then
! 		     y1=max(yinj-vmax_inj_upper,real(yborders(procyind(proc))))
! 	   	     y2=min(yinj,real(yborders(procyind(proc)+1)))
! 			 ydrift=>ydrift_upper
! 			 InjDen=>InjDenUpper
! 		 end if
! 	   	 if(side.eq.2) then
! 		     y1=max(yinj,real(yborders(procyind(proc))))
! 	   	     y2=min(yinj+vmax_inj_lower,real(yborders(procyind(proc)+1)))
! 			 ydrift=>ydrift_lower
! 			 InjDen=>InjDenLower
! 		 end if
!
!
! 	   	 y1=y1-yborders(procyind(proc))+3 !switch to the local cordinate
! 	   	 y2=y2-yborders(procyind(proc))+3
! 		 dy=max(0.0_psn,y2-y1)
! 		 nprtl_new= epc*dy*(xmax-xmin)*(zmax-zmin)
! 		 if(nprtl_new.eq.0) return
! 		 if(nprtl_new.lt.1000) call GetIntPoissonDist(real(real(nprtl_new),psn),nprtl_new)
! 		 if(used_prtl_arr_size+2*nprtl_new.gt.prtl_arr_size) call ReshapePrtlArr(prtl_arr_size+3*nprtl_new)
!
! 		 do n=1,nprtl_new
! 	 		call random_number(r1)
! 	 	    call random_number(r2)
! 	 		call random_number(r3)
! 	 		xlocal=xmin+r1*(xmax-xmin)
! 			ylocal=y1+r2*(y2-y1)
! 			zlocal=zmin+r3*(zmax-zmin)
!
! 			xglobal=xlocal+xborders(procxind(proc))-3
! 			yglobal=ylocal+yborders(procyind(proc))-3
! #ifndef twoD
! 			zglobal=zlocal+zborders(proczind(proc))-3
! #else
!             zglobal=zlocal
! #endif
! 			call ydrift(xglobal,yglobal,zglobal,vx,vy,vz)
! 			if(side.eq.1) dely=ylocal-(yinj_local + vy*c)
! 			if(side.eq.2) dely=(yinj_local + vy*c)-ylocal
!
! 			if(dely.ge.0) then
! 				call random_number(rnd_acpt)
! 				if(rnd_acpt.le.InjDen(xglobal,yglobal,zglobal)) then
! 					Temp=InjTempIon(xglobal,yglobal,zglobal)
! 					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
! 					call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
! 					call GetPrtlTag(tag,1)
! 					call InsertParticleAt(used_prtl_arr_size+1,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,0.0_psn)
! 					Temp=InjTempElc(xglobal,yglobal,zglobal)
! 					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
! 					call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
! 					call GetPrtlTag(tag,2)
! 					call InsertParticleAt(used_prtl_arr_size+2,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,0.0_psn)
! 					used_prtl_arr_size=used_prtl_arr_size+2
! 					np=np+2
! 			    end if
! 			end if
! 		 end do
! 	end subroutine InjectNewPrtlY




!
! subroutine RefPrtlFluidFrameUpperY(yinj) !reflect the backward-moving particles in the fluid frame
! 	real(psn) :: yinj,y1,vx,vy,vz,x0,y0,z0
! 	real(psn) :: ox,oy,oz
! 	integer :: n
!
! 	if( (yinj-yborders(procyind(proc))+3).gt.my-1) return
!
! 	ox=xborders(procxind(proc))-3
! 	oy=yborders(procyind(proc))-3
! 	oz=zborders(proczind(proc))-3
!
! 	do n=1,used_prtl_arr_size
! 		call ydrift_upper(xp(n)+ox,yp(n)+oy,zp(n)+oz,vx,vy,vz)
! 		y1= yinj+vy*c-yborders(procyind(proc))+3 ! instantaneous location of the drifted boundary layer
! 		if((qp(n).ne.0).and.(yp(n).gt.y1)) then
! 			call ReflectPrtlMovingFrameY(vy,up(n),vp(n),wp(n))
! 			x0=xp(n)
! 			y0=yp(n)
! 			z0=zp(n)
! 			yp(n)=2.0_psn*y1-y0
! 			call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
! 		end if
! 	end do
! end subroutine RefPrtlFluidFrameUpperY
!
! subroutine RefPrtlFluidFrameLowerY(yinj) !reflect the backward-moving particles in the fluid frame
! 	real(psn) :: yinj,y1,vx,vy,vz,x0,y0,z0
! 	real(psn) :: ox,oy,oz
! 	integer :: n
!
! 	if( (yinj-yborders(procyind(proc))+3).lt.2) return
!
! 	ox=xborders(procxind(proc))-3
! 	oy=yborders(procyind(proc))-3
! 	oz=zborders(proczind(proc))-3
!
! 	do n=1,used_prtl_arr_size
! 		call ydrift_lower(xp(n)+ox,yp(n)+oy,zp(n)+oz,vx,vy,vz)
! 		y1= yinj+vy*c-yborders(procyind(proc))+3 ! instantaneous location of the drifted boundary layer
! 		if((qp(n).ne.0).and.(yp(n).lt.y1)) then
! 			call ReflectPrtlMovingFrameY(vy,up(n),vp(n),wp(n))
! 			x0=xp(n)
! 			y0=yp(n)
! 			z0=zp(n)
! 			yp(n)=2.0_psn*y1-y0
! 			call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
! 		end if
! 	end do
! end subroutine RefPrtlFluidFrameLowerY
!
!
!
! subroutine ReflectPrtlMovingFrameY(Drift,ugamma,vgamma,wgamma)
! 	real(psn) :: Drift,DriftGamma,DriftBeta
! 	real(psn) :: ugamma,vgamma,wgamma
! 	real(psn) :: gamma
!
! 	if(abs(Drift).lt.1.0) then
!          DriftBeta=Drift
!          DriftGamma=1.0_psn/sqrt((1.0_psn-Drift)*(1.0_psn+Drift))
!     else
!          DriftGamma=abs(Drift)
!          DriftBeta=sqrt((Drift-1.0_psn)*(Drift+1.0_psn))/Drift
!     end if
!
! 	gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
! 	DriftBeta=-DriftBeta
! 	vgamma=DriftGamma*vgamma + DriftGamma*DriftBeta*gamma!in the drifting frame
!
! 	vgamma=-vgamma
! 	gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
! 	DriftBeta=-DriftBeta
! 	vgamma=DriftGamma*vgamma + DriftGamma*DriftBeta*gamma! back to the simualtion frame
! end subroutine ReflectPrtlMovingFrameY



! 	subroutine ClearOldPrtlLeft(xinj)
! 		real(dbpsn) :: xinj
! 		real(psn) :: xinj_local
! 		integer :: n,count
!
!    	    xinj_local=xinj-xborders(procxind(proc))+3
! 		if(xinj_local.lt.1) return
! 		count=0
! 		do n=1,used_prtl_arr_size
! 			if(qp(n).eq.0) cycle
! 			if(xp(n).le.xinj_local) then
! 				call DeletePrtl(n)
! 			    count=count+1
! 			end if
! 		end do
! 		np=np-count
! 	end subroutine ClearOldPrtlLeft
! 	subroutine ClearOldPrtlRight(xinj)
! 		real(dbpsn) :: xinj
! 		real(dbpsn) :: xinj_local
! 		integer :: n,count
!
!    	    xinj_local=xinj-xborders(procxind(proc))+3
! 		if(xinj_local.gt.mx) return
!
! #ifdef gpu
!         call ClearOldPrtlRightGPU(xinj_local)
! #endif
!
! 		count=0
! 		do n=1,used_prtl_arr_size
! 			if(qp(n).eq.0) cycle
! 			if(xp(n).gt.xinj_local) then
! 				call DeletePrtl(n)
! 				count=count+1
! 			end if
! 		end do
! 		np=np-count
! 	end subroutine ClearOldPrtlRight
!
! 	subroutine ClearOldPrtlBottom(yinj)
! 		real(dbpsn) :: yinj
! 		real(psn)   :: yinj_local
! 		integer :: n,count
!
!    	    yinj_local=yinj-yborders(procyind(proc))+3
! 		if(yinj_local.lt.1) return
! 		count=0
! 		do n=1,used_prtl_arr_size
! 			if(qp(n).eq.0) cycle
! 			if(yp(n).le.yinj_local) then
! 				call DeletePrtl(n)
! 			    count=count+1
! 			end if
! 		end do
! 		np=np-count
! 	end subroutine ClearOldPrtlBottom
! 	subroutine ClearOldPrtlTop(yinj)
! 		real(dbpsn) :: yinj
! 		real(psn)   :: yinj_local
! 		integer :: n,count
!
!    	    yinj_local=yinj-yborders(procyind(proc))+3
! 		if(yinj_local.gt.my) return
! 		count=0
! 		do n=1,used_prtl_arr_size
! 			if(qp(n).eq.0) cycle
! 			if(yp(n).ge.yinj_local) then
! 				call DeletePrtl(n)
! 				count=count+1
! 			end if
! 		end do
! 		np=np-count
! 	end subroutine ClearOldPrtlTop
	
end module injector
