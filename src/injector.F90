module injector 
	use parameters
	use vars
	use deposit
	use prob_dist
	use memory
	use initialise
	implicit none 
	
	integer :: nInjPrtl=0
	
	type InflowFldType
		real(psn) :: Attenuate  
		procedure(fld_ext), pointer, nopass  :: MagFld =>null()
		procedure(fld_ext), pointer, nopass :: Drift =>null()	
	end type InflowFldType
	
	type InjPrtlType
		integer :: side=0 ! 1: x,left; 2: x,right; 3: y,bottom; 4: y,top  
		integer :: flvr1,flvr2
		integer :: dist_type1, dist_type2! 1: thermal; 2: isotropic, use PDF table;  
		real(psn) :: vmax
		
		procedure(func_ext), pointer, nopass :: Den =>null()
		procedure(fld_ext), pointer, nopass :: Drift =>null()
		procedure(func_ext), pointer, nopass :: Temp1 =>null()
		procedure(func_ext), pointer, nopass :: Temp2 =>null()
		procedure(func_ext), pointer, nopass :: SpeedDist1 =>null()
		procedure(func_ext), pointer, nopass :: SpeedDist2 =>null()
		
		integer :: TableSize=10000
		real(psn), dimension(:), allocatable :: Table1, PDF_Table1
		real(psn), dimension(:), allocatable :: Table2, PDF_Table2
	
	end type InjPrtlType
	
	type(InjPrtlType), dimension(:), allocatable :: InjPrtl
	type(InflowFldType), dimension(6) :: InflowFld ! 1: x,left; 2: x,right; 3: y,bottom; 4: y,top 
	real(psn), dimension(6) :: inflowBC_speed = 0
	
	abstract interface
	    function func_ext(x,y,z)
	 	     import :: psn
			 real(psn) :: x,y,z
			 real(psn) :: func_ext
	 	end function
		subroutine fld_ext(x,y,z,fx,fy,fz)
			import :: psn
			real(psn) :: x,y,z,fx,fy,fz
		end subroutine
	end interface 
	
contains 
!------------------------------------------	
! Particles
!------------------------------------------	
	
    subroutine InflowBC_Prtl
		integer :: n 
		real(psn) :: pos
		
		do n=1,nInjPrtl
			pos = inflow_pos(InjPrtl(n)%side)

			if(InjPrtl(n)%side.eq.1) call ClearOldPrtlLeft(pos)
			if(InjPrtl(n)%side.eq.2) call ClearOldPrtlRight(pos)
			if(InjPrtl(n)%side.eq.3) call ClearOldPrtlBottom(pos)
			if(InjPrtl(n)%side.eq.4) call ClearOldPrtlTop(pos)
			
			call InjectNewPrtlPair(InjPrtl(n))
			
		end do
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
	
	subroutine ClearOldPrtlLeft(xinj)
		real(psn) :: xinj,xinj_local
		integer :: n,count 
		
   	    xinj_local=xinj-xborders(procxind(proc))+3
		if(xinj_local.lt.1) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(xp(n).le.xinj_local) then 
				call DeletePrtl(n)
			    count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlLeft
	subroutine ClearOldPrtlRight(xinj)
		real(psn) :: xinj,xinj_local
		integer :: n,count 
		
   	    xinj_local=xinj-xborders(procxind(proc))+3
		if(xinj_local.gt.mx) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(xp(n).ge.xinj_local) then 
				call DeletePrtl(n)
				count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlRight
	
	subroutine ClearOldPrtlBottom(yinj)
		real(psn) :: yinj,yinj_local
		integer :: n,count 
		
   	    yinj_local=yinj-yborders(procyind(proc))+3
		if(yinj_local.lt.1) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(yp(n).le.yinj_local) then 
				call DeletePrtl(n)
			    count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlBottom
	subroutine ClearOldPrtlTop(yinj)
		real(psn) :: yinj,yinj_local
		integer :: n,count 
		
   	    yinj_local=yinj-yborders(procyind(proc))+3
		if(yinj_local.gt.my) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(yp(n).ge.yinj_local) then 
				call DeletePrtl(n)
				count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlTop
	
	
	subroutine InjectNewPrtlPair(inj)
		type(InjPrtlType) :: inj
		real(psn) :: x1,x2,y1,y2,z1,z2,q1,q2
		real(psn) :: xglobal,yglobal,zglobal
		integer :: nprtl_new,n,tag
		real(psn) :: rnd_acpt,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
		real(psn) :: vx,vy,vz
		real(psn) :: Temp

        q1=FlvrCharge(inj%flvr1)
		q2=FlvrCharge(inj%flvr2)
		x1 = xmin; x2 = xmax; y1 = ymin; y2 = ymax; z1 = zmin; z2 = zmax 
		call SetInjDomain(inj,nprtl_new,x1,x2,y1,y2,z1,z2)
		 
		 do n=1,nprtl_new
			 
		    call GetNewPrtlPos(x1,x2,y1,y2,z1,z2,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
			call inj%Drift(xglobal,yglobal,zglobal,vx,vy,vz)
			
			call random_number(rnd_acpt)
			if(rnd_acpt.le.inj%Den(xglobal,yglobal,zglobal)) then	
				
				if(inj%dist_type1.eq.1) then 
					Temp=inj%Temp1(xglobal,yglobal,zglobal)
					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
				else if(inj%dist_type1.eq.2) then 
					call GetIsoVelGammaTable(inj%TableSize,inj%Table1,inj%PDF_Table1,ugamma,vgamma,wgamma) 
			    end if  
				call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
				call GetPrtlTag(tag,inj%flvr1)
				call InsertParticleAt(used_prtl_arr_size+1,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,q1,tag,inj%flvr1,0.0_psn)
				
				
				if(inj%dist_type2.eq.1) then 
					Temp=inj%Temp2(xglobal,yglobal,zglobal)
					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
				else if(inj%dist_type2.eq.2) then 
					call GetIsoVelGammaTable(inj%TableSize,inj%Table2,inj%PDF_Table2,ugamma,vgamma,wgamma) 
			    end if   
				call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
				call GetPrtlTag(tag,inj%flvr2)
				call InsertParticleAt(used_prtl_arr_size+2,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,q2,tag,inj%flvr2,0.0_psn)
				
				used_prtl_arr_size=used_prtl_arr_size+2
				np=np+2
				
			 end if 
		 end do 
	end subroutine InjectNewPrtlPair
	
	subroutine SetInjDomain(inj,nprtl_new,x1,x2,y1,y2,z1,z2)
		type(InjPrtlType) :: inj
		real(psn) :: x1,x2,y1,y2,z1,z2
		real(psn) :: pos, vmax
		integer   :: nprtl_new
		integer   :: side
		
		pos = inflow_pos(inj%side)
		vmax = inj%vmax
		side = inj%side
		
	   	if(side.eq.1) then 
		     x1=max(pos-vmax,real(xborders(procxind(proc)))) -xborders(procxind(proc))+3
	   	     x2=min(pos,real(xborders(procxind(proc)+1))) -xborders(procxind(proc))+3
		end if
		if(side.eq.2) then 
		     x1=max(pos,real(xborders(procxind(proc)))) -xborders(procxind(proc))+3
	   	     x2=min(pos+vmax,real(xborders(procxind(proc)+1))) -xborders(procxind(proc))+3	 
		end if
	   	if(side.eq.3) then 
		     y1=max(pos-vmax,real(yborders(procyind(proc)))) -yborders(procyind(proc))+3
	   	     y2=min(pos,real(yborders(procyind(proc)+1))) -yborders(procyind(proc))+3
		end if
		if(side.eq.4) then 
		     y1=max(pos,real(yborders(procyind(proc)))) -yborders(procyind(proc))+3
	   	     y2=min(pos+vmax,real(yborders(procyind(proc)+1))) -yborders(procyind(proc))+3	 
		end if

		nprtl_new= epc*max(0.0_psn,x2-x1)*max(0.0_psn,y2-y1)*max(0.0_psn,z2-z1)
		if(nprtl_new.eq.0) return 
		if(nprtl_new.lt.1000) call GetIntPoissonDist(real(real(nprtl_new),psn),nprtl_new)
		if(used_prtl_arr_size+2*nprtl_new.gt.prtl_arr_size) call ReshapePrtlArr(prtl_arr_size+3*nprtl_new)
		
	end subroutine SetInjDomain
	
	subroutine GetNewPrtlPos(x1,x2,y1,y2,z1,z2,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
		real(psn) :: x1,x2,y1,y2,z1,z2
		real(psn) :: xlocal,ylocal,zlocal,xglobal,yglobal,zglobal
		real(psn) :: r1,r2,r3
	 		call random_number(r1)
	 	    call random_number(r2)
	 		call random_number(r3)
	 		xlocal=x1+r1*(x2-x1)
			ylocal=y1+r2*(y2-y1)
			zlocal=z1+r3*(z2-z1)
			
			xglobal=xlocal+xborders(procxind(proc))-3
			yglobal=ylocal+yborders(procyind(proc))-3
#ifndef twoD			
			zglobal=zlocal+zborders(proczind(proc))-3
#else
            zglobal=zlocal
#endif 	

	end subroutine GetNewPrtlPos
			

	!------------------------------------------------------------
	!subroutine to genrate unique tag for particles, warning:: tagging start from 1 if simulation is restarted
	!may not work in some cases when large number of particles are tagged 
	!------------------------------------------------------------
	subroutine GetPrtlTag(tag,FlvID)
		implicit none
		integer :: tag,FlvID
	
		if(TagCounter(FlvID).gt.psave_ratio) TagCounter(FlvID)=1
		if(TagCounter(FlvID).eq.1) then  
			tag=CurrentTagID(FlvID)
		    if(mod(CurrentTagID(FlvID),NtagProcLen).eq.0) CurrentTagID(FlvID)=CurrentTagID(FlvID)+NtagProcLen*(nproc-1)
	        CurrentTagID(FlvID)=CurrentTagID(FlvID)+1
	    else
			tag=0
		end if	
		TagCounter(FlvID)=TagCounter(FlvID)+1
	end subroutine GetPrtlTag

!------------------------------------------	
! Fields
!------------------------------------------	

	subroutine InflowBC_SetFld(flw,side)
		type(InflowFldType) :: flw
		real(psn):: b_x,b_y,b_z,vx,vy,vz,x0,y0,z0,xg,yg,zg
		real(psn):: pos
		integer  :: side
		logical  :: outside
	    integer :: i,j,k
		
		!enlarge the domain size if needed
		if(inflowBC_speed(side).ne.0)  call EnlargeDomain(side)
		
		pos = inflow_pos(side)

		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif	
	
		 do k=1,mz
	          do j=1,my
				  do i=1,mx
				  
					  xg = real(i+0.5_psn,psn)+x0; yg = real(j,psn)+y0; zg = real(k,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ex(i,j,k)= vz * b_y - vy * b_z 
					  end if

					  xg = real(i,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ey(i,j,k)= vx * b_z - vz * b_x
					  end if 

					  xg = real(i,psn)+x0; yg = real(j,psn)+y0; zg = real(k+0.5+psn,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call flw%Drift(xg,yg,zg,vx,vy,vz)
						  Ez(i,j,k)= vy * b_x - vx * b_y
					  end if 
			  
					  xg = real(i,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k+0.5_psn,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=b_x
					  end if 

					  xg = real(i+0.5_psn,psn)+x0; yg = real(j,psn)+y0; zg = real(k+0.5_psn,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=b_y
					  end if 

					  xg = real(i+0.5_psn,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k,psn)+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=b_z
					  end if 
		  
				  end do
			  end do  
	     end do
		 
		 if(flw%Attenuate.gt.0) call AttenuateFld(flw,side) 
		 
		 if(inflowBC_speed(side).ne.0) call updateInflowBC_position(side) ! update the injector position, assuming no further actions are required
			
	end subroutine InflowBC_SetFld
	
	logical function exterior(xg,yg,zg,pos,side)
		real(psn) :: xg,yg,zg,pos
		integer :: side
		exterior = .false.
		if(side.eq.1 .and. xg.le. pos) exterior = .true. 
        if(side.eq.2 .and. xg.ge. pos) exterior = .true. 
		if(side.eq.3 .and. yg.le. pos) exterior = .true. 
	    if(side.eq.4 .and. yg.ge. pos) exterior = .true. 
	end function exterior
	
	!------------------------------------------------------------------------
	! Attenuation of the EM field componnets
	!------------------------------------------------------------------------
	subroutine AttenuateFld(flw,side)
		type(InflowFldType) :: flw
		integer  :: side
		real(psn):: b_x,b_y,b_z,vx,vy,vz,x0,y0,z0,xg,yg,zg, f0
		real(psn):: pos, dl
		real(psn) :: AttenFactor=0.9_psn
	    integer :: i,j,k
		
		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif	
        pos = inflow_pos(side)
		dl = flw%Attenuate
		
		do k=1,mz
	          do j=1,my
	               do i=1,mx
					   
 					  xg = real(i+0.5_psn,psn)+x0; yg = real(j,psn)+y0; zg = real(k,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then 
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0 = vz * b_y - vy * b_z 
						  Ex(i,j,k)=(Ex(i,j,k)-f0)*AttenFactor+f0
 					  end if

 					  xg = real(i,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0 = vx * b_z - vz * b_x
						  Ey(i,j,k)=(Ey(i,j,k)-f0)*AttenFactor+f0
 					  end if 

 					  xg = real(i,psn)+x0; yg = real(j,psn)+y0; zg = real(k+0.5+psn,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then 
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0= vy * b_x - vx * b_y
						  Ez(i,j,k)=(Ez(i,j,k)-f0)*AttenFactor+f0
 					  end if 
			  
 					  xg = real(i,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k+0.5_psn,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=(Bx(i,j,k)-b_x)*AttenFactor+b_x
 					  end if 

 					  xg = real(i+0.5_psn,psn)+x0; yg = real(j,psn)+y0; zg = real(k+0.5_psn,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=(By(i,j,k)-b_y)*AttenFactor+b_y
 					  end if 

 					  xg = real(i+0.5_psn,psn)+x0; yg = real(j+0.5_psn,psn)+y0; zg = real(k,psn)+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=(Bz(i,j,k)-b_z)*AttenFactor+b_z
 					  end if 
		
	               end do
	          end do
	    end do
	end subroutine AttenuateFld
	
	logical function atnn_domain(xg,yg,zg,pos,side,dl)
		real(psn) :: xg,yg,zg,pos,dl
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
	
	real(psn) function inflow_pos(side)
		integer :: side
		inflow_pos = 0 
		if(side.eq.1) inflow_pos = BC_Xmin_Prtl
		if(side.eq.2) inflow_pos = BC_Xmax_Prtl
		if(side.eq.3) inflow_pos = BC_Ymin_Prtl
		if(side.eq.4) inflow_pos = BC_Ymax_Prtl
    end function inflow_pos
	
	! update the injector position, assuming no further actions are required
	subroutine updateInflowBC_position(side)
		integer :: side
		real(psn) :: shift
		
		if(inflowBC_speed(side).eq.0) return
		shift = inflowBC_speed(side)*c

		select case (side)
		    case(1)
				BC_Xmin_Prtl = BC_Xmin_Prtl + shift 
		    case(2) 
				BC_Xmax_Prtl = BC_Xmax_Prtl + shift
		    case(3) 
				BC_Ymin_Prtl = BC_Ymin_Prtl + shift
		    case(4)
				BC_Ymax_Prtl = BC_Ymax_Prtl + shift
		end select

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
		
	end subroutine updateInflowBC_position
	
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
		call Enlarge_EM_Fld_X(inc)
	
		mx = mx + inc
	
		call InitPrtlBoundaries
		call ReshapeAuxFld
		call ReshapeShortMoverFldArr
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

	
end module injector