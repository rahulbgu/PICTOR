module bc
    use parameters
	use vars
	use communication
	use mem_prtl
	use init_prtl
	use bc_open
	use bc_adjust_speed
	use bc_attenuate
	use bc_inflow
	use bc_curved_surf
	use bc_reflect
	use subdomains
#ifdef gpu
    use bc_gpu 
#endif	
	implicit none 
contains 	
!---------------------------------------------------------------------------------
!  Subroutines to implement BCs from the setup file 
!---------------------------------------------------------------------------------
    subroutine SetBC_Fld(Axis,Position,Type,Side,AttnThickness,AttnScale)
		 character (len=*) :: Axis,Type,Side
		 real(dbpsn)  :: Position
		 integer :: side_ind
		 real(psn), optional :: AttnThickness, AttnScale
		 
		 side_ind = SideInd(Axis,Side)
		 
		 if(.not.restart) bc_face(side_ind)%pos_fld=Position
		 bc_face(side_ind)%type_fld=Type
		 call PeriodicOff(side_ind)
		 
		 if(Type.eq.'attn') then 
			 attn_em_vacumm = .true.
	 		 bc_face(side_ind)%attn_thickness = 10.0 !default
	 		 bc_face(side_ind)%attn_scale = 4.0 ! default
			 if(present(AttnThickness)) bc_face(side_ind)%attn_thickness = real(AttnThickness,psn)
			 if(present(AttnScale)) bc_face(side_ind)%attn_scale = real(AttnScale,psn)
	 		 call CheckAttnBCPositioning(side_ind)
		 end if

	end subroutine SetBC_Fld	 		

    subroutine SetBC_Prtl(Axis,Position,Type,Side)
		 character (len=*) :: Axis,Type,Side
		 real(dbpsn)  :: Position
		 integer      :: side_ind
		 
		 side_ind = SideInd(Axis,Side)
		 
		 if(.not.restart) bc_face(side_ind)%pos_prtl=Position
		 bc_face(side_ind)%type_prtl=Type
		 call PeriodicOff(side_ind)
	
	end subroutine SetBC_Prtl
	
	subroutine SetBC_InflowPrtl(Axis,Position,Side,InjSpeed,Flvr1,Flvr2,Density,Temperature1,Temperature2,SpeedDist1,SpeedDist2,DriftVelocity,DriftVmax,Vmax1,Vmax2,Multiplicity1,Multiplicity2,pid,Placement,MaxDisplacement)
		character (len=*) :: Axis,Side
		integer, optional :: Flvr1,Flvr2
		integer :: side_ind
		real(dbpsn)  :: Position
		real(psn), optional :: Vmax1 !the maximum particle speed in the plasma frame for Flvr1, default is c 
		real(psn), optional :: Vmax2 ! max. speed for Flvr2; Vmax1 and Vmax2 are only used when speed dist is provided
		real(psn), optional :: InjSpeed
		real(psn), optional :: DriftVmax !the maximum drift speed, helps improve the efficiency of prtl injection process
		procedure(scalar_global), optional :: Temperature1
		procedure(scalar_global), optional :: Temperature2
		procedure(func1D),        optional :: SpeedDist1
		procedure(func1D),        optional :: SpeedDist2
		procedure(scalar_global)           :: Density 
		procedure(vector_global), optional :: DriftVelocity
		integer,                  optional :: Multiplicity1, Multiplicity2 ! num. of prtl. initialised simultaneously at same location 
		real(psn) :: vmax
		
		character (len=4), optional :: Placement
		real(psn), optional :: MaxDisplacement	
		integer, dimension(:), optional :: pid
		integer, dimension(:), allocatable :: pid_this
		integer :: n, i1, i2

		
		inflowBC = .true.
		side_ind = SideInd(Axis,Side) 
		
		call SetBC_Prtl(Axis,Position,'iflw',Side)
		call PeriodicOff(side_ind)
		
		if(.not.allocated(InjPrtl)) then 
			 allocate(InjPrtl(8)) ! change here or allocate manually if a larger size is required 
			 nInjPrtl=0
		end if 
		
		nInjPrtl = nInjPrtl + 1
		
		InjPrtl(nInjPrtl)%side  = side_ind		
		InjPrtl(nInjPrtl)%Den => Density
		InjPrtl(nInjPrtl)%Drift => DriftVelocity
		InjPrtl(nInjPrtl)%prtl_placement = 'rndm' ! default is random placement
		if(present(Placement)) InjPrtl(nInjPrtl)%prtl_placement = Placement
		InjPrtl(nInjPrtl)%max_displacement = 0.0_psn 
		if(present(MaxDisplacement)) InjPrtl(nInjPrtl)%max_displacement = MaxDisplacement
		
		if(present(Flvr1)) then 
			allocate(InjPrtl(nInjPrtl)%pid(1))
			used_PSP_list_size = used_PSP_list_size +1
			i1 = used_PSP_list_size
			call SetPhaseSpaceProperty(PSP_list(i1), Flvr1,Density,Temperature1,SpeedDist1,Vmax1,DriftVelocity,Multiplicity1)
			InjPrtl(nInjPrtl)%pid = (/ i1/)
			
			if(present(Flvr2)) then 
				if(allocated(InjPrtl(nInjPrtl)%pid)) deallocate(InjPrtl(nInjPrtl)%pid)
				allocate(InjPrtl(nInjPrtl)%pid(2))
				used_PSP_list_size = used_PSP_list_size +1
				i2 = used_PSP_list_size
				call SetPhaseSpaceProperty(PSP_list(i2), Flvr2,Density,Temperature2,SpeedDist2,Vmax2,DriftVelocity,Multiplicity2)
				InjPrtl(nInjPrtl)%pid = (/ i1, i2 /)
			end if 
			
		end if 
		
		if(present(pid)) then 
			if(allocated(InjPrtl(nInjPrtl)%pid)) deallocate(InjPrtl(nInjPrtl)%pid)
			allocate(InjPrtl(nInjPrtl)%pid(size(pid)))
			do n=1,size(pid)
				InjPrtl(nInjPrtl)%pid(n) = pid(n)
			end do
		end if
		
		if(present(DriftVelocity)) then 
			bc_face(side_ind)%flw%Drift => DriftVelocity
		else 
			bc_face(side_ind)%flw%Drift => ZeroVecFld
		end if 
		
		bc_face(side_ind)%flw%vmax = 1.0
		if(present(DriftVmax)) bc_face(side_ind)%flw%vmax = DriftVmax		
		
		if(.not.restart) then 
			bc_face(side_ind)%speed = 0
			if(present(InjSpeed)) bc_face(side_ind)%speed = InjSpeed
		end if

	end subroutine SetBC_InflowPrtl
	
	subroutine ZeroVecFld(x,y,z,vx,vy,vz)
		real(dbpsn) :: x,y,z
		real(psn)   :: vx,vy,vz
		vx = 0.0_psn
		vy = 0.0_psn
		vz = 0.0_psn
	end subroutine ZeroVecFld 
	
	subroutine SetBC_InflowFld(Axis,Position,Side,InjSpeed,DriftVelocity,MagFld,DriftVmax,AttnThickness,AttnScale)
		character (len=*) :: Axis,Side
		integer :: side_ind
		real(dbpsn)  :: Position
		real(psn), optional :: DriftVmax !the maximum drift speed, helps improve the efficiency of prtl injection process
	    real(psn), optional :: AttnThickness, AttnScale
		real(psn), optional :: InjSpeed
		procedure(vector_global), optional :: MagFld
		procedure(vector_global), optional :: DriftVelocity
	
		inflowBC = .true.	
		side_ind = SideInd(Axis,Side) 
		
		call SetBC_Fld(Axis,Position,'iflw',Side)
		call PeriodicOff(side_ind)

		if(present(MagFld)) then 
			bc_face(side_ind)%flw%MagFld => MagFld
		else 
			bc_face(side_ind)%flw%MagFld => ZeroVecFld
		end if
		
		if(present(DriftVelocity)) then 
			bc_face(side_ind)%flw%Drift => DriftVelocity
		else 
			bc_face(side_ind)%flw%Drift => ZeroVecFld
		end if
		
		
		bc_face(side_ind)%flw%vmax = 1.0
		if(present(DriftVmax)) bc_face(side_ind)%flw%vmax = DriftVmax
		
		if(.not.restart) then 
			if(present(InjSpeed)) bc_face(side_ind)%speed = InjSpeed
	    end if 
		
		bc_face(side_ind)%attn_thickness = 1.0 ! default 
		bc_face(side_ind)%attn_scale = 4.0 ! default
	    if(present(AttnThickness)) bc_face(side_ind)%attn_thickness = real(AttnThickness,psn)
	    if(present(AttnScale)) bc_face(side_ind)%attn_scale = real(AttnScale,psn)
		call CheckAttnBCPositioning(side_ind)
	    
	end subroutine SetBC_InflowFld
	
	subroutine SetBoundarySpeed(Axis,Side,Speed)
		character (len=*) :: Axis,Side
		real(psn) :: Speed
		integer :: side_ind
		
		side_ind = SideInd(Axis,Side) 
		bc_face(side_ind)%speed = Speed
		
	end subroutine SetBoundarySpeed
	
	
	integer function SideInd(Axis,Side)
		 character (len=*) :: Axis,Side
		 SideInd=0
		 if(Axis.eq.'x') then
			 if(Side.eq.'lower') SideInd=1 
			 if(Side.eq.'upper') SideInd=2 
		 end if 
		 if(Axis.eq.'y') then
			 if(Side.eq.'lower') SideInd=3 
			 if(Side.eq.'upper') SideInd=4 
		 end if 
		 if(Axis.eq.'z') then
			 if(Side.eq.'lower') SideInd=5 
			 if(Side.eq.'upper') SideInd=6 
		 end if 
	end function SideInd

!---------------------------------------------------------------------------------
!  Essential rotuines to enforece BCs from the main loop 
!---------------------------------------------------------------------------------	
	subroutine EnforceBC_PostMovDep

		if(inflowBC) call InflowBC_Prtl
				
		if(bc_face(4)%type_prtl.eq.'refl') call RefBC_Prtl_Top(bc_face(4)%pos_prtl)
		if(bc_face(3)%type_prtl.eq.'refl') call RefBC_Prtl_Bottom(bc_face(3)%pos_prtl)	
		if(bc_face(2)%type_prtl.eq.'refl') call RefBC_Prtl_Right(bc_face(2)%pos_prtl)
		if(bc_face(1)%type_prtl.eq.'refl') call RefBC_Prtl_Left(bc_face(1)%pos_prtl)	
		
		select case (bc_face(1)%type_prtl)
			case('refl')
				call RefBC_Prtl_Left(bc_face(1)%pos_prtl)
			case('remv')
				call RemovePrtl_Left( bc_pos_local(1, bc_face(1)%pos_prtl), qp,flvp,tagp,xp,yp,zp,up,vp,wp)	
		end select
		
		select case (bc_face(2)%type_prtl)
			case('refl')
				call RefBC_Prtl_Right(bc_face(2)%pos_prtl)
			case('remv')
				call RemovePrtl_Right( bc_pos_local(2, bc_face(2)%pos_prtl), qp,flvp,tagp,xp,yp,zp,up,vp,wp)
		end select
		
		select case (bc_face(3)%type_prtl)
			case('refl')
				call RefBC_Prtl_Bottom(bc_face(3)%pos_prtl)
			case('remv')
				call RemovePrtl_Left( bc_pos_local(3, bc_face(3)%pos_prtl), qp,flvp,tagp,yp,xp,zp,up,vp,wp)
		end select
		
		select case (bc_face(4)%type_prtl)
			case('refl')
				call RefBC_Prtl_Top(bc_face(4)%pos_prtl)
			case('remv')
				call RemovePrtl_Right( bc_pos_local(4, bc_face(4)%pos_prtl), qp,flvp,tagp,yp,xp,zp,up,vp,wp)	
		end select
		
		select case (bc_face(5)%type_prtl)
			case('remv')
				call RemovePrtl_Left( bc_pos_local(5, bc_face(5)%pos_prtl), qp,flvp,tagp,yp,xp,zp,up,vp,wp)
		end select
		
		select case (bc_face(6)%type_prtl)
			case('remv')
				call RemovePrtl_Right( bc_pos_local(6, bc_face(6)%pos_prtl), qp,flvp,tagp,yp,xp,zp,up,vp,wp)	
		end select
		
		
		
		if(curved_bc) call PrtlBC_curved_surf
	end subroutine EnforceBC_PostMovDep

	subroutine EnforceBC_Final
				
		
		select case (bc_face(1)%type_fld)
			case('cond')
				call CondBC_Fld_Left(bc_face(1)%pos_fld)
			case('open')
			    call OpenBC_Fld_Left(bc_face(1)%pos_fld)
			case('iflw')
				call Attenuate_EM_Flow(1)		
	    end select
		
		
		select case (bc_face(2)%type_fld)
			case('cond')
				call CondBC_Fld_Right(bc_face(2)%pos_fld)
			case('open')
			    call OpenBC_Fld_Right(bc_face(2)%pos_fld)
			case('iflw')
				call Attenuate_EM_Flow(2)	   
		end select
		
		
		select case (bc_face(3)%type_fld)
			case('cond')
				call CondBC_Fld_Bottom(bc_face(3)%pos_fld)
			case('iflw')
				call Attenuate_EM_Flow(3)		
	    end select
		
		
		select case (bc_face(4)%type_fld)
			case('cond')
				call CondBC_Fld_Top(bc_face(3)%pos_fld)
			case('iflw')
				call Attenuate_EM_Flow(4)		
	    end select
		
		
		if( attn_em_vacumm ) call Attenuate_EM_Vacuum
		
		if(bc_face(2)%type_prtl.eq.'open') call OpenBC_Prtl_Right(bc_face(2)%pos_prtl, 1.0_psn)
		if(bc_face(1)%type_prtl.eq.'open') call OpenBC_Prtl_Left(bc_face(1)%pos_prtl, 1.0_psn)
		
		if(curved_bc) call SetExteriorEfld_USC 
		
		call updateBC_Pos
	
	end subroutine EnforceBC_Final
	
	subroutine updateBC_Pos
		integer :: n
		do n=1,6
			bc_face(n)%pos_fld = bc_face(n)%pos_fld + bc_face(n)%speed*c
			bc_face(n)%pos_prtl = bc_face(n)%pos_prtl + bc_face(n)%speed*c
		end do 
	end subroutine updateBC_Pos
	
	real(psn) function bc_pos_local(side, pos)
	    integer :: side
		real(dbpsn) :: pos
		
		select case (side)
			case(1,2)
				bc_pos_local = pos - xborders(procxind) +3
			case(3,4)
				bc_pos_local = pos - yborders(procyind) +3
			case(5,6)
				bc_pos_local = pos - zborders(proczind) +3
		end select
	
	end function bc_pos_local
	
	subroutine CheckAttnBCPositioning(side)
		integer :: side
		logical :: bc_outside
	
		bc_outside = .false.
		
		select case (side)
			case(1)
			    if( bc_face(side)%pos_fld - bc_face(side)%attn_thickness .lt. xborders(0) ) bc_outside = .true. 
			case(2)
				if( bc_face(side)%pos_fld + bc_face(side)%attn_thickness .gt. xborders(nSubDomainsX) ) bc_outside = .true. 
			case(3)
				if( bc_face(side)%pos_fld - bc_face(side)%attn_thickness .lt. yborders(0) ) bc_outside = .true. 
			case(4)
				if( bc_face(side)%pos_fld + bc_face(side)%attn_thickness .gt. yborders(nSubDomainsY) ) bc_outside = .true. 
#ifndef twoD				
			case(5)
				if( bc_face(side)%pos_fld - bc_face(side)%attn_thickness .lt. zborders(0) ) bc_outside = .true. 
			case(6)
				if( bc_face(side)%pos_fld + bc_face(side)%attn_thickness .gt. zborders(nSubDomainsZ) ) bc_outside = .true. 
#endif				
		end select 
		
		if(bc_outside) then
			 print*,'error in setting boundary condition at side = ', side
			 STOP 'The domain where attenuation is applied (attn_thickness) extends beyond the simulation box. Change relative position of the boundary'
		end if
	
	    if(bc_face(side)%type_fld.eq.'iflw') then 
			if( abs(bc_face(side)%pos_fld-bc_face(side)%pos_prtl) .lt. curr_filter +1) then 
				
				if(proc.eq.0) then 
					print*,'Warning !!! Inflow BC issue: fld and prtl boundaries appear to be too close'
					print*,'If the filter is applied, current may be allowed to be filtered into the domain between prtl and fld boundary'
				end if 
			end if
		end if 
		
		
	end subroutine CheckAttnBCPositioning
	
	
	
!---------------------------------------------------------------------------------
!  BC  for EM Fld
!---------------------------------------------------------------------------------
subroutine OpenBC_Fld_Right(xind)
	real(dbpsn) :: xind
	integer :: xind_local,i
	xind_local=ceiling(xind-xborders(procxind)+3)
	if((xind_local.ge.1).and.(xind_local.le.mx-1)) then
! 		do i=xind_local+1,mx
! 			Ey(i,:,:)=Ey(xind_local,:,:)
! 			Ez(i,:,:)=Ez(xind_local,:,:)
! 			Bx(i,:,:)=Bx(xind_local,:,:)
! 		end do
!
! 		do i=xind_local,mx
! 			Ex(i,:,:)=Ex(xind_local-1,:,:)
! 			By(i,:,:)=0.0_psn
! 			Bz(i,:,:)=0.0_psn
! 		end do
		
		do i=xind_local+1,mx
			Ey(i,:,:)=0.5*Ey(i,:,:)
			Ez(i,:,:)=0.5*Ez(i,:,:)
			Bx(i,:,:)=Bx(xind_local,:,:)
		end do 

		do i=xind_local,mx
			Ex(i,:,:)=0.5*Ex(i,:,:)
			By(i,:,:)=0.5*By(i,:,:)
			Bz(i,:,:)=0.5*Bz(i,:,:)
		end do 		
	end if 
end subroutine OpenBC_Fld_Right

subroutine OpenBC_Fld_Left(xind)
	real(dbpsn) :: xind
	integer :: xind_local,i
	xind_local=ceiling(xind-xborders(procxind)+3)
	if((xind_local.ge.2).and.(xind_local.le.mx)) then
! 		do i=1,xind_local-1
! 			Ey(i,:,:)=Ey(xind_local,:,:)
! 			Ez(i,:,:)=Ez(xind_local,:,:)
! 			Bx(i,:,:)=Bx(xind_local,:,:)
!
! 			Ex(i,:,:)=Ex(xind_local,:,:)
! 			By(i,:,:)=By(xind_local,:,:)
! 			Bz(i,:,:)=Bz(xind_local,:,:)
!
! 			By(i,:,:)=0.0_psn
! 			Bz(i,:,:)=0.0_psn
! 		end do
		do i=1,xind_local-1
			Ey(i,:,:)=0.5*Ey(i,:,:)
			Ez(i,:,:)=0.5*Ez(i,:,:)
			Bx(i,:,:)=Bx(xind_local,:,:)
			
			Ex(i,:,:)=0.5*Ex(i,:,:)
			By(i,:,:)=0.5*By(i,:,:)
			Bz(i,:,:)=0.5*Bz(i,:,:)
		end do 
	end if 
end subroutine OpenBC_Fld_Left


!---------------------------------------------------------------------------------
!  BC  for Prtl
!---------------------------------------------------------------------------------

!---- Open Boundary condition
!The following version works only for two species(ion + electron)
subroutine OpenBC_Prtl_Right(xmax,dx)
	real(dbpsn) :: xmax
	real(psn)   ::  xmax_local, dx
	integer :: n,count
	
	!real(psn), dimension(my,mz):: den1,vel1,den2,vel2 
	real(psn), dimension(my,mz):: den,vel
	
	integer :: j,k
	real(psn) :: v0,gamma
	
! 	den1=1.0_psn
! 	den2=1.0_psn
! 	vel1=0.0_psn
! 	vel2=0.0_psn
	
	den=1.0_psn
	vel=0.0_psn 
	
	xmax_local=xmax-xborders(procxind)+3
	if(xmax_local.lt.mx-1) then 
		count=0
		do n=1,used_prtl_arr_size !remove old particles
			if((qp(n).ne.0).and.(xp(n).gt.xmax_local)) then 
				call DeletePrtl(n) 
			    count=count+1
			end if
		end do
		np=np-count
		
		
		do n=1,used_prtl_arr_size
			if((qp(n).ne.0).and.(xp(n).gt.xmax_local-dx).and.(xp(n).lt.xmax_local)) then
				call OpenBC_MeanPrtlSpeed(n,den,vel)
				!if(flvp(n).eq.1) call OpenBC_MeanPrtlSpeed(n,den1,vel1)
				!if(flvp(n).eq.2) call OpenBC_MeanPrtlSpeed(n,den2,vel2)
			end if 
		end do 
		
		vel=vel/den
		
		do n=1,used_prtl_arr_size
			if((qp(n).ne.0).and.(xp(n).gt.xmax_local-dx).and.(xp(n).lt.xmax_local)) then
				j=yp(n)
#ifdef twoD
                k=1
#else 
               	k=zp(n)
#endif	
! 				if(flvp(n).eq.1) v0=vel1(j,k)
! 				if(flvp(n).eq.2) v0=vel2(j,k)
                v0 = vel(j,k)
				if(v0.gt.0) v0=0.0
				gamma=sqrt(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n))
                call InsertNewPrtl(2*xmax_local-xp(n),yp(n),zp(n),up(n)-gamma*v0,vp(n),wp(n),qp(n),tagp(n),flvp(n),var1p(n))
			end if
		end do
	end if
	
end subroutine OpenBC_Prtl_Right

subroutine OpenBC_Prtl_Left(xmin,dx)
	real(dbpsn) :: xmin
	real(psn)   :: xmin_local,dx
	integer :: n,count
	!real(psn), dimension(my,mz):: den1,vel1,den2,vel2
	real(psn), dimension(my,mz):: den,vel
	integer :: j,k
	real(psn) :: v0,gamma
	
! 	den1=1.0
! 	den2=1.0
! 	vel1=0.0
! 	vel2=0.0

	den=1.0
	vel=0.0
	
	xmin_local=xmin-xborders(procxind)+3
	if(xmin_local.gt.1) then 
		count=0
		do n=1,used_prtl_arr_size !remove old particles
			if((qp(n).ne.0).and.(xp(n).lt.xmin_local)) then 
				call DeletePrtl(n)
	            count=count+1	
			end if 
		end do
		np=np-count
		
		
		do n=1,used_prtl_arr_size
			if((qp(n).ne.0).and.(xp(n).lt.xmin_local+dx).and.(xp(n).gt.xmin_local)) then
				call OpenBC_MeanPrtlSpeed(n,den,vel)
				!if(flvp(n).eq.1) call OpenBC_MeanPrtlSpeed(n,den1,vel1)
				!if(flvp(n).eq.2) call OpenBC_MeanPrtlSpeed(n,den2,vel2)
			end if 
		end do 
		vel=vel/den
		
		
		do n=1,used_prtl_arr_size
			if((qp(n).ne.0).and.(xp(n).lt.xmin_local+dx).and.(xp(n).gt.xmin_local)) then
				j=yp(n)
#ifdef twoD
                k=1
#else 
               	k=zp(n)
#endif				
				v0=vel(j,k)
! 				if(flvp(n).eq.1) v0=vel1(j,k)
! 				if(flvp(n).eq.2) v0=vel2(j,k)
				if(v0.lt.0) v0=0.0
				gamma=sqrt(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n))
				call InsertNewPrtl(2*xmin_local-xp(n),yp(n),zp(n),up(n)-gamma*v0,vp(n),wp(n),qp(n),tagp(n),flvp(n),var1p(n))
			end if
		end do
		
		
	end if

end subroutine OpenBC_Prtl_Left

!cell-by-cell mean speed of particles which are near the boundary layer and within the physical domain 
subroutine OpenBC_MeanPrtlSpeed(n,den,vel)
	integer :: n
	integer :: j,k
	real(psn), dimension(my,mz) :: den,vel
	real(psn) :: vx,invg
	
	if(flvp(n).le.2) then 
		j=yp(n)
#ifdef twoD
        k=1
#else 
        k=zp(n)
#endif		
        invg=1.0_psn/sqrt(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n))
        vx=up(n)*invg
		den(j,k)=den(j,k)+1.0
		vel(j,k)=vel(j,k)+vx			
	end if
end subroutine OpenBC_MeanPrtlSpeed
	

end module bc