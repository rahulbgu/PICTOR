module bc
    use parameters
	use vars
	use memory
	use communication
	use injector
	use bc_curved_surf
#ifdef gpu
    use bc_gpu 
#endif	
	implicit none 
contains 	
!---------------------------------------------------------------------------------
!  Subroutines to implement BCs from the setup file 
!---------------------------------------------------------------------------------
    subroutine SetBC_Fld(Axis,Position,Type,Side)
		 character (len=*) :: Axis,Type,Side
		 real(psn)  :: Position
		 if(Axis.eq.'x') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Xmin_Fld=Position
				 BC_Xmin_Fld_Type=Type
				 call PeriodicX_Off
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Xmax_Fld=Position
				 BC_Xmax_Fld_Type=Type
				 call PeriodicX_Off
			 end if 
		 end if 
		 if(Axis.eq.'y') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Ymin_Fld=Position
				 BC_Ymin_Fld_Type=Type
				 call PeriodicY_Off
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Ymax_Fld=Position
				 BC_Ymax_Fld_Type=Type
				 call PeriodicY_Off
			 end if 
		 end if 
		 if(Axis.eq.'z') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Zmin_Fld=Position
				 BC_Zmin_Fld_Type=Type
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Zmax_Fld=Position
				 BC_Zmax_Fld_Type=Type
			 end if 
		 end if 
	end subroutine SetBC_Fld	 		

    subroutine SetBC_Prtl(Axis,Position,Type,Side)
		 character (len=*) :: Axis,Type,Side
		 real(psn)  :: Position, dx
		 if(Axis.eq.'x') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Xmin_Prtl=Position
				 BC_Xmin_Prtl_Type=Type
				 call PeriodicX_Off
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Xmax_Prtl=Position
				 BC_Xmax_Prtl_Type=Type
				 call PeriodicX_Off
			 end if 
		 end if 
		 if(Axis.eq.'y') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Ymin_Prtl=Position
				 BC_Ymin_Prtl_Type=Type
				 call PeriodicY_Off
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Ymax_Prtl=Position
				 BC_Ymax_Prtl_Type=Type
				 call PeriodicY_Off
			 end if 
		 end if 
		 if(Axis.eq.'z') then
			 if(Side.eq.'lower') then 
				 if(.not.restart) BC_Zmin_Prtl=Position
				 BC_Zmin_Prtl_Type=Type
			 end if 
			 if(Side.eq.'upper') then 
				 if(.not.restart) BC_Zmax_Prtl=Position
				 BC_Zmax_Prtl_Type=Type
			 end if 
		 end if 
	end subroutine SetBC_Prtl
	
	subroutine SetBC_InflowPrtl(Axis,Position,Side,InjSpeed,Flvr1,Flvr2,Density,Temperature1,Temperature2,SpeedDist1,SpeedDist2,DriftVelocity,Vmax1,Vmax2)
		character (len=*) :: Axis,Side
		integer :: Flvr1,Flvr2
		integer :: side_ind
		real(psn)  :: Position
		real(psn), optional :: Vmax1 !the maximum particle speed in the plasma frame for Flvr1, default is c 
		real(psn), optional :: Vmax2 ! max. speed for Flvr2; Vmax1 and Vmax2 are only used when speed dist is provided
		real(psn), optional :: InjSpeed
		real(psn), external, optional :: Temperature1
		real(psn), external, optional :: Temperature2
		real(psn), external, optional :: SpeedDist1
		real(psn), external, optional :: SpeedDist2
		real(psn), external  :: Density 
		real(psn) :: vmax
		interface 
			subroutine DriftVelocity(x,y,z,fx,fy,fz)
				import :: psn
				real(psn) :: x,y,z,fx,fy,fz
			end subroutine
	    end interface 
		
		inflowBC = .true.
		
		call SetBC_Prtl(Axis,Position,'iflw',Side)
		
		if(Axis.eq.'x') then
			if(Side.eq.'lower') side_ind = 1
			if(Side.eq.'upper') side_ind = 2   
			call PeriodicX_Off
		end if  
		if(Axis.eq.'y') then
			if(Side.eq.'lower') side_ind = 3
			if(Side.eq.'upper') side_ind = 4  
			call PeriodicY_Off  
		end if 
		
		
		if(.not.allocated(InjPrtl)) then 
			 allocate(InjPrtl(8)) ! change here or allocate manually if a larger size is required 
			 nInjPrtl=0
		end if 
		
		nInjPrtl = nInjPrtl + 1
		
		InjPrtl(nInjPrtl)%side  = side_ind
		InjPrtl(nInjPrtl)%flvr1 = Flvr1 
		InjPrtl(nInjPrtl)%flvr2 = Flvr2 
		
		InjPrtl(nInjPrtl)%Den => Density
		InjPrtl(nInjPrtl)%Drift => DriftVelocity
		
		if(present(Temperature1)) then 
			InjPrtl(nInjPrtl)%dist_type1 = 1
			InjPrtl(nInjPrtl)%Temp1 => Temperature1
		end if 
		
		if(present(Temperature2)) then 
			InjPrtl(nInjPrtl)%dist_type2 = 1
			InjPrtl(nInjPrtl)%Temp2 => Temperature2
		end if 
		
		if(present(SpeedDist1)) then 
			InjPrtl(nInjPrtl)%dist_type1 = 2
			InjPrtl(nInjPrtl)%SpeedDist1 => SpeedDist1
			vmax = 1.0 ! v/c 
			if(present(Vmax1)) vmax = Vmax1
			call InitPDFTable(InjPrtl(nInjPrtl)%TableSize,InjPrtl(nInjPrtl)%Table1,InjPrtl(nInjPrtl)%PDF_Table1,SpeedDist1,vmax)
		end if 
		
		if(present(SpeedDist2)) then 
			InjPrtl(nInjPrtl)%dist_type2 = 2
			InjPrtl(nInjPrtl)%SpeedDist2 => SpeedDist2
			vmax = 1.0  
			if(present(Vmax2)) vmax = Vmax2
			call InitPDFTable(InjPrtl(nInjPrtl)%TableSize,InjPrtl(nInjPrtl)%Table2,InjPrtl(nInjPrtl)%PDF_Table2,SpeedDist2,vmax)
		end if 
			

		
		InflowFld(side_ind)%Drift => DriftVelocity
		
		inflowBC_speed(side_ind) = 0
		if(present(InjSpeed)) inflowBC_speed(side_ind) = InjSpeed

	end subroutine SetBC_InflowPrtl
	
	subroutine SetBC_InflowFld(Axis,Position,Side,InjSpeed,DriftVelocity,MagFld,DriftVmax,Attenuate)
		character (len=*) :: Axis,Side
		integer :: side_ind
		real(psn)  :: Position
		real(psn), optional :: DriftVmax !the maximum drift speed, helps improve the efficiency of prtl injection process
		real(psn), optional :: Attenuate
		real(psn), optional :: InjSpeed
		
		interface 
			subroutine MagFld(x,y,z,fx,fy,fz)
				import :: psn
				real(psn) :: x,y,z,fx,fy,fz
			end subroutine
			subroutine DriftVelocity(x,y,z,fx,fy,fz)
				import :: psn
				real(psn) :: x,y,z,fx,fy,fz
			end subroutine
	    end interface
		
		inflowBC = .true.
		
		call SetBC_Fld(Axis,Position,'iflw',Side)
		
		if(Axis.eq.'x') then
			if(Side.eq.'lower') side_ind = 1
			if(Side.eq.'upper') side_ind = 2   
			call PeriodicX_Off
		end if  
		if(Axis.eq.'y') then
			if(Side.eq.'lower') side_ind = 3
			if(Side.eq.'upper') side_ind = 4  
			call PeriodicY_Off  
		end if 
		
		InflowFld(side_ind)%MagFld => MagFld
		InflowFld(side_ind)%Drift => DriftVelocity
		
		InflowFld(side_ind)%vmax = 1.0
		if(present(DriftVmax)) InflowFld(side_ind)%vmax = DriftVmax
		
		inflowBC_speed(side_ind) = 0
		if(present(InjSpeed)) inflowBC_speed(side_ind) = InjSpeed
		
		InflowFld(side_ind)%Attenuate = 0 
		if(present(Attenuate)) InflowFld(side_ind)%Attenuate = Attenuate
	end subroutine SetBC_InflowFld
	
	subroutine PeriodicX_Off
	    if(procxind(proc).eq.0) lproc=MPI_PROC_NULL
	    if(procxind(proc).eq.nSubDomainsX-1) rproc=MPI_PROC_NULL
	end subroutine PeriodicX_Off	
	subroutine PeriodicY_Off
	    if(procyind(proc).eq.0) bproc=MPI_PROC_NULL
	    if(procyind(proc).eq.nSubDomainsY-1) tproc=MPI_PROC_NULL
	end subroutine PeriodicY_Off	
!---------------------------------------------------------------------------------
!  Essential rotuines to enforece BCs from the main loop 
!---------------------------------------------------------------------------------	
	subroutine EnforceBC_PostMovDep

		if(inflowBC) call InflowBC_Prtl
				
		if(BC_Ymax_Prtl.gt.0.and.BC_Ymax_Prtl_Type.eq.'refl') call RefBC_Prtl_Top(BC_Ymax_Prtl)
		if(BC_Ymin_Prtl.gt.0.and.BC_Ymin_Prtl_Type.eq.'refl') call RefBC_Prtl_Bottom(BC_Ymin_Prtl)	
		if(BC_Xmax_Prtl_Type.eq.'refl') call RefBC_Prtl_Right(BC_Xmax_Prtl)
		if(BC_Xmin_Prtl_Type.eq.'refl') call RefBC_Prtl_Left(BC_Xmin_Prtl)	
		
		if(curved_bc) call PrtlBC_curved_surf
	end subroutine EnforceBC_PostMovDep

	subroutine EnforceBC_Final
				
		
		select case (BC_Xmin_Fld_Type)
			case('cond')
				call CondBC_Fld_Left(BC_Xmin_Fld)
			case('open')
			    call OpenBC_Fld_Left(BC_Xmin_Fld)
			case('iflw')
				call InflowBC_Fld(InflowFld(1),1)		
	    end select
		
		
		select case (BC_Xmax_Fld_Type)
			case('cond')
				call CondBC_Fld_Right(BC_Xmax_Fld)
			case('open')
			     call OpenBC_Fld_Right(BC_Xmax_Fld)
			case('iflw')
				call  InflowBC_Fld(InflowFld(2),2)	    
		end select
		
		
		select case (BC_Ymin_Fld_Type)
			case('cond')
				call CondBC_Fld_Bottom(BC_Ymin_Fld)
			case('iflw')
				call InflowBC_Fld(InflowFld(3),3)		
	    end select
		
		
		select case (BC_Ymax_Fld_Type)
			case('cond')
				call CondBC_Fld_Top(BC_Ymax_Fld)
			case('iflw')
				call InflowBC_Fld(InflowFld(4),4)		
	    end select
		
		
		
		
		
		if(BC_Xmax_Prtl.gt.0.and.BC_Xmax_Prtl_Type.eq.'open') call OpenBC_Prtl_Right(BC_Xmax_Prtl, 1.0_psn)
		if(BC_Xmin_Prtl.gt.0.and.BC_Xmin_Prtl_Type.eq.'open') call OpenBC_Prtl_Left(BC_Xmin_Prtl, 1.0_psn)
		
		if(curved_bc) call SetExteriorEfld_USC 
	
	end subroutine EnforceBC_Final
	

!---------------------------------------------------------------------------------
!
!                       Conducting BC
!
!---------------------------------------------------------------------------------

subroutine CondBC_Fld_Top(yind)
	real(dbpsn) :: yind
	integer :: yind_local
	yind_local=ceiling(yind-yborders(procyind(proc))+3)
	if(yind_local.gt.my) return
    yind_local = max(1,yind_local)
	
	Ex(:,yind_local:my,:)=0.0_psn
	Ez(:,yind_local:my,:)=0.0_psn
	
#ifdef gpu
    call CondBC_Fld_Top_GPU(yind_local)
#endif		
end subroutine CondBC_Fld_Top

subroutine CondBC_Fld_Bottom(yind)
	real(dbpsn) :: yind
	integer :: yind_local
	yind_local=floor(yind-yborders(procyind(proc))+3)
	if(yind_local.lt.1) return
    yind_local = min(my,yind_local)

	Ex(:,1:yind_local,:)=0.0_psn
	Ez(:,1:yind_local,:)=0.0_psn

#ifdef gpu
    call CondBC_Fld_Bottom_GPU(yind_local)
#endif		
end subroutine CondBC_Fld_Bottom

subroutine CondBC_Fld_Right(xind)
	real(dbpsn) :: xind
	integer :: xind_local
	xind_local=ceiling(xind-xborders(procxind(proc))+3)
	if(xind_local.gt.mx) return
    xind_local = max(1,xind_local)
	
	Ey(xind_local:mx,:,:)=0.0_psn
	Ez(xind_local:mx,:,:)=0.0_psn
	
#ifdef gpu
    !call CondBC_Fld_Right_GPU(xind_local) ! not implemented yet
#endif		
end subroutine CondBC_Fld_Right

subroutine CondBC_Fld_Left(xind)
	real(dbpsn) :: xind
	integer :: xind_local
	xind_local=floor(xind-xborders(procxind(proc))+3) 
	if(xind_local.lt.1) return
    xind_local = min(mx,xind_local)

	Ey(1:xind_local,:,:)=0.0_psn
	Ez(1:xind_local,:,:)=0.0_psn

#ifdef gpu
    call CondBC_Fld_Left_GPU(xind_local)
#endif		
end subroutine CondBC_Fld_Left

subroutine RefBC_Prtl_Top(ymax)
	real(dbpsn) :: ymax
	real(psn)  :: ymax_local
	integer :: n,PrtlWall_local
	
	ymax_local=ymax-yborders(procyind(proc))+3
	PrtlWall_local=ymax_local

#ifdef gpu
    call RefBC_Prtl_Top_GPU(ymax_local,PrtlWall_local)
	return
#endif	
	
	if(ymax_local.lt.my-1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).gt.ymax_local)) then 
				 vp(n)=-vp(n)
                 yp(n)=ymax_local-(yp(n)-ymax_local)
			end if 
     	end do
	end if
	

	if(PrtlWall_local.ge.2.and.PrtlWall_local.le.my-2) then
		Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local+1,:)
	 	Jy(:,PrtlWall_local-1,:)=Jy(:,PrtlWall_local-1,:)-Jy(:,PrtlWall_local,:)
	 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local+1,:)
		Jx(:,PrtlWall_local+1,:)=0.0_psn
		Jy(:,PrtlWall_local,:)=0.0_psn
		Jz(:,PrtlWall_local+1,:)=0.0_psn
    end if
	
end subroutine RefBC_Prtl_Top
subroutine RefBC_Prtl_Bottom(ymin)
	real(dbpsn) :: ymin
	real(psn)   :: ymin_local
	integer :: n,PrtlWall_local
	
	ymin_local=ymin-yborders(procyind(proc))+3
    PrtlWall_local=ymin_local
	
#ifdef gpu
    call RefBC_Prtl_Bottom_GPU(ymin_local,PrtlWall_local)
	return
#endif	
		
	if(ymin_local.gt.1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(yp(n).lt.ymin_local)) then 
				vp(n)=-vp(n)
				yp(n)=ymin_local+(ymin_local-yp(n))
			end if 
		end do 
	end if
	

	if(PrtlWall_local.ge.3.and.PrtlWall_local.le.my-2) then
		Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local-1,:)
	 	Jy(:,PrtlWall_local,:)=Jy(:,PrtlWall_local,:)-Jy(:,PrtlWall_local-1,:)
	 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local-1,:)
		Jx(:,PrtlWall_local-1,:)=0.0_psn
		Jy(:,PrtlWall_local-1,:)=0.0_psn
		Jz(:,PrtlWall_local-1,:)=0.0_psn
	end if
end subroutine RefBC_Prtl_Bottom

subroutine RefBC_Prtl_Right(xmax)
	real(dbpsn) :: xmax
	real(psn)   :: xmax_local
	integer :: n,PrtlWall_local
	
	xmax_local=xmax-xborders(procxind(proc))+3
	PrtlWall_local=xmax_local

#ifdef gpu
    !call RefBC_Prtl_Right_GPU(ymax_local,PrtlWall_local) ! not implemented 
	!return
#endif	
	
	if(xmax_local.lt.mx-1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(xp(n).gt.xmax_local)) then 
				 up(n)=-up(n)
                 xp(n)=xmax_local-(xp(n)-xmax_local)
			end if 
     	end do
	end if
	

	if(PrtlWall_local.ge.2.and.PrtlWall_local.le.mx-2) then
		Jy(PrtlWall_local,:,:)=Jy(PrtlWall_local,:,:)+Jy(PrtlWall_local+1,:,:)
	 	Jx(PrtlWall_local-1,:,:)=Jx(PrtlWall_local-1,:,:)-Jx(PrtlWall_local,:,:)
	 	Jz(PrtlWall_local,:,:)=Jz(PrtlWall_local,:,:)+Jz(PrtlWall_local+1,:,:)
		Jy(PrtlWall_local+1,:,:)=0.0_psn
		Jx(PrtlWall_local,:,:)=0.0_psn
		Jz(PrtlWall_local+1,:,:)=0.0_psn
    end if
	
end subroutine RefBC_Prtl_Right

subroutine RefBC_Prtl_Left(xmin)
	real(dbpsn) :: xmin
	real(psn)   :: xmin_local
	integer :: n,PrtlWall_local
	
	! in the new scheme prtl reflecting boundary is in the middle of a cell; this is to simply incorporate zero normal flux 
	xmin_local=xmin-xborders(procxind(proc))+3 + 0.5  
    PrtlWall_local=xmin_local
	
#ifdef gpu
    call RefBC_Prtl_Left_GPU(xmin_local,PrtlWall_local) ! not implemented
	return
#endif	
		
	if(xmin_local.gt.1) then 
		do n=1,used_prtl_arr_size !particles
			if((qp(n).ne.0).and.(xp(n).lt.xmin_local)) then 
				up(n)=-up(n)
				xp(n)=xmin_local+(xmin_local-xp(n))
			end if 
		end do 
	end if
	

	if(PrtlWall_local.ge.3.and.PrtlWall_local.le.mx-2) then
		
		Jy(PrtlWall_local+1,:,:)=Jy(PrtlWall_local+1,:,:)+Jy(PrtlWall_local,:,:)
	 	Jx(PrtlWall_local,:,:)=0
	 	Jz(PrtlWall_local+1,:,:)=Jz(PrtlWall_local+1,:,:)+Jz(PrtlWall_local,:,:)
		Jy(PrtlWall_local,:,:)=0.0_psn
		Jz(PrtlWall_local,:,:)=0.0_psn
		
! 		Jy(PrtlWall_local,:,:)=Jy(PrtlWall_local,:,:)+Jy(PrtlWall_local-1,:,:)
! 	 	Jx(PrtlWall_local,:,:)=Jx(PrtlWall_local,:,:)-Jx(PrtlWall_local-1,:,:)
! 	 	Jz(PrtlWall_local,:,:)=Jz(PrtlWall_local,:,:)+Jz(PrtlWall_local-1,:,:)
! 		Jy(PrtlWall_local-1,:,:)=0.0_psn
! 		Jx(PrtlWall_local-1,:,:)=0.0_psn
! 		Jz(PrtlWall_local-1,:,:)=0.0_psn
	end if
end subroutine RefBC_Prtl_Left
	
	
	
!---------------------------------------------------------------------------------
!  BC  for EM Fld
!---------------------------------------------------------------------------------
subroutine OpenBC_Fld_Right(xind)
	real(dbpsn) :: xind
	integer :: xind_local,i
	xind_local=ceiling(xind-xborders(procxind(proc))+3)
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
	xind_local=ceiling(xind-xborders(procxind(proc))+3)
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
	
	xmax_local=xmax-xborders(procxind(proc))+3
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
	
	xmin_local=xmin-xborders(procxind(proc))+3
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