module bc_reflect
	use parameters
	use vars
contains 

!---------------------------------------------------------------------------------
!
!                       Conducting/Reflecting BC for EM Field
!
!---------------------------------------------------------------------------------

subroutine CondBC_Fld_Top(yind)
	real(dbpsn) :: yind
	integer :: yind_local
	yind_local=ceiling(yind-yborders(procyind)+3)
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
	yind_local=floor(yind-yborders(procyind)+3)
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
	xind_local=ceiling(xind-xborders(procxind)+3)
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
	xind_local=floor(xind-xborders(procxind)+3) 
	if(xind_local.lt.1) return
    xind_local = min(mx,xind_local)

	Ey(1:xind_local,:,:)=0.0_psn
	Ez(1:xind_local,:,:)=0.0_psn

#ifdef gpu
    call CondBC_Fld_Left_GPU(xind_local)
#endif		
end subroutine CondBC_Fld_Left


!---------------------------------------------------------------------------------
!
!                       Reflecting BC for particles
!
!---------------------------------------------------------------------------------


subroutine RefBC_Prtl_Top(ymax)
	real(dbpsn) :: ymax
	real(psn)  :: ymax_local
	integer :: n,PrtlWall_local

	ymax_local=ymax-yborders(procyind)+3
	PrtlWall_local=ymax_local

#ifdef gpu
    call RefBC_Prtl_Top_GPU(ymax_local,PrtlWall_local)
	return
#endif	

	if(ymax_local.lt.my-1) call RefPrtlRight(flvp,yp,vp,ymax_local)

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

	ymin_local=ymin-yborders(procyind)+3
    PrtlWall_local=ymin_local

#ifdef gpu
    call RefBC_Prtl_Bottom_GPU(ymin_local,PrtlWall_local)
	return
#endif	
	
	if(ymin_local.gt.1) call RefPrtlLeft(flvp,yp,vp,ymin_local)


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

	xmax_local=xmax-xborders(procxind)+3! - 0.5 ! 0.5: refl. boundary is in the middle of a cell; normal flux excatly zero; assuming static boundary 

	if(xmax_local.lt.mx-1) call RefPrtlRight(flvp,xp,up,xmax_local)
	
	call ReflCurrentRight(xmax)

end subroutine RefBC_Prtl_Right

subroutine RefBC_Prtl_Left(xmin)
	real(dbpsn) :: xmin
	real(psn)   :: xmin_local
	
	xmin_local=xmin-xborders(procxind)+3
	
	if(xmin_local.gt.1) call RefPrtlLeft(flvp,xp,up,xmin_local)
	
	call ReflCurrentLeft(xmin)
	
end subroutine RefBC_Prtl_Left


!---------------------------------------------------------------------------------
!
!                       Specular reflection of particles
!
!---------------------------------------------------------------------------------

subroutine RefPrtlLeft(flv,x,u,xmin)
	integer, dimension(:) :: flv
	real(psn), dimension(:) :: x,u
	real(psn) :: xmin
	integer :: n
	
	do n=1,used_prtl_arr_size
		if((flv(n).ne.0).and.(x(n).lt.xmin)) then 
			u(n)=-u(n)
			x(n)=2.0_psn*xmin-x(n)
		end if
	end do

end subroutine RefPrtlLeft

subroutine RefPrtlRight(flv,x,u,xmax)
	integer, dimension(:) :: flv
	real(psn), dimension(:) :: x,u
	real(psn) :: xmax
	integer :: n
	
	do n=1,used_prtl_arr_size
		if((flv(n).ne.0).and.(x(n).gt.xmax)) then 
			u(n)=-u(n)
			x(n)=2.0_psn*xmax-x(n)
		end if
	end do

end subroutine RefPrtlRight

!---------------------------------------------------------------------------------
!
!         Reflection of current outside the boundaries
!
!---------------------------------------------------------------------------------

subroutine ReflCurrentLeft(pos)
	 real(dbpsn) :: pos
	 integer :: ind_local
 
	 ind_local = ceiling( pos -xborders(procxind)+3 )
	 if(ind_local.ge.2.and.ind_local.le.mx) then 
		 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local-1,:,:)
		 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local-1,:,:)
		 Jy(ind_local-1,:,:) = 0 
		 Jz(ind_local-1,:,:) = 0
	 end if 

	 ind_local = ceiling(pos-0.5_psn -xborders(procxind)+3)
	 if(ind_local.ge.2.and.ind_local.le.mx) then 
		 Jx(ind_local,:,:) = Jx(ind_local,:,:) - Jx(ind_local-1,:,:)
		 Jx(ind_local-1,:,:) = 0
	 end if 
end subroutine ReflCurrentLeft

subroutine ReflCurrentRight(pos)
	 real(dbpsn) :: pos
	 integer :: ind_local
 
	 ind_local = floor(pos-xborders(procxind)+3)
	 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
		 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local+1,:,:)
		 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local+1,:,:)
		 Jy(ind_local+1,:,:) = 0 
		 Jz(ind_local+1,:,:) = 0
	 end if 

 
	 ind_local = floor(bc_face(2)%pos_prtl-0.5_psn -xborders(procxind)+3)
	 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
		 Jx(ind_local,:,:) = Jx(ind_local,:,:) - Jx(ind_local+1,:,:)
		 Jx(ind_local+1,:,:) = 0
	 end if 
end subroutine ReflCurrentRight

	
end module bc_reflect