module cyl_bc_cpu
    use parameters
    use vars
	use cyl_vars
	implicit none 
contains 
!---------------------------------------------------------------------------------
!  BC  for EM Fld
!---------------------------------------------------------------------------------

	subroutine ConductingOuterBC_Fld(rmax)
		real(psn) :: rmax
		integer   :: rmax_local
		integer  :: i
		rmax_local=rmax-rborders(procxind(proc))+3
		rmax_local=max(rmax_local,1)
		if(rmax_local.le.mx) then 
			Ey(rmax_local:mx,:,:)=0.0_psn
			Ez(rmax_local:mx,:,:)=0.0_psn
		end if 
	end subroutine ConductingOuterBC_Fld

	subroutine ConductingInnerBC_Fld(rmin) !Note: place the conducting wall a few cells away from the particle wall
		real(psn) :: rmin
		integer   :: rmin_local
		rmin_local=rmin-rborders(procxind(proc))+3
		rmin_local=min(rmin_local,mx)
		if(rmin_local.ge.1) then 
			Ey(1:rmin_local,:,:)=0.0_psn
			Ez(1:rmin_local,:,:)=0.0_psn
		end if 
	end subroutine ConductingInnerBC_Fld

! 	subroutine AxisCurrentBC
!   		 integer :: n,ind,j,k
!   		 real(psn) :: phi,wt
! 		 if(procxind(proc).ne.0) return
!
! ! 		wt=0
! ! 		do n=0,ny-1
! ! 			phi= 0.5*dtheta + n*dtheta
! ! 			wt=wt+abs(sin(phi))
! ! 		end do
! ! 		wt=1.0_psn/wt
!
! #ifdef twoD
!         do k=1,1
! #else
!         do k=1,mz
! #endif
! !            do j=1,my
! ! 				do n=0,ny-1
! ! 					phi= 0.5*dtheta + n*dtheta
! ! 					ind=j+n
! ! 					if(ind.gt.my-3) ind= ind - ny
! ! 					Jy(4,j,k)=Jy(4,j,k)-Jx(3,ind,k)*sin(phi)*wt
! ! 				end do
! ! 		   end do
!
! 		   do j=1,my
! 			ind=j+ny/2
! 			if(ind.gt.my-3) ind= ind - ny
! 			Jx(4,j,k)=Jx(4,j,k)+0.5_psn*Jx(3,ind,k)
! 			Jx(4,ind,k)=Jx(4,ind,k)-0.5_psn*Jx(3,ind,k)
! 			!Jy(4,ind,k)=Jy(4,ind,k)-Jy(3,j,k)
! 			Jy(4,ind,k)=Jy(4,ind,k)+Jy(3,j,k)
! 			Jz(4,ind,k)=Jz(4,ind,k)+Jz(3,j,k)
! 		   end do
!         end do
!
! 		Jx(3,:,:)=0.0_psn
! 		Jy(3,:,:)=0.0_psn
! 		Jz(3,:,:)=0.0_psn
! 	end subroutine AxisCurrentBC

    subroutine AxisCurrentBC
		integer :: j,k,jp,jm,jpp
		if(procxind(proc).ne.0) return

! #ifdef twoD
!          	do k=1,1
! #else
!          	do k=1,mz
! #endif
!             	do j=1,my
! 					jp=j + ny/4
! 					jm=j - ny/4
! 					jpp = j + ny/2
! 					if(jp.gt.my-3) jp= jp - ny
! 					if(jpp.gt.my-3) jpp= jpp - ny
! 					if(jm.lt.3) jm = jm + ny
!
! 					Jy(3,jp,k) = Jy(3,jp,k) - 0.5_psn*Jx(3,j,k)
! 					Jy(3,jm,k) = Jy(3,jm,k) + 0.5_psn*Jx(3,j,k)
!
! 					Jy(3,jpp,k) = Jy(3,jpp,k) - Jy(2,jpp,k)
! 					Jz(3,jpp,k) = Jz(3,jpp,k) + Jz(2,jpp,k)
!
! 				end do
! 			end do
           
		    Jy(4,:,:)=Jy(4,:,:)-Jy(3,:,:)
            Jz(4,:,:)=Jz(4,:,:)-Jz(3,:,:)

			Jx(3,:,:) = 0.0_psn
			Jy(2,:,:) = 0.0_psn
			Jz(2,:,:) = 0.0_psn


	end subroutine AxisCurrentBC 
	
	 
	 subroutine UpdateFldAxis
		 integer :: j,k,jp,jm,jpp
! 		 if((inc_axis.eqv..true.).and.(procxind(proc).eq.0)) then
! 			 Ex(3,:,:)=0.0_psn
! 			 Ex(2,:,:)=-Ex(4,:,:)
!
! 			 Ey(3,:,:)=-Ey(4,:,:)
! 			 Ez(3,:,:)=-Ez(4,:,:)
!
! 			 Bx(3,:,:)=-Bx(4,:,:)
! 			 By(3,:,:)=0.0_psn
! 			 By(2,:,:)=-By(4,:,:)
! 			 !Bz(2,:,:)=2*Bz(3,:,:)-Bz(4,:,:)
! 			 Bz(2,:,:)=Bz(3,:,:)
! 	     end if
	 end subroutine UpdateFldAxis
	 

!---------------------------------------------------------------------------------
!  BC  for Prtl
!---------------------------------------------------------------------------------
	subroutine RefInnerBC_Prtl(rmin) 
		real(psn) :: rmin,rmin_local
		integer   :: ind,n
	    rmin_local=rmin-rborders(procxind(proc))+3
		ind=rmin_local
		if(rmin_local.lt.0) return
	
		do n=1,used_prtl_arr_size
			if((xp(n).lt.rmin_local).and.(qp(n).ne.0)) then
				xp(n)=rmin_local+(rmin_local-xp(n))
				up(n)=-up(n)
			end if
		end do 
		do n=1,used_test_prtl_arr_size
			if((xtp(n).lt.rmin_local).and.(qtp(n).ne.0)) then
				xtp(n)=rmin_local+(rmin_local-xtp(n))
				utp(n)=-utp(n)
			end if
		end do 
		if(ind.ge.3.and.ind.le.mx-2) then
			!To ensure that the current is deposited on right place for reflected particles
			!Jx(ind,:,:)=Jx(ind,:,:)-Jx(ind-1,:,:)
			! factor of f? 
		 	Jy(ind+1,:,:)=Jy(ind+1,:,:)+Jy(ind,:,:)
		 	Jz(ind+1,:,:)=Jz(ind+1,:,:)+Jz(ind,:,:)
			Jx(ind,:,:)=0.0_psn
			Jy(ind,:,:)=0.0_psn
			Jz(ind,:,:)=0.0_psn
		end if
	
	end subroutine RefInnerBC_Prtl

	subroutine RefOuterBC_Prtl(rmax) !rmax must be a whole number
		real(psn) :: rmax,rmax_local,f
		integer   :: ind,n
	    !rmax_local=rmax-rborders(procxind(proc))+xmin !prtl boundary lies in the middle of the cell
		rmax_local=BC_Rmax_Prtl-rborders(procxind(proc))+3
		ind=rmax_local
		if(rmax_local.gt.mx) return
	
		do n=1,used_prtl_arr_size
			if((xp(n).gt.rmax_local).and.(qp(n).ne.0)) then
				xp(n)=rmax_local-(xp(n)-rmax_local)
				up(n)=-up(n)
			end if
		end do 
		do n=1,used_test_prtl_arr_size
			if((xtp(n).gt.rmax_local).and.(qtp(n).ne.0)) then
				xtp(n)=rmax_local-(xtp(n)-rmax_local)
				utp(n)=-utp(n)
			end if
		end do 
		call CurrentBC_Outer

	end subroutine RefOuterBC_Prtl
	
	subroutine CurrentBC_Outer
		integer :: ind
		real(psn) :: f
		ind= BC_Rmax_Prtl-rborders(procxind(proc))+xmin
		!ind= BC_Rmax_Prtl-xborders(procxind(proc))+xmin
		if(ind.ge.2.and.ind.le.mx-2) then
			f=(BC_Rmax_Prtl+0.5_psn)/(BC_Rmax_Prtl-0.5_psn) 
			!To ensure that the current is deposited on right place for reflected particles
			!Jx(ind-1,:,:)=Jx(ind-1,:,:)-Jx(ind,:,:)
		 	Jy(ind,:,:)=Jy(ind,:,:)+Jy(ind+1,:,:)!*f
		 	Jz(ind,:,:)=Jz(ind,:,:)+Jz(ind+1,:,:)!*f 
			Jx(ind,:,:)=0.0_psn
			Jy(ind+1,:,:)=0.0_psn
			Jz(ind+1,:,:)=0.0_psn
		end if
	end subroutine CurrentBC_Outer
	
	subroutine CurrentBC_Outer1
		integer :: ind
		real(psn) :: f
		ind= BC_Rmax_Prtl-xborders(procxind(proc))+3
		!ind= BC_Rmax_Prtl-xborders(procxind(proc))+xmin
		if(ind.ge.2.and.ind.le.mx-1) then
			!To ensure that the current is deposited on right place for reflected particles
			!Jx(ind-1,:,:)=Jx(ind-1,:,:)-Jx(ind,:,:)
		 	Jy(ind-1,:,:)=Jy(ind-1,:,:)+Jy(ind+1,:,:)!*f
		 	Jz(ind-1,:,:)=Jz(ind-1,:,:)+Jz(ind+1,:,:)!*f 
			Jx(ind-1,:,:)=Jx(ind-1,:,:)-Jx(ind,:,:)
			Jx(ind,:,:)=0.0_psn
			Jy(ind+1,:,:)=0.0_psn
			Jz(ind+1,:,:)=0.0_psn
		end if
	end subroutine CurrentBC_Outer1
	
end module cyl_bc_cpu