module bc_surf_cylinder
	use parameters
	use vars
	use mem_prtl
	implicit none
	
	integer :: bc_cyl_axis ! axis of the cylinder
	real(dbpsn) :: bc_cyl_x0 , bc_cyl_y0, bc_cyl_z0 
	real(dbpsn) :: bc_cyl_R_fld ! radius of the cylinder for EM field attn./refl.
	real(dbpsn) :: bc_cyl_R_prtl ! radius of the cylinder for particle boundary 
	
	character (len=4) :: bc_cyl_type_fld, bc_cyl_type_prtl ! type of boundary condition
	real(psn) :: attn_cyl_scale, attn_cyl_thickness 
	procedure(scalar_global_dp), pointer :: distanceFromAxis
contains 
	
	subroutine SetBC_Cylinder(axis,x0,y0,z0,R,fld_type,prtl_type,attn_scale,attn_thickness)
		integer :: axis
		real(dbpsn) :: x0,y0,z0,R
		character (len=4), optional :: fld_type, prtl_type
		
		if(present(fld_type)) bc_cyl_type_fld = fld_type
		if(present(prtl_type))  bc_cyl_type_prtl = prtl_type
		if(present(attn_scale)) attn_cyl_scale = attn_scale
		if(present(attn_thickness)) attn_cyl_thickness = attn_thickness
		
		bc_cyl_x0 = x0
	    bc_cyl_y0 = y0
	    bc_cyl_z0 = z0
	    bc_cyl_R_fld = R
		
		bc_cyl_axis = axis
		
		if(axis.eq.0) distanceFromAxis => distanceFromAxis_x
		if(axis.eq.1) distanceFromAxis => distanceFromAxis_y
		if(axis.eq.2) distanceFromAxis => distanceFromAxis_z
		
		
	end subroutine SetBC_Cylinder
	
	subroutine EnforceBC_Fld_Cylinder
		if(bc_cyl_type_fld.eq.'attn') call BC_Attn_Fld_Cylinder
	end subroutine EnforceBC_Fld_Cylinder
	
	subroutine EnforceBC_Prtl_Cylinder
		if(bc_cyl_type_prtl.eq.'remv') call BC_RemovePrtl_Cylinder
	end subroutine EnforceBC_Prtl_Cylinder
	
	!----------------------------------------------------------------------------------------------------------
	! Attenuate electric and magentic field outside the cylidercal boundary; include planar faces as well if appropiate
	!----------------------------------------------------------------------------------------------------------
	subroutine BC_Attn_Fld_Cylinder
		real(dbpsn) :: x0,y0,z0,xg,yg,zg
		integer :: i,j,k
		
		x0=xborders(procxind)-3
		y0=yborders(procyind)-3
		z0=zborders(proczind)-3
#ifdef twoD
        z0=0
#endif
		
		do k=1,mz
	          do j=1,my
	               do i=1,mx
					  
! 			   		  inside =  i+xborders(procxind)-3 .gt. attn_vac_pos_ind(1) .and. i+xborders(procxind)-3 .lt. attn_vac_pos_ind(2) .and. &
! 			   		            j+yborders(procyind)-3 .gt. attn_vac_pos_ind(3) .and. j+yborders(procyind)-3 .lt. attn_vac_pos_ind(4) .and. &
! 			   				    k+zborders(proczind)-3 .gt. attn_vac_pos_ind(5) .and. k+zborders(proczind)-3 .lt. attn_vac_pos_ind(6)
!
! 					  if(inside) cycle		
					   
 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+z0;
					  Ex(i,j,k)=Ex(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)
					  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  Ey(i,j,k)=Ey(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)
 	 
 					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  Ez(i,j,k)=Ez(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)
			  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
                      Bx(i,j,k)=Bx(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)

 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
                      By(i,j,k)=By(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)

 					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
                      Bz(i,j,k)=Bz(i,j,k)*AttenFactorVacuum_Cyl(xg,yg,zg)
		
	               end do
	          end do
	    end do
	
	end subroutine BC_Attn_Fld_Cylinder
	
	function AttenFactorVacuum_Cylinder(x,y,z)
		real(dbpsn) :: x,y,z
		real(dbpsn) :: d,d1,d2
		real(psn) :: AttenFactorVacuum_Cylinder
		
		AttenFactorVacuum_Cylinder = 0.0_psn
		d = distanceFromAxis(x,y,z) - bc_cyl_R_fld
		
		!the other two boundaries orthogonal to the axis
		d1 = distanceFromTop(x,y,z,bc_cyl_axis)
		d2 = distanceFromBottom(x,y,z,bc_cyl_axis)
		d = max(d,d1)
		d = max(d,d2)
		
		if( d .gt. 0) AttenFactorVacuum_Cylinder = (d/attn_cyl_scale)**3
		  
	    AttenFactorVacuum_Cylinder = min(1.0_psn,AttenFactorVacuum_Cylinder)
        AttenFactorVacuum_Cylinder = max(0.0_psn,AttenFactorVacuum_Cylinder)
	   
	    AttenFactorVacuum_Cylinder = 1.0_psn - AttenFactorVacuum_Cylinder
		
	end function AttenFactorVacuum_Cylinder
	
	function distanceFromAxis_x(x,y,z)
		real(dbpsn) :: x,y,z
		real(dbpsn) :: distanceFromAxis_x
		distanceFromAxis_x = sqrt( (y-bc_cyl_y0)**2 + (z-bc_cyl_z0)**2)
	end function distanceFromAxis_x
	
	function distanceFromAxis_y(x,y,z)
		real(dbpsn) :: x,y,z
		real(dbpsn) :: distanceFromAxis_y
		distanceFromAxis_y = sqrt( (x-bc_cyl_x0)**2 + (z-bc_cyl_z0)**2) 
	end function distanceFromAxis_y
	
	function distanceFromAxis_z(x,y,z)
		real(dbpsn) :: x,y,z
		real(dbpsn) :: distanceFromAxis_z
		distanceFromAxis_z = sqrt( (x-bc_cyl_x0)**2 + (y-bc_cyl_y0)**2) 
	end function distanceFromAxis_z
	
	function distanceFromTop(x,y,z,axis)
		integer :: axis
		real(dbpsn) :: x,y,z
		
		distanceFromTop = 0 
		select case (axis)
	    	case(0)
				distanceFromTop = x-bc_face(2*axis+2)%pos_fld
			case(1)
			    distanceFromTop = y-bc_face(2*axis+2)%pos_fld
			case(2)
			    distanceFromTop = z-bc_face(2*axis+2)%pos_fld
		end select
	end function distanceFromTop
	
	function distanceFromBottom(x,y,z,axis)
		integer :: axis
		real(dbpsn) :: x,y,z
		
		distanceFromBottom = 0 
		select case (axis)
	    	case(0)
				distanceFromBottom = bc_face(2*axis+1)%pos_fld -x
			case(1)
			    distanceFromBottom = bc_face(2*axis+1)%pos_fld -y
			case(2)
			    distanceFromBottom = bc_face(2*axis+1)%pos_fld -z
		end select
	end function distanceFromBottom
	
	!----------------------------------------------------------------------------------------------------------
	! Remove particles outside the cylinder
	!----------------------------------------------------------------------------------------------------------
	subroutine BC_RemovePrtl_Cylinder
		real(dbpsn) :: xg,yg,zg , ox,oy,oz, d 
		integer :: n, count
		
		ox=xborders(procxind)-3
		oy=yborders(procyind)-3
		oz=zborders(proczind)-3
		
		count = 0
		do n=1,used_prtl_arr_size
			if(flv(n).eq.0) cycle
			xg = xp(n) + ox 
			yg = yp(n) + oy
			zg = zp(n) + oz
			d = distanceFromAxis(xg,yg,zg) -bc_cyl_R_prtl 
			if(d.gt.0) then 
				call DeletePrtl(n)
				count = count +1 
			end if
			
		end do 
		np=np -count
	end subroutine BC_RemovePrtl_Cylinder
		

end module bc_surf_cylinder