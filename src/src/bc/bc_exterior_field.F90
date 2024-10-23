module bc_exterior_field
	use parameters
	use vars
contains 
	
	!----------------------------------------------------------------------------------------------
	!  Set elecric and magentic fields corresponding to a flow profile 
	!----------------------------------------------------------------------------------------------
	subroutine Set_BC_EM_FlowField(side)
		integer :: side
		call SetMagFld(bc_face(side)%flw%MagFld, bc_face(side)%pos_fld, side)
		call SetFlowElcFld(bc_face(side)%flw%MagFld, bc_face(side)%flw%Drift, bc_face(side)%pos_fld, side)
	end subroutine Set_BC_EM_FlowField
	
	!----------------------------------------------------------------------------------------------
	!  Set magentic field corresponding to a flow profile 
	!----------------------------------------------------------------------------------------------
	subroutine SetMagFld(MagFld,pos,side)
		procedure(vector_global) :: MagFld
		real(dbpsn):: pos
		integer  :: side
		real(psn):: b_x,b_y,b_z
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
					  
					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=b_x
					  end if 

					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=b_y
					  end if 

					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=b_z
					  end if 
  
				  end do
			  end do  
	     end do
		
	end subroutine SetMagFld
	
	!----------------------------------------------------------------------------------------------
	!  Set electric field corresponding to a flow profile 
	!----------------------------------------------------------------------------------------------
	subroutine SetFlowElcFld(MagFld,Drift,pos,side)
		procedure(vector_global) :: MagFld, Drift
		real(dbpsn):: pos
		integer  :: side
		real(psn)::  b_x,b_y,b_z,vx,vy,vz
		real(dbpsn) :: xg,yg,zg, x0,y0,z0
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
		  
					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call Drift(xg,yg,zg,vx,vy,vz)
						  Ex(i,j,k)= vz * b_y - vy * b_z 
					  end if

					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  if(exterior(xg,yg,zg,pos,side)) then
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call Drift(xg,yg,zg,vx,vy,vz)
						  Ey(i,j,k)= vx * b_z - vz * b_x
					  end if 

					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  if(exterior(xg,yg,zg,pos,side)) then 
						  call MagFld(xg,yg,zg,b_x,b_y,b_z)
						  call Drift(xg,yg,zg,vx,vy,vz)
						  Ez(i,j,k)= vy * b_x - vx * b_y
					  end if 
  
				  end do
			  end do  
	     end do
		
	end subroutine SetFlowElcFld
	
	logical function exterior(xg,yg,zg,pos,side)
		real(dbpsn) :: xg,yg,zg,pos
		integer :: side
		exterior = .false.
		if(side.eq.1 .and. xg.le. pos) exterior = .true. 
        if(side.eq.2 .and. xg.ge. pos) exterior = .true. 
		if(side.eq.3 .and. yg.le. pos) exterior = .true. 
	    if(side.eq.4 .and. yg.ge. pos) exterior = .true. 
		if(side.eq.5 .and. zg.le. pos) exterior = .true. 
	    if(side.eq.6 .and. zg.ge. pos) exterior = .true. 
	end function exterior
	

end module bc_exterior_field