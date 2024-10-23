module bc_attenuate
	use parameters
	use vars
#ifdef gpu
    use fields_gpu
#endif
    logical :: attn_em_vacumm = .false.
	real(dbpsn), dimension(6) :: attn_vac_pos	= (/ real(-huge(1),dbpsn), real(huge(1),dbpsn), real(-huge(1),dbpsn), real(huge(1),dbpsn), real(-huge(1),dbpsn), real(huge(1),dbpsn) /)
    integer, dimension(6) :: attn_vac_pos_ind = (/ -huge(1), huge(1), -huge(1), huge(1), -huge(1), huge(1) /)
contains 
	
	
	!------------------------------------------------------------------------
	! Attenuation of the EM field componnets in the presense of a flow field
	!------------------------------------------------------------------------
	subroutine Attenuate_EM_Flow(side)
		type(FlowFldType) :: flw
		integer  :: side
		real(psn):: b_x,b_y,b_z,vx,vy,vz, f0
		real(dbpsn) :: x0,y0,z0,xg,yg,zg
		real(dbpsn):: pos, dl, scale
	    integer :: i,j,k
		
		
		x0=xborders(procxind)-3
		y0=yborders(procyind)-3
		z0=zborders(proczind)-3
#ifdef twoD
        z0=0
#endif	
        flw = bc_face(side)%flw 
		pos = bc_face(side)%pos_fld
		dl = bc_face(side)%attn_thickness
		scale = bc_face(side)%attn_scale
		
		if(dl.eq.0) return
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
						  Ex(i,j,k)=(Ex(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,scale)+f0
 					  end if

 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0 = vx * b_z - vz * b_x
						  Ey(i,j,k)=(Ey(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,scale)+f0
 					  end if 

 					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then 
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
 						  call flw%Drift(xg,yg,zg,vx,vy,vz)
 						  f0= vy * b_x - vx * b_y
						  Ez(i,j,k)=(Ez(i,j,k)-f0)*AttenFactor(xg,yg,zg,pos,side,scale)+f0
 					  end if 
			  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bx(i,j,k)=(Bx(i,j,k)-b_x)*AttenFactor(xg,yg,zg,pos,side,scale)+b_x
 					  end if 

 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  By(i,j,k)=(By(i,j,k)-b_y)*AttenFactor(xg,yg,zg,pos,side,scale)+b_y
 					  end if 

 					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
 					  if(atnn_domain(xg,yg,zg,pos,side,dl)) then
 						  call flw%MagFld(xg,yg,zg,b_x,b_y,b_z)
						  Bz(i,j,k)=(Bz(i,j,k)-b_z)*AttenFactor(xg,yg,zg,pos,side,scale)+b_z
 					  end if 
		
	               end do
	          end do
	    end do
		
#ifdef gpu
	    call SendFullDomainEMFldtoGPU
#endif			
	
	end subroutine Attenuate_EM_Flow
	
	
	
	function AttenFactor(x,y,z,pos,side,dl)
		integer   :: side
		real(dbpsn) :: x,y,z,pos,dl
		real(psn) :: AttenFactor
		AttenFactor = 0.0_psn
		select case (side)
			case(1,2)
				AttenFactor = (abs(x-pos)/dl)**3
			case(3,4)
			    AttenFactor = (abs(y-pos)/dl)**3
			case(5,6)
				AttenFactor = (abs(z-pos)/dl)**3
		end select
		AttenFactor = min(1.0_psn,AttenFactor)
		AttenFactor = max(0.0_psn,AttenFactor)
		AttenFactor = 1.0_psn - AttenFactor
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
	! Attenuate EM field to zero beyond the boundaries 
	!------------------------------------------------------------------------
	
	subroutine Attenuate_EM_Vacuum
		real(dbpsn) :: x0,y0,z0,xg,yg,zg
		integer :: i,j,k , n
		logical :: inside
		
		
		!update position of the attenuating boundaries
		do n=1,6
			if(bc_face(n)%type_fld.eq.'attn') then 
				attn_vac_pos(n) = bc_face(n)%pos_fld
			    if(mod(n,2).eq.0) attn_vac_pos_ind(n) = ceiling(attn_vac_pos(n)) + 1
				if(mod(n,2).ne.0) attn_vac_pos_ind(n) = floor(attn_vac_pos(n)) - 1
			end if
		end do
		
		inside =  1+xborders(procxind)-3 .gt. attn_vac_pos_ind(1) .and. mx+xborders(procxind)-3 .lt. attn_vac_pos_ind(2) .and. & 
		          1+yborders(procyind)-3 .gt. attn_vac_pos_ind(3) .and. my+yborders(procyind)-3 .lt. attn_vac_pos_ind(4) .and. & 
				  1+zborders(proczind)-3 .gt. attn_vac_pos_ind(5) .and. mz+zborders(proczind)-3 .lt. attn_vac_pos_ind(6) 
		
		if(inside) return ! the entire eubdomain is in the interior and does not need to updated any fld
		
		
		x0=xborders(procxind)-3
		y0=yborders(procyind)-3
		z0=zborders(proczind)-3
#ifdef twoD
        z0=0
#endif

#ifdef gpu
		call RecvFullDomainEMFldFromGPU
#endif		
		
		do k=1,mz
	          do j=1,my
	               do i=1,mx
					  
			   		  inside =  i+xborders(procxind)-3 .gt. attn_vac_pos_ind(1) .and. i+xborders(procxind)-3 .lt. attn_vac_pos_ind(2) .and. & 
			   		            j+yborders(procyind)-3 .gt. attn_vac_pos_ind(3) .and. j+yborders(procyind)-3 .lt. attn_vac_pos_ind(4) .and. & 
			   				    k+zborders(proczind)-3 .gt. attn_vac_pos_ind(5) .and. k+zborders(proczind)-3 .lt. attn_vac_pos_ind(6)  
					  
					  if(inside) cycle			
					   
 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+z0;
					  Ex(i,j,k)=Ex(i,j,k)*AttenFactorVacuum(xg,yg,zg)
					  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
					  Ey(i,j,k)=Ey(i,j,k)*AttenFactorVacuum(xg,yg,zg)
 	 
 					  xg = i+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
					  Ez(i,j,k)=Ez(i,j,k)*AttenFactorVacuum(xg,yg,zg)
			  
 					  xg = i+x0; yg = j+0.5_dbpsn+y0; zg = k+0.5_dbpsn+z0;
                      Bx(i,j,k)=Bx(i,j,k)*AttenFactorVacuum(xg,yg,zg)

 					  xg = i+0.5_dbpsn+x0; yg = j+y0; zg = k+0.5_dbpsn+z0;
                      By(i,j,k)=By(i,j,k)*AttenFactorVacuum(xg,yg,zg)

 					  xg = i+0.5_dbpsn+x0; yg = j+0.5_dbpsn+y0; zg = k+z0;
                      Bz(i,j,k)=Bz(i,j,k)*AttenFactorVacuum(xg,yg,zg)
		
	               end do
	          end do
	    end do
		
#ifdef gpu
	    call SendFullDomainEMFldtoGPU
#endif	
		
	end subroutine Attenuate_EM_Vacuum
	
	function AttenFactorVacuum(x,y,z)
		real(dbpsn) :: x,y,z
		real(psn) :: AttenFactorVacuum
		
		AttenFactorVacuum = 0.0_psn
		
        if(x.lt.attn_vac_pos(1)) then                    
             if((y.le.(attn_vac_pos(4)+attn_vac_pos(1)-x)).and.(y.ge.(attn_vac_pos(3)-attn_vac_pos(1)+x))) then  
                 if(z.le.(attn_vac_pos(6)+attn_vac_pos(1)-x).and.(z.ge.(attn_vac_pos(5)-attn_vac_pos(1)+x))) then
					 
					 AttenFactorVacuum = (abs(x-attn_vac_pos(1))/bc_face(1)%attn_scale)**3
							 
                 end if						                              
            end if               
        else if(x.gt.attn_vac_pos(2)) then
            if((y.le.(attn_vac_pos(4)-attn_vac_pos(2)+x)).and.(y.ge.(attn_vac_pos(3)+attn_vac_pos(2)-x))) then						
                 if((z.le.(attn_vac_pos(6)-attn_vac_pos(2)+x)).and.(z.ge.(attn_vac_pos(5)+attn_vac_pos(2)-x))) then
					 
					 AttenFactorVacuum = (abs(x-attn_vac_pos(2))/bc_face(2)%attn_scale)**3
							    
                 end if 						                                                      
            end if     
        end if
        if(y.lt.attn_vac_pos(3)) then          
            if((y.lt.(attn_vac_pos(3)-attn_vac_pos(1)+x)).and.(y.lt.(attn_vac_pos(3)+attn_vac_pos(2)-x))) then						
                 if((z.ge.(attn_vac_pos(5)-attn_vac_pos(3)+y)).and.(z.le.(attn_vac_pos(6)+attn_vac_pos(3)-y))) then
					 
					  AttenFactorVacuum = (abs(y-attn_vac_pos(3))/bc_face(3)%attn_scale)**3
							 
                 end if						
            end if
        else if(y.gt.attn_vac_pos(4)) then
            if((y.gt.(attn_vac_pos(4)+attn_vac_pos(1)-x)).and.(y.gt.(attn_vac_pos(4)-attn_vac_pos(2)+x))) then						
                  if((z.ge.(attn_vac_pos(5)+attn_vac_pos(4)-y)).and.(z.le.(attn_vac_pos(6)-attn_vac_pos(4)+y))) then	
					 
					 AttenFactorVacuum = (abs(y-attn_vac_pos(4))/bc_face(4)%attn_scale)**3						 
							 
                 end if						
            end if
       end if
       
#ifndef twoD               
       if(z.lt.attn_vac_pos(5)) then          
            if((z.lt.(attn_vac_pos(5)-attn_vac_pos(1)+x)).and.(z.lt.(attn_vac_pos(5)+attn_vac_pos(2)-x))) then
                  if((z.lt.(attn_vac_pos(5)-attn_vac_pos(3)+y)).and.(z.lt.(attn_vac_pos(5)+attn_vac_pos(4)-y))) then					 
					
                	 AttenFactorVacuum = (abs(z-attn_vac_pos(5))/bc_face(5)%attn_scale)**3	
					 
                 end if
            end if
       else if(z.gt.attn_vac_pos(6)) then
            if(z.gt.(attn_vac_pos(6)+attn_vac_pos(1)-x).and.(z.gt.(attn_vac_pos(6)-attn_vac_pos(2)+x))) then
                 if((z.gt.(attn_vac_pos(6)+attn_vac_pos(3)-y)).and.(z.gt.(attn_vac_pos(6)-attn_vac_pos(4)+y))) then
					 
					 AttenFactorVacuum = (abs(z-attn_vac_pos(6))/bc_face(6)%attn_scale)**3

                 end if
            end if
       end if
#endif	

	   AttenFactorVacuum = min(1.0_psn,AttenFactorVacuum)
       AttenFactorVacuum = max(0.0_psn,AttenFactorVacuum)
	   
	   AttenFactorVacuum = 1.0_psn - AttenFactorVacuum
		
	end function AttenFactorVacuum
		
	
	
	

end module bc_attenuate