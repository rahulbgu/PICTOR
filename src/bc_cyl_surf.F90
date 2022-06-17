module bc_cyl_surf
	use parameters
	use vars
	real(psn) :: x0,y0,r,r2,r_prtl,r2_prtl
	real(psn), dimension(:,:), allocatable :: lx,ly,arz
	integer, dimension(:,:), allocatable :: bcell
	integer :: ndiv=100
	procedure(vector_global),pointer   :: MagFld
contains 
	
	!------------------------------------------------------------------------------
	! muti-cell conducting material
	! E Fld witing the material at r>r_c is gradually brought down to zero
	!------------------------------------------------------------------------------
	subroutine ThickCondCyl_OuterBC(x0,y0,Rc,delr)
		real(psn) :: x0,y0,Rc,delr
		real(psn) :: ox,oy,oz
		real(psn) :: r
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		oz=zborders(proczind(proc))-3
#ifdef twoD
        oz=0
#endif	
		
		do k=1,mz			
			do j=1,my
				do i=1,mx
					r=xyzToR(real(i+0.5_psn,psn)+ox,real(j,psn)+oy,real(k,psn)+oz,x0,y0)
					if(r.gt.Rc) call attenuateEfld(Ex(i,j,k),r,Rc,delr)
				    
					r=xyzToR(real(i,psn)+ox,real(j+0.5_psn,psn)+oy,real(k,psn)+oz,x0,y0)
					if(r.gt.Rc) call attenuateEfld(Ey(i,j,k),r,Rc,delr)
					
				    r=xyzToR(real(i,psn)+ox,real(j,psn)+oy,real(k+0.5_psn,psn)+oz,x0,y0)
					if(r.gt.Rc) call attenuateEfld(Ez(i,j,k),r,Rc,delr)
				
				end do 
			end do 
		end do
		
	end subroutine ThickCondCyl_OuterBC
	
	subroutine attenuateEfld(fld,r,Rc,delr)
		real(psn) :: fld, r, Rc, delr
		real(psn) :: f 
		f=1.0_psn - (r-Rc)/delr
		f= max(0.0_psn,f)
		f=0.0_psn
		fld = fld*f  
	end subroutine attenuateEFld
	
	function xyzToR(x,y,z,x0,y0)
		real(psn) :: x,y,z,x0,y0
		real(psn) :: xyzToR
		
		xyzToR=sqrt((x-x0)**2+(y-y0)**2)
	end function xyzToR
	
	!------------------------------------------------------------------------------
	! Conformal FTDT method
	!------------------------------------------------------------------------------
	subroutine EnforceBC_cyl_surf
		call SetFldBC_cyl_surf
		call SetPrtlBC_cyl_surf
	end subroutine EnforceBC_cyl_surf
	
	subroutine SetBC_cyl_surf(v1,v2,v3,v4,mag_setup)
		real(psn) :: v1,v2,v3,v4
		interface 
			subroutine mag_setup(x,y,z,fx,fy,fz)
				import :: psn
				real(psn) :: x,y,z,fx,fy,fz
			end subroutine
		end interface 
		
		x0=v1
		y0=v2
		r=v3
		r2=v3**2
		r_prtl=v4
		r2_prtl=v4**2
		MagFld=>mag_setup
	    call InitBC_LenArea
	end subroutine SetBC_cyl_surf
	
	subroutine InitBC_LenArea
	    integer :: i,j,k
		real(psn) :: ox,oy, xc,yc
		real(psn), dimension(mx,my) :: dist
		
		if(allocated(lx)) then
			deallocate(lx,ly,arz, bcell) 
		end if 
			allocate(lx(mx,my),ly(mx,my),arz(mx,my),bcell(mx,my))
    	lx=1.0_psn
		ly=1.0_psn
		arz=0.0_psn
		bcell=0
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		
	    !distance^2 of all cell corners from the center 
		dist=0
		do j=1,my
			do i=1,mx
				dist(i,j) = dist2_cyl( real(i,psn)+ox, real(j,psn)+oy )
			end do 
		end do 
	
		do j=1,my-1
			do i=1,mx-1
				!if( dist(i,j).gt.0 .and. dist(i,j+1).gt.0 .and. dist(i+1,j).gt.0 .and. dist(i+1,j+1).gt.0 ) continue
				!if( dist(i,j).lt.0 .and. dist(i,j+1).lt.0 .and. dist(i+1,j).lt.0 .and. dist(i+1,j+1).lt.0 ) continue 
			    
				bcell(i,j)=1
				xc = real(i+0.5_psn) +ox -x0
				yc = real(j+0.5_psn) +oy -y0
				
				if(yc.gt.0 .and. abs(xc).lt.abs(yc)) call area_intx(i,j,ox,oy,1)
				if(yc.lt.0 .and. abs(xc).lt.abs(yc)) call area_intx(i,j,ox,oy,-1)
				if(xc.gt.0 .and. abs(xc).ge.abs(yc)) call area_inty(i,j,ox,oy,1)
				if(xc.lt.0 .and. abs(xc).ge.abs(yc)) call area_inty(i,j,ox,oy,-1)
				
				call edge_lenx(i,j,real(i,psn)+ox,real(i+1,psn)+ox, real(j,psn)+oy)
				call edge_leny(i,j,real(j,psn)+oy,real(j+1,psn)+oy, real(i,psn)+ox)
			end do 
		end do 
		
	end subroutine InitBC_LenArea
	
	subroutine edge_lenx(i,j,x1,x2,y)
		integer :: i,j
		real(psn) :: x1,x2,y,p1,p2
		if(r2 - (y-y0)**2 .lt.0) then 
			lx(i,j)=0
			return
		end if 
		p1= max(x1, x0 -sqrt(r2 - (y-y0)**2))
		p2= min(x2, x0 +sqrt(r2 - (y-y0)**2)) 
		
		lx(i,j)=p2-p1
		lx(i,j)=max(0.0,lx(i,j))
		lx(i,j)=min(1.0,lx(i,j))		
	end subroutine edge_lenx
	
	subroutine edge_leny(i,j,y1,y2,x)
		integer :: i,j
		real(psn) :: y1,y2, x,p1,p2
		if(r2 - (x-x0)**2 .lt.0) then 
		    ly(i,j)=0
			return
		end if
		p1= max(y1, y0 -sqrt(r2 - (x-x0)**2))
		p2= min(y2, y0 +sqrt(r2 - (x-x0)**2))  
		
		ly(i,j)=p2-p1
		ly(i,j)=max(0.0,ly(i,j))
		ly(i,j)=min(1.0,ly(i,j))	
	end subroutine edge_leny
	
	subroutine area_intx(i,j,ox,oy,ysign)
		integer :: i,j,n, ysign
		real(psn) :: ox,oy,x1,x2,y1,xi,yi,dx,dy
		
		x1=real(i,psn)+ox
		x2=real(i+1,psn)+ox
		y1=real(j,psn)+oy
		dx= 1.0_psn/ndiv
		
		if(ysign.lt.0) y1=y1+1
		
		do n=1,ndiv
			xi= x1 + (n-1)*dx + dx/2.0_psn
			
			if(r2 - (xi-x0)**2 .lt.0) cycle 
			yi= y0 + ysign*sqrt(r2 - (xi-x0)**2)
			dy= (yi - y1)*ysign
		
			if(dy.gt.1.0_psn) dy=1.0_psn
			if(dy.lt.0.0_psn) dy=0.0_psn
			arz(i,j)= arz(i,j) + dx*dy
		end do
	end subroutine area_intx
	
	subroutine area_inty(i,j,ox,oy,xsign)
		integer :: i,j,n, xsign
		real(psn) :: ox,oy,y1,y2,x1,xi,yi,dx,dy
		
		y1=real(j,psn)+oy
		y2=real(j+1,psn)+oy
		x1=real(i,psn)+ox
		dy= 1.0_psn/ndiv
		
		if(xsign.lt.0) x1=x1+1
		
		do n=1,ndiv
			yi= y1 + (n-1)*dy + dy/2.0_psn
			
			if(r2 - (yi-y0)**2 .lt.0) cycle		
			xi= x0 + xsign*sqrt(r2 - (yi-y0)**2)
			dx= (xi - x1)*xsign
	
			if(dx.gt.1.0_psn) dx=1.0_psn
			if(dx.lt.0.0_psn) dx=0.0_psn
			arz(i,j)= arz(i,j) + dx*dy
		end do
	end subroutine area_inty
	
	subroutine SetFldBC_cyl_surf
		integer :: i,j,k
		real(psn) :: ox,oy,oz
		real(psn) :: b_x,b_y,b_z
		real(psn), dimension(mx,my,mz) :: Ex_cond, Ey_cond 
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		oz=zborders(proczind(proc))-3
#ifdef twoD
        oz=0
#endif		
		
		do k=1,mz			
			do j=2,my-1
				do i=2,mx-1
					if(dist2_cyl(real(i+0.5_psn,psn)+ox,real(j,psn)+oy).ge.0) Ex_cond(i,j,k)=radialEx(i,j,k,ox,oy) 
					if(dist2_cyl(real(i,psn)+ox,real(j+0.5_psn,psn)+oy).ge.0) Ey_cond(i,j,k)=radialEy(i,j,k,ox,oy) 
					
					if(dist2_cyl(real(i+0.5_psn,psn)+ox,real(j,psn)+oy).ge.0) Ex(i,j,k)=Ex_cond(i,j,k)
					if(dist2_cyl(real(i,psn)+ox,real(j+0.5_psn,psn)+oy).ge.0) Ey(i,j,k)=Ey_cond(i,j,k)
					
					if(dist2_cyl(real(i,psn)+ox,real(j,psn)+oy).ge.0) Ez(i,j,k)=0

! 					if(dist2_cyl(real(i,psn)+ox,real(j+0.5_psn,psn)+oy).ge.0) then
! 						call MagFld(real(i,psn)+ox,real(j+0.5_psn,psn)+oy,real(k+0.5_psn)+oz,b_x,b_y,b_z)
! 						Bx(i,j,k)=b_x
! 					end if
!
! 					if(dist2_cyl(real(i+0.5_psn,psn)+ox,real(j,psn)+oy).ge.0) then
! 						call MagFld(real(i+0.5_psn,psn)+ox,real(j,psn)+oy,real(k+0.5_psn)+oz,b_x,b_y,b_z)
! 						By(i,j,k)=b_y
! 					end if
!
! 					if(dist2_cyl(real(i+0.5_psn,psn)+ox,real(j+0.5_psn,psn)+oy).ge.0) then
! 						call MagFld(real(i+0.5_psn,psn)+ox,real(j+0.5_psn,psn)+oy,real(k)+oz,b_x,b_y,b_z)
! 						Bz(i,j,k)=b_z
! 					end if

! 					if(bcell(i,j).eq.1) then
! 						if(dist2_cyl(real(i+0.5_psn,psn)+ox,real(j+0.5_psn,psn)+oy).lt.0) then
! 							Bz(i,j,k)= Bz(i,j,k) -fldc*( (ly(i+1,j)-1.0_psn)*Ey(i+1,j,k) - (ly(i,j)-1.0_psn)*Ey(i,j,k)) + fldc*( (lx(i,j+1)-1.0_psn)*Ex(i,j+1,k) - (lx(i,j)-1.0_psn)*Ex(i,j,k))
! 						    !Bz(i,j,k)= Bz(i,j,k)/arz(i,j)
! 						end if
! 					end if
				end do  
			end do 
		end do 
		
		Jx(:,:,1)=lx(:,:)
		Jy(:,:,1)=ly(:,:)
		Jz(:,:,1)=arz(:,:)
	end subroutine SetFldBC_cyl_surf
	
	function radialEx(i,j,k,ox,oy)
		integer :: i,j,k
		real(psn) :: ox,oy,ct,st,e_x,e_y,e_r,xg,yg,rg
		real(psn) :: radialEx
		
		xg=i+0.5_psn+ox -x0
		yg=j+oy -y0
		rg=xg**2 + yg**2
		ct=xg/rg
		st=yg/rg
		e_y=0.25_psn*(Ey(i,j,k)+Ey(i,j-1,k)+Ey(i+1,j,k)+Ey(i+1,j-1,k))
		e_x=Ex(i,j,k)
		e_r= e_x*ct + e_y*st
		radialEx=e_r*ct
	end function radialEx
	
	function radialEy(i,j,k,ox,oy)
		integer :: i,j,k
		real(psn) :: ox,oy,ct,st,e_x,e_y,e_r,xg,yg,rg
		real(psn) :: radialEy
		
		xg=i+ox -x0
		yg=j+0.5_psn+oy -y0
		rg=xg**2 + yg**2
		ct=xg/rg
		st=yg/rg
		e_x=0.25_psn*(Ex(i-1,j,k)+Ex(i,j,k)+Ex(i-1,j+1,k)+Ex(i,j+1,k))
		e_y=Ex(i,j,k)
		e_r= e_x*ct + e_y*st
		radialEy=e_r*st
	end function radialEy
	
	
	function dist2_cyl(x,y)
		real(psn) :: x,y
		real(psn) :: dist2_cyl
		dist2_cyl= ((x-x0)**2+(y-y0)**2) -r2
	end function dist2_cyl
	
	
	subroutine SetPrtlBC_cyl_surf
		integer :: n 
		real(psn) :: ox,oy
		real(psn) :: vr,vt,ct,st,rp 
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			rp=(xp(n)+ox -x0)**2 + (yp(n)+oy -y0)**2
			if( rp .gt. r2_prtl) then
				ct=xp(n)/rp
				st=yp(n)/rp
				
				vr= up(n)*ct + vp(n)*st
				vt=-up(n)*st + vp(n)*ct
				
				vr=-vr
				up(n)= vr*ct - vt*st
				vp(n)= vr*st + vt*ct 
			end if
		end do 
	end subroutine SetPrtlBC_cyl_surf
	
	
		
end module bc_cyl_surf