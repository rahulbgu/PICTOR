module bc_curved_surf
	! the current implementation includes only a conducting cylinderical surface as an outer boundary
	use parameters
	use vars
	use deposit
	use usc_CFDTD
#ifdef gpu
    use var_gpu
    use usc_CFDTD_gpu
#endif	
	implicit none 
	character (len=3) :: surf_type
	real(psn) :: x0, y0, Rc
	
contains 

		
	!------------------------------------------------------------------------------
	! Cylinder
	!------------------------------------------------------------------------------
	subroutine Set_CondBC_Cylinder(cx,cy,r)
		real(psn) :: cx,cy,r
		x0=cx
		y0=cy
		Rc=r
		curved_bc=.true.
		surf_type='cyl'
		call InitCurvedBC
	end subroutine Set_CondBC_Cylinder

	 
	 subroutine InitCurvedBC
		
		if(.not.curved_bc) return
		
        call Init_USC_CFDTD
		
		if(surf_type.eq.'cyl') then !only the cylinder is implemented
		     call Init_LenArea_cyl 
		end if
		
		call InitWeight_USC
		
#ifdef gpu
        call Init_USC_CFDTD_gpu
#endif
	 end subroutine InitCurvedBC
	 	 
 	 subroutine PrtlBC_curved_surf
        if(surf_type.eq.'cyl') then
#ifdef gpu
        	call PrtlBC_cyl_gpu(x0,y0,Rc+1)
			return
#endif		
 			!call PrtlBC_cyl(x0,y0,Rc+1)	
			call PrtlBC_cyl(x0,y0,Rc+1.5)			
		end if

 	 end subroutine PrtlBC_curved_surf
	 	 

 	!---------------------------------------------------------------
	!
	!              Cylinder 
	!
	!---------------------------------------------------------------
	
	!---------------------------------------------------------------
 	! Calculate the edge length and area of the facets that lie within the domain
 	!---------------------------------------------------------------
	subroutine Init_LenArea_cyl
		integer :: i,j,k
		real(psn) :: ox,oy,xc,yc
			
			
		ox=xborders(procxind)-3
		oy=yborders(procyind)-3
			
		do j=1,my
			do i=1,mx
				
				!main grid where e is at the center of grid edges
				xc = real(i) +ox -x0
				yc = real(j) +oy -y0
				
				if(yc.ge.0 .and. abs(xc).lt.abs(yc)) call area_intx_cyl(b_arz(i,j,1),real(i,psn)+ox,real(i+1,psn)+ox,real(j,psn)+oy,1)
				if(yc.le.0 .and. abs(xc).lt.abs(yc)) call area_intx_cyl(b_arz(i,j,1),real(i,psn)+ox,real(i+1,psn)+ox,real(j,psn)+oy,-1)
				if(xc.ge.0 .and. abs(xc).ge.abs(yc)) call area_inty_cyl(b_arz(i,j,1),real(j,psn)+oy,real(j+1,psn)+oy,real(i,psn)+ox,1)
				if(xc.le.0 .and. abs(xc).ge.abs(yc)) call area_inty_cyl(b_arz(i,j,1),real(j,psn)+oy,real(j+1,psn)+oy,real(i,psn)+ox,-1)
				
				!if(b_arz(i,j,1).lt.0.01) b_arz(i,j,1)=0.0_psn
				
				call calc_lx_cyl(e_lx(i,j,1),real(i,psn)+ox,real(i+1,psn)+ox, real(j,psn)+oy)
				call calc_ly_cyl(e_ly(i,j,1),real(j,psn)+oy,real(j+1,psn)+oy, real(i,psn)+ox)
			    call calc_lz_cyl(e_lz(i,j,1),real(i,psn)+ox,real(j,psn)+oy)
								
			end do 
		end do
		
#ifndef twoD 
!copy the calculated area and legnth to other k=const planes 
		do k=2,mz 
			do j=1,my
				do i=1,mx
					b_arz(i,j,k) = b_arz(i,j,1)
					e_lx(i,j,k) = e_lx(i,j,1)
					e_ly(i,j,k) = e_ly(i,j,1)
					e_lz(i,j,k) = e_lz(i,j,1)
				end do 
			end do 
		end do
#endif 		
			
		b_arx=e_ly
		b_ary=e_lx
		
		
	end subroutine Init_LenArea_cyl

	subroutine calc_lx_cyl(lx,x1,x2,y)
		real(psn) :: lx, x1,x2,y,p1,p2,r2
		
		r2=Rc**2
		if(r2 - (y-y0)**2 .lt.0) then 
			lx=0.0_psn
			return
		end if 
		
		p1= max(x1, x0 -sqrt(r2 - (y-y0)**2))
		p2= min(x2, x0 +sqrt(r2 - (y-y0)**2)) 
		
		lx=p2-p1
		lx=max(0.0_psn,lx)
		lx=min(1.0_psn,lx)		
	end subroutine calc_lx_cyl
	
	subroutine calc_ly_cyl(ly,y1,y2,x)
		real(psn) :: ly, y1,y2, x,p1,p2,r2
		
		r2=Rc**2
		if(r2 - (x-x0)**2 .lt.0) then 
		    ly=0.0_psn
			return
		end if
		
		p1= max(y1, y0 -sqrt(r2 - (x-x0)**2))
		p2= min(y2, y0 +sqrt(r2 - (x-x0)**2))  
		
		ly=p2-p1
		ly=max(0.0_psn,ly)
		ly=min(1.0_psn,ly)	
	end subroutine calc_ly_cyl
	
	subroutine calc_lz_cyl(lz,x,y)
		real(psn) :: lz, x,y
		if((x-x0)**2+(y-y0)**2>Rc**2) then 
			lz=0.0_psn
		else 
			lz=1.0_psn
		end if
	end subroutine calc_lz_cyl
	
	subroutine area_intx_cyl(ar,x1,x2,y1,ysign)
		integer :: n, ysign
		integer :: ndiv=1000
		real(psn) :: ar,x1,x2,y1,xi,yi,dx,dy,r2
		
		dx= 1.0_psn/ndiv
		r2=Rc**2
		
		if(ysign.lt.0) y1=y1+1
		
		ar=0.0_psn
		do n=1,ndiv
			xi= x1 + (n-1)*dx + dx/2.0_psn
			
			if(r2 - (xi-x0)**2 .lt.0) cycle 
			yi= y0 + ysign*sqrt(r2 - (xi-x0)**2)
			dy= (yi - y1)*ysign
		
			if(dy.gt.1.0_psn) dy=1.0_psn
			if(dy.lt.0.0_psn) dy=0.0_psn
			ar= ar + dx*dy
		end do
	end subroutine area_intx_cyl
	
	subroutine area_inty_cyl(ar,y1,y2,x1,xsign)
		integer :: n, xsign
		integer :: ndiv=1000
		real(psn) :: ar,y1,y2,x1,xi,yi,dx,dy,r2
		
		dy= 1.0_psn/ndiv
		r2=Rc**2
		
		if(xsign.lt.0) x1=x1+1
		
		ar=0.0_psn
		do n=1,ndiv
			yi= y1 + (n-1)*dy + dy/2.0_psn
			
			if(r2 - (yi-y0)**2 .lt.0) cycle		
			xi= x0 + xsign*sqrt(r2 - (yi-y0)**2)
			dx= (xi - x1)*xsign
	
			if(dx.gt.1.0_psn) dx=1.0_psn
			if(dx.lt.0.0_psn) dx=0.0_psn
			ar= ar + dx*dy
		end do
	end subroutine area_inty_cyl
	
	!--------------------------------------------------------------
	! BC for particles; deposit current due the mirror image and swap if the original particle is outside the domain
	!--------------------------------------------------------------
	subroutine PrtlBC_cyl(x0,y0,Rc)
		real(psn) :: x0,y0,Rc
		real(psn) :: ox,oy
		real(psn) :: im_xf,im_yf,im_zf
		real(psn) :: xpi,ypi,zpi,g,f
		real(psn) :: vr,vt,ct,st,rp,rb
		real(psn) :: xb, yb, zb, eps 
		integer :: n
		
		ox=xborders(procxind)-3
		oy=yborders(procyind)-3
		eps = epsilon(x0) 
		
		do n=1,used_prtl_arr_size
			if(flvp(n).eq.0) cycle
			rp=sqrt((xp(n)+ox -x0)**2 + (yp(n)+oy -y0)**2)
			if(rp.lt.Rc) cycle
			
			!find initial positon by moving particle back in time
			g=c/sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
			xpi=xp(n)-up(n)*g
			ypi=yp(n)-vp(n)*g
			zpi=zp(n)-wp(n)*g
			
			!image of the final position
			ct=(xp(n)+ox -x0)/rp
			st=(yp(n)+oy -y0)/rp
			im_xf=(Rc+(Rc-rp))*ct -ox+x0
			im_yf=(Rc+(Rc-rp))*st -oy+y0
			im_zf=zp(n)
			
			!if(t.gt.16) print*,'before reflection',im_xi,im_xf,im_yi,im_yf,'xp yp',xp(n),yp(n),'up vp',up(n),vp(n),'xpi ypi',xpi,ypi

			
			!if(rp.lt.Rc)  call DepositCurrentPIC(im_xi,im_yi,im_zi,im_xf,im_yf,im_zf,qp(n))  
	
			
			
			! replace the particle with it's image if the original particle is outside the domain
			
				

			call circle_line_intersection(Rc,xpi+ox-x0,ypi+oy-y0,xp(n)+ox-x0,yp(n)+oy-y0,xb,yb,eps)
			rb=sqrt( xb**2 + yb**2)
			ct=xb/rb
			st=yb/rb
			
			f = sqrt((xp(n)-xpi)**2 + (yp(n)-ypi)**2) 
			g = sqrt((xb+x0-ox -xpi)**2 + (yb+y0-oy -ypi)**2) 
			if(f.ge.eps) then 
				g = min(g/f,1.0_psn)
			else
				g=0.0
			endif
			zb = zpi + g*(zp(n)-zpi) 
			
			
			vr= up(n)*ct + vp(n)*st
			vt=-up(n)*st + vp(n)*ct

			vr=-vr
			up(n)= vr*ct - vt*st
			vp(n)= vr*st + vt*ct

			xb = xb - ox + x0
			yb = yb - oy + y0  

			
			!!!Warning :: 3D need to take care of z as wll 

			!if(t.gt.33) print*,'xb,yb'xb,yb
			!reverse of the old trajectory
            call DepositCurrentPIC(xp(n),yp(n),zp(n),xpi,ypi,zpi,qp(n)) 
			
		    !new trajetory
			call DepositCurrentPIC(xpi,ypi,zpi,xb,yb,zb,qp(n)) 
			call DepositCurrentPIC(xb,yb,zb,im_xf,im_yf,im_zf,qp(n)) 
			
			xp(n)=im_xf
			yp(n)=im_yf
				
				
			
			

			
		end do
	end subroutine PrtlBC_cyl
	
	subroutine circle_line_intersection(R,xi,yi,xf,yf,xc,yc,eps)
		real(psn) :: R,xi,yi,xf,yf,xc,yc
		real(psn) :: dsq, A , B , C , AB2, xm, ym , m
		real(psn) :: eps
		
		!default result 
		xc = xi
		yc = yi
		


		A = yi - yf
		B = xf - xi
		C  = xi*yf - xf*yi 
		AB2 = A*A + B*B
		
		if(A.le.eps .and. B.le.eps) return
		if(AB2.le.eps) return
		
		xm = - A*C/AB2
		ym = - B*C/AB2

		dsq = R*R - (C*C/AB2)
		
		if(dsq.le.eps) return
		
		m = sqrt(dsq/AB2)
		
		if(-B*(xf-xm)+A*(yf-ym).ge.0) then
			xc = xm - B*m
			yc = ym + A*m
		else
			xc = xm + B*m
			yc = ym - A*m
		end if 
		
		!Failsafe: dispacement/step does not exceed 0.5 in the explicit scheme and (xc,yc) is in between i and f
		if( ((xi-xc)**2+(yi-yc)**2).gt.0.25_psn .or. ((xf-xc)**2+(yf-yc)**2).gt.0.25_psn ) then 
			xc = xi
			yc = yi
		end if 
		
		!if(t.gt.32) print*,xi,xf,yi,yf,'final' ,sqrt(xi*xi+yi*yi),sqrt(xf*xf+yf*yf),sqrt(xm*xm+ym*ym),R, m, dsq
				
	end subroutine circle_line_intersection

	
end module bc_curved_surf