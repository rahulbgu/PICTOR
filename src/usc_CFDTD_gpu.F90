module usc_CFDTD_gpu
	use parameters
	use vars 
	use var_gpu
	use deposit_gpu
contains
	 
	 subroutine Init_USC_CFDTD_gpu
 		allocate(e_lx_gpu(mx,my,mz),e_ly_gpu(mx,my,mz),e_lz_gpu(mx,my,mz),b_arx_gpu(mx,my,mz),b_ary_gpu(mx,my,mz),b_arz_gpu(mx,my,mz))
		e_lx_gpu = e_lx; e_ly_gpu = e_ly; e_lz_gpu = e_lz; b_arx_gpu = b_arx; b_ary_gpu = b_ary; b_arz_gpu = b_arz;  
		
 		allocate(usc_db1_gpu(mx,my,mz),usc_db2_gpu(mx,my,mz))
 		usc_db1_gpu=usc_db1; usc_db2_gpu=usc_db2; 
 		
		allocate(usc_fdtd_bx_gpu(mx,my,mz),usc_fdtd_by_gpu(mx,my,mz),usc_fdtd_bz_gpu(mx,my,mz))
 		usc_fdtd_bx_gpu=usc_fdtd_bx; usc_fdtd_by_gpu=usc_fdtd_by; usc_fdtd_bz_gpu=usc_fdtd_bz;
		
        allocate(usc_wtr_bx_gpu(size(usc_wtr_bx,1),9), usc_wtr_by_gpu(size(usc_wtr_by,1),9), usc_wtr_bz_gpu(size(usc_wtr_bz,1),9)) ! r= row 
        allocate(usc_wtc_bx_gpu(size(usc_wtc_bx,1),9), usc_wtc_by_gpu(size(usc_wtc_by,1),9), usc_wtc_bz_gpu(size(usc_wtc_bz,1),9)) ! c= column
		
		usc_wtr_bx_gpu=usc_wtr_bx; usc_wtr_by_gpu=usc_wtr_by; usc_wtr_bz_gpu=usc_wtr_bz;  
		usc_wtc_bx_gpu=usc_wtc_bx; usc_wtc_by_gpu=usc_wtc_by; usc_wtc_bz_gpu=usc_wtc_bz;  
	 end subroutine Init_USC_CFDTD_gpu
	 
	 subroutine DeallocateCurvedBCVars_gpu
         deallocate(e_lx_gpu,e_ly_gpu,e_lz_gpu, b_arx_gpu, b_ary_gpu, b_arz_gpu) 
         deallocate(usc_db1_gpu,usc_db2_gpu) 
		 deallocate(usc_fdtd_bx_gpu,usc_fdtd_by_gpu,usc_fdtd_bz_gpu)
		 deallocate(usc_wtr_bx_gpu,usc_wtr_by_gpu,usc_wtr_bz_gpu)
		 deallocate(usc_wtc_bx_gpu,usc_wtc_by_gpu,usc_wtc_bz_gpu)
	 end subroutine DeallocateCurvedBCVars_gpu
 	!------------------------------------------------------------------------------
 	! Modified FDTD for the USC method: GPU version
 	!------------------------------------------------------------------------------
	 subroutine UpdateBfieldHalf_USC_gpu
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db1_gpu,mx,my,mz)
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db2_gpu,mx,my,mz)
		 
		 call Calc_dbx_gpu<<<grid,tBlock>>>(usc_db1_gpu,Ey_gpu,Ez_gpu,e_ly_gpu,e_lz_gpu,mx,my,mz,fld_halfc)
		 call gather_dbx_bc_cells_gpu<<<grid,tBlock>>>(usc_db1_gpu,usc_db2_gpu,usc_wtr_bx_gpu,usc_fdtd_bx_gpu,mx,my,mz)
		 call gather_dbx_bc_cells_gpu<<<grid,tBlock>>>(usc_db2_gpu,usc_db1_gpu,usc_wtc_bx_gpu,usc_fdtd_bx_gpu,mx,my,mz)
		 call add_usc_db1_gpu<<<grid,tBlock>>>(Bx_gpu,usc_db1_gpu,mx,my,mz)
		 
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db1_gpu,mx,my,mz)
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db2_gpu,mx,my,mz)
		 
		 call Calc_dby_gpu<<<grid,tBlock>>>(usc_db1_gpu,Ex_gpu,Ez_gpu,e_lx_gpu,e_lz_gpu,mx,my,mz,fld_halfc)
		 call gather_dby_bc_cells_gpu<<<grid,tBlock>>>(usc_db1_gpu,usc_db2_gpu,usc_wtr_by_gpu,usc_fdtd_by_gpu,mx,my,mz)
		 call gather_dby_bc_cells_gpu<<<grid,tBlock>>>(usc_db2_gpu,usc_db1_gpu,usc_wtc_by_gpu,usc_fdtd_by_gpu,mx,my,mz)
		 call add_usc_db1_gpu<<<grid,tBlock>>>(By_gpu,usc_db1_gpu,mx,my,mz)
		 
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db1_gpu,mx,my,mz)
		 call reset_db_gpu<<<grid,tBlock>>>(usc_db2_gpu,mx,my,mz)
		 
		 call Calc_dbz_gpu<<<grid,tBlock>>>(usc_db1_gpu,Ex_gpu,Ey_gpu,e_lx_gpu,e_ly_gpu,mx,my,mz,fld_halfc)
		 call gather_dbz_bc_cells_gpu<<<grid,tBlock>>>(usc_db1_gpu,usc_db2_gpu,usc_wtr_bz_gpu,usc_fdtd_bz_gpu,mx,my,mz)
		 call gather_dbz_bc_cells_gpu<<<grid,tBlock>>>(usc_db2_gpu,usc_db1_gpu,usc_wtc_bz_gpu,usc_fdtd_bz_gpu,mx,my,mz)
		 call add_usc_db1_gpu<<<grid,tBlock>>>(Bz_gpu,usc_db1_gpu,mx,my,mz)		 
	 end subroutine UpdateBfieldHalf_USC_gpu
	 
	 subroutine SetExteriorEfld_USC_gpu
		 call SetExteriorEfld_USC_kernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,e_lx_gpu,e_ly_gpu,e_lz_gpu,mx,my,mz)
	 end subroutine SetExteriorEfld_USC_gpu
	 
	 attributes(global) subroutine SetExteriorEfld_USC_kernel(Ex,Ey,Ez,e_lx,e_ly,e_lz,mx,my,mz)
	 	 integer, value :: mx,my,mz
		 real, dimension(mx,my,mz) :: Ex,Ey,Ez,e_lx,e_ly,e_lz 
		 integer :: i,j,k
	
		 i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		 j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		 k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
         k=1
#endif
         if(i.le.mx .and. j.le.my .and. k.le.mz) then 
			  if(e_lx(i,j,k).eq.0) Ex(i,j,k)=0.0
			  if(e_ly(i,j,k).eq.0) Ey(i,j,k)=0.0
			  if(e_lz(i,j,k).eq.0) Ez(i,j,k)=0.0
		 end if 
	 
	 end subroutine SetExteriorEfld_USC_kernel 
	 
	 attributes(global) subroutine add_usc_db1_gpu(arr,db1,mx,my,mz)
		 real, dimension(mx,my,mz) :: arr, db1
		 integer, value :: mx,my,mz 
		 integer :: i,j,k
	
		 i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		 j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		 k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
         k=1
#endif	
         if(i.le.mx .and. j.le.my .and. k.le.mz) arr(i,j,k)=arr(i,j,k)+db1(i,j,k)
		 
	 end subroutine add_usc_db1_gpu 
	 attributes(global) subroutine reset_db_gpu(arr,mx,my,mz)
		 real, dimension(mx,my,mz) :: arr
		 integer, value :: mx,my,mz 
		 integer :: i,j,k
	
		 i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		 j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		 k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
         k=1          
#endif		 
         if(i.le.mx .and. j.le.my .and. k.le.mz) arr(i,j,k) = 0.0	
	 end subroutine reset_db_gpu 
	 
	 attributes(global) subroutine Calc_dbx_gpu(arr,Ey,Ez,e_ly,e_lz,mx,my,mz,fld_halfc)
	 	real, dimension(mx,my,mz) :: arr, Ey,Ez,e_ly,e_lz
		real, value :: fld_halfc
		integer, value :: mx,my,mz 
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		
#ifdef twoD
        k=1
		if(i.gt.mx .or. j.gt.my-1) return
        arr(i,j,k)=-fld_halfc*(e_lz(i,j+1,k)*Ez(i,j+1,k)-e_lz(i,j,k)*Ez(i,j,k))
#else
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
        if(i.gt.mx .or. j.gt.my-1 .or. k.gt.mz-1) return
		arr(i,j,k)=-fld_halfc*(e_lz(i,j+1,k)*Ez(i,j+1,k)-e_lz(i,j,k)*Ez(i,j,k))+&
					fld_halfc*(e_ly(i,j,k+1)*Ey(i,j,k+1)-e_ly(i,j,k)*Ey(i,j,k))
#endif		
     end subroutine Calc_dbx_gpu
	 
  	attributes(global) subroutine gather_dbx_bc_cells_gpu(arr_in,arr_out,wt,usc_fdtd_bx,mx,my,mz)
  		real, dimension(mx,my,mz) :: arr_in, arr_out
		integer, dimension(mx,my,mz) :: usc_fdtd_bx
  		real, dimension(:,:)   :: wt
 		integer, value :: mx,my,mz
  		integer :: i,j,k,n
		
  		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
  		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
  		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
        k=1 
 		if(i.gt.mx .or. j.lt.2 .or. j.gt.my-1) return	
#else 
 	    if(i.gt.mx .or. j.lt.2 .or. j.gt.my-1 .or. k.lt.2 .or. k.gt.mz-1) return	
#endif			
 		if(usc_fdtd_bx(i,j,k).gt.0) then 
 			n = usc_fdtd_bx(i,j,k)
#ifdef twoD		
	 		arr_out(i,j,k) = wt(n,1)*arr_in(i,j,k) + wt(n,2)*arr_in(i,j-1,k) + wt(n,3)*arr_in(i,j+1,k)
#else		
	 		arr_out(i,j,k) = wt(n,1)*arr_in(i,j,k) + wt(n,2)*arr_in(i,j,k-1) + wt(n,3)*arr_in(i,j,k+1) + wt(n,4)*arr_in(i,j-1,k) + wt(n,5)*arr_in(i,j+1,k)+&
	 		                 wt(n,6)*arr_in(i,j-1,k-1) + wt(n,7)*arr_in(i,j-1,k+1) + wt(n,8)*arr_in(i,j+1,k-1) + wt(n,9)*arr_in(i,j+1,k+1)
#endif
 		end if 
		
 	end subroutine gather_dbx_bc_cells_gpu
	 
	 attributes(global) subroutine Calc_dby_gpu(arr,Ex,Ez,e_lx,e_lz,mx,my,mz,fld_halfc)
	 	real, dimension(mx,my,mz) :: arr, Ex,Ez,e_lx,e_lz
		real, value :: fld_halfc
		integer, value :: mx,my,mz 
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		
#ifdef twoD
        k=1
        if(i.gt.mx-1 .or. j.gt.my) return
		arr(i,j,k)=fld_halfc*(e_lz(i+1,j,k)*Ez(i+1,j,k)-e_lz(i,j,k)*Ez(i,j,k))
#else
        k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
		if(i.gt.mx-1 .or. j.gt.my .or. k.gt.mz-1) return
		arr(i,j,k)=-fld_halfc*(e_lx(i,j,k+1)*Ex(i,j,k+1)-e_lx(i,j,k)*Ex(i,j,k))+&
				    fld_halfc*(e_lz(i+1,j,k)*Ez(i+1,j,k)-e_lz(i,j,k)*Ez(i,j,k))
#endif	
	 end subroutine Calc_dby_gpu 
	 
 	attributes(global) subroutine gather_dby_bc_cells_gpu(arr_in,arr_out,wt,usc_fdtd_by,mx,my,mz)
 		real, dimension(mx,my,mz) :: arr_in, arr_out
		integer, dimension(mx,my,mz) :: usc_fdtd_by
 		real, dimension(:,:)   :: wt
		integer, value :: mx,my,mz
 		integer :: i,j,k,n
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
        k=1 
		if(i.lt.2 .or. i.gt.mx-1 .or. j.gt.my) return	
#else 
	    if(i.lt.2 .or. i.gt.mx-1 .or. j.gt.my .or. k.lt.2 .or. k.gt.mz-1) return	
#endif			
		if(usc_fdtd_by(i,j,k).gt.0) then 
			n = usc_fdtd_by(i,j,k)			
#ifdef twoD		
	 		arr_out(i,j,k) = wt(n,1)*arr_in(i,j,k) + wt(n,2)*arr_in(i-1,j,k) + wt(n,3)*arr_in(i+1,j,k)
#else		
	 		arr_out(i,j,k) = wt(n,1)*arr_in(i,j,k) + wt(n,2)*arr_in(i-1,j,k) + wt(n,3)*arr_in(i+1,j,k) + wt(n,4)*arr_in(i,j,k-1) + wt(n,5)*arr_in(i,j,k+1)+&
	 		                 wt(n,6)*arr_in(i-1,j,k-1) + wt(n,7)*arr_in(i+1,j,k-1) + wt(n,8)*arr_in(i-1,j,k+1) + wt(n,9)*arr_in(i+1,j,k+1)
#endif
		end if 
		
	end subroutine gather_dby_bc_cells_gpu
	 
	 attributes(global) subroutine Calc_dbz_gpu(arr,Ex,Ey,e_lx,e_ly,mx,my,mz,fld_halfc)
	 	real, dimension(mx,my,mz) :: arr, Ex,Ey,e_lx,e_ly
		real, value :: fld_halfc
		integer, value :: mx,my,mz 
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
        k=1 
#endif		
        if(i.gt.mx-1 .or. j.gt.my-1 .or. k.gt.mz) return
        arr(i,j,k)=-fld_halfc*(e_ly(i+1,j,k)*Ey(i+1,j,k)-e_ly(i,j,k)*Ey(i,j,k))+&
	                fld_halfc*(e_lx(i,j+1,k)*Ex(i,j+1,k)-e_lx(i,j,k)*Ex(i,j,k))
	end subroutine Calc_dbz_gpu
	
	attributes(global) subroutine gather_dbz_bc_cells_gpu(arr_in,arr_out,wt,usc_fdtd_bz,mx,my,mz)
		real, dimension(mx,my,mz) :: arr_in, arr_out
		integer, dimension(mx,my,mz) :: usc_fdtd_bz
		real, dimension(:,:)   :: wt
		integer, value :: mx,my,mz
		integer :: i,j,k,n
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#ifdef twoD
        k=1 
#endif			
		if(i.lt.2 .or. i.gt.mx-1 .or. j.lt.2 .or. j.gt.my-1 .or. k.gt.mz) return
		if(usc_fdtd_bz(i,j,k).gt.0) then 
			n=usc_fdtd_bz(i,j,k)
			arr_out(i,j,k)=wt(n,1)*arr_in(i,j,k) + wt(n,2)*arr_in(i-1,j,k) + wt(n,3)*arr_in(i+1,j,k) + wt(n,4)*arr_in(i,j-1,k) + wt(n,5)*arr_in(i,j+1,k)+&
			               wt(n,6)*arr_in(i-1,j-1,k) + wt(n,7)*arr_in(i+1,j-1,k) + wt(n,8)*arr_in(i-1,j+1,k) + wt(n,9)*arr_in(i+1,j+1,k) 
		end if 
    end subroutine gather_dbz_bc_cells_gpu
	 
	!-------------------------------------------------------------
	! Boundary condition for the particles
	!------------------------------------------------------------- 

	subroutine PrtlBC_cyl_gpu(x0,y0,Rc)
		real(psn) :: x0,y0,Rc
		integer   :: kc, indi,indf
		real(psn) :: eps, ox,oy
	    
		eps = epsilon(x0) 
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			call PrtlBC_cyl_gpu_kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,qp_gpu, indi,indf,x0,y0,Rc,ox,oy,eps, Jx_wide_gpu,Jy_wide_gpu,Jz_wide_gpu,mx,my,mz,Jwidth_gpu, qi,c)
		end do 
		
		call ReduceJwide
		
	end subroutine PrtlBC_cyl_gpu

	attributes(global) subroutine PrtlBC_cyl_gpu_kernel(xp,yp,zp,up,vp,wp,qp, indi,indf,x0,y0,Rc,ox,oy,eps ,Jx,Jy,Jz,mx,my,mz,jwidth, qi,c)
		real, dimension(:) :: xp,yp,zp,up,vp,wp,qp
		integer, value :: indi,indf
		real, value :: qi,c,x0,y0,Rc,ox,oy, eps
		integer, value :: mx,my,mz,jwidth
		real, dimension(mx,my,mz,jwidth) :: Jx,Jy,Jz
		real :: im_xf,im_yf,im_zf
		real :: xpi,ypi,zpi,g,f
		real :: vr,vt,ct,st,rp,rb
		real :: xb, yb, zb
		integer        :: n
		
        n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
        if(n.gt.indf) return
	
		if(qp(n).eq.0) return
		rp=sqrt((xp(n)+ox -x0)**2 + (yp(n)+oy -y0)**2)
		if(rp.lt.Rc) return
	
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
	

		call circle_line_intersection_gpu(Rc,xpi+ox-x0,ypi+oy-y0,xp(n)+ox-x0,yp(n)+oy-y0,xb,yb,eps)
		rb=sqrt( xb**2 + yb**2)
		ct=xb/rb
		st=yb/rb
		
		f = sqrt((xp(n)-xpi)**2 + (yp(n)-ypi)**2) 
		g = sqrt((xb+x0-ox -xpi)**2 + (yb+y0-oy -ypi)**2) 
		if(f.ge.eps) then 
			g = min(g/f,1.0)
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
		
		
		
		!reverse of the old trajectory
        call DepositCurrentPIC_GPU(xp(n),yp(n),zp(n),xpi,ypi,zpi,qp(n),n,qi, Jx,Jy,Jz,mx,my,mz,jwidth) 
		
	    !new trajetory
		call DepositCurrentPIC_GPU(xpi,ypi,zpi,xb,yb,zb,qp(n),n,qi, Jx,Jy,Jz,mx,my,mz,jwidth) 
		call DepositCurrentPIC_GPU(xb,yb,zb,im_xf,im_yf,im_zf,qp(n),n,qi, Jx,Jy,Jz,mx,my,mz,jwidth) 
		
		xp(n)=im_xf
		yp(n)=im_yf
	
	
    end subroutine PrtlBC_cyl_gpu_kernel
	
	
	attributes(device) subroutine circle_line_intersection_gpu(R,xi,yi,xf,yf,xc,yc,eps)
		real :: R,xi,yi,xf,yf,xc,yc
		real :: dsq, A , B , C , AB2, xm, ym , m
		real :: eps
		
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
		if( ((xi-xc)**2+(yi-yc)**2).gt.0.25 .or. ((xf-xc)**2+(yf-yc)**2).gt.0.25 ) then 
			xc = xi
			yc = yi
		end if 
				
	end subroutine circle_line_intersection_gpu

end module usc_CFDTD_gpu