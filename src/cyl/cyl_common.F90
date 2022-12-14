module cyl_common
	use parameters
	use vars 
	implicit none
		 			 
#ifdef gpu
     type(dim3)         :: grid_ax, tBlock_ax
	 real(psn), dimension(:), allocatable, device :: Bz_ax_gpu
     real(psn), dimension(:,:), allocatable, device :: Ey_ax_gpu
	 real(psn), dimension(:,:), allocatable        :: Ey_ax_host	
#endif

contains 
!-----------------------------------------------------------------------------------------------
! The cylinderical version braches out of the main rectangular version by defining/redefining/constraining
! necessary variables. This module contains routines that are appropiately nested inside the routines in the rectangular version 
!-----------------------------------------------------------------------------------------------
	
	subroutine InitParam_cyl
		dtheta=2.0_psn*pi/ny
		ax_perm_area=4.0_psn/ny
		inv_dtheta=1.0_psn/dtheta
		rshift = -0.5_psn
		
		if(inc_axis) indepLBaxis=2 !indepLBaxis must be along the z-axis, if the axis in included
		
		if(inc_axis .and. nSubDomainsY.ne.1) STOP 'nSubDomainsY must be 1 if the cylinderical axis is included' 
		
	end subroutine InitParam_cyl
	
#ifdef gpu
    subroutine InitAll_cyl_gpu 
		allocate(Bz_ax_gpu(mz),Ey_ax_gpu(my,mz),Ey_ax_host(my,mz))
        tBlock_ax = dim3 (1 ,64 ,4)
        grid_ax = dim3(1, ceiling(real(my)/tBlock_ax%y), ceiling(real(mz)/tBlock_ax%z)) 
	end subroutine InitAll_cyl_gpu 
#endif 

    !number of electrons in the subdomain Nelc, assuming uniform distribution, is redefined	
	subroutine Nelc_cyl
		 real(psn) :: rmin,rmax
		 
		 rmin=xborders(procxind) + rshift
		 rmax=xborders(procxind+1) + rshift
		 if(procxind.eq.0) rmin=0
		 
		 Nelc=0.5_psn*epc*dtheta*(rmax**2-rmin**2)*(ymax-ymin)*(zmax-zmin)

	end subroutine Nelc_cyl
	
	!if weights are used, weight of each particles is proportional to r\dtheta if num. of prtl. ib cells is independent of r  
	subroutine azimuthal_weight_prtl(wt,x)
		real(psn) :: wt
		real(dbpsn) :: x, r
		r=max(x+rshift,0.0_dbpsn)
		wt = wt * r * dtheta
	end subroutine azimuthal_weight_prtl

   	
	

end module cyl_common