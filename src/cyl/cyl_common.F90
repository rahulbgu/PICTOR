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
		
#ifndef twoD		
		if(inc_axis) indepLBaxis=2 !indepLBaxis must be along the z-axis, if the axis in included
#endif		
		
		if(inc_axis .and. nSubDomainsY.ne.1) STOP 'nSubDomainsY must be 1 if the cylinderical axis is included' 
		if(grid_dy.ne.1.0_psn)  STOP ' grid_dy must be set to 1.0_psn if the cylinderical grid is used' 
		
	end subroutine InitParam_cyl
	
#ifdef gpu
    subroutine InitAll_cyl_gpu 
		allocate(Bz_ax_gpu(mz),Ey_ax_gpu(my,mz),Ey_ax_host(my,mz))
        tBlock_ax = dim3 (1 ,64 ,4)
        grid_ax = dim3(1, ceiling(real(my)/tBlock_ax%y), ceiling(real(mz)/tBlock_ax%z)) 
	end subroutine InitAll_cyl_gpu 
#endif 

    !number of electrons in the subdomain Nelc, assuming uniform distribution
	integer function Nelc_uniform_cyl(x1,x2,y1,y2,z1,z2)
		 real(dbpsn) :: x1,x2,y1,y2,z1,z2 !range of domain; global cordinate
		 real(dbpsn) :: rmin,rmax
		 
		 rmin=x1 + rshift
		 rmax=x2 + rshift
		 rmin=max(rmin,0.0_dbpsn)
		 
		 Nelc_uniform_cyl=0.5_psn*epc*dtheta*(rmax**2-rmin**2)*(y2-y1)*(z2-z1)

	end function Nelc_uniform_cyl
	
	
	!if weights are used, weight of each particles is proportional to r\dtheta if num. of prtl. ib cells is independent of r  
	subroutine azimuthal_weight_prtl(wt,x)
		real(psn) :: wt
		real(dbpsn) :: x, r
		r=max(x+rshift,0.0_dbpsn)
		wt = wt * r * dtheta
	end subroutine azimuthal_weight_prtl

	!a random xglobal between x1 and x2; r1 is a random number
	subroutine random_xglobal_cyl(r1,xglobal,x1,x2)
		real(dbpsn) :: xglobal, x1,x2, r1
		real(dbpsn) :: rmin,rmax
		
	    rmin=x1 + rshift
	    rmax=x2 + rshift
	    rmin=max(rmin,0.0_dbpsn)
		
		xglobal =  sqrt(rmin**2 + r1*(rmax**2 - rmin**2)) -rshift
		
	end subroutine random_xglobal_cyl


   	
	

end module cyl_common