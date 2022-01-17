module cyl_vars
     use parameters
	 use vars
#ifdef gpu
     use cudafor
#endif	 
	 implicit none 
	 
	 real(psn) :: dtheta,inv_dtheta,ax_perm_area
	 real(psn), dimension(:), allocatable :: Bz_ax
	 
	 real(psn), dimension (0:nSubDomainsX) :: rborders 
	 integer :: Nelc_uniform_cyl
	 
	 integer, dimension(:), allocatable :: rinp_count_axis,rintp_count_axis,rpcross_axis,rtpcross_axis
	 type(particle), dimension(:,:), allocatable :: rinp_axis,routp_axis
		 
	 real(psn) :: BC_Rmin_Prtl, BC_Rmax_Prtl
	 real(psn)   :: BC_Rmin_Fld, BC_Rmax_Fld
	 character (len=4) :: BC_Rmin_Fld_Type, BC_Rmax_Fld_Type
	 character (len=4) :: BC_Rmin_Prtl_Type, BC_Rmax_Prtl_Type

#ifdef gpu
     type(dim3)         :: grid_ax, tBlock_ax
	 real(psn), dimension(:), allocatable, device :: Bz_ax_gpu
     real(psn), dimension(:,:), allocatable, device :: Ey_ax_gpu
	 real(psn), dimension(:,:), allocatable        :: Ey_ax_host	
#endif	 
	 
end module cyl_vars