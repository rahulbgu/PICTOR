module var_gpu
	use parameters
	use cudafor
	implicit none 
    
	integer, parameter :: NthreadsGPU=256
    real(psn),dimension(:,:,:), allocatable, device :: Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: Jx_gpu,Jy_gpu,Jz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: buffJx_gpu,buffJy_gpu,buffJz_gpu           
    real(psn),dimension(:,:,:), allocatable         :: buffJx_host,buffJy_host,buffJz_host
    real(psn), dimension(:,:,:), allocatable, device, target :: TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu
	
	integer, parameter :: Jwidth_gpu=16 
	real(psn),dimension(:,:,:,:), allocatable, device :: Jx_wide_gpu,Jy_wide_gpu,Jz_wide_gpu
	
	integer, dimension(:), allocatable, device :: cell_count_gpu
	integer :: Ncell_gpu
	
	!varaibles used to communicated Fld Data Between CPU and GPU
	integer, device :: mx_gpu,my_gpu,mz_gpu 
	type(dim3)         :: grid, tBlock
	type(dim3)         :: tBlock_gpu_YZedge,tGrid_gpu_YZedge
    type(dim3)         :: tBlock_gpu_ZXedge,tGrid_gpu_ZXedge
	type(dim3)         :: tBlock_gpu_XYedge,tGrid_gpu_XYedge
	type(dim3)         :: tBlock_gpu_YZedge1,tGrid_gpu_YZedge1	
	
	
	!-----------------------------------------------------------------------------------------------------
	! Particles
	!-----------------------------------------------------------------------------------------------------
	real(psn), dimension(:), allocatable, device :: flvrqm_gpu 
    !prtl arr size, host variables 
	integer, parameter :: Nchunk_prtl_gpu=64 ! total number of chunks
	integer, parameter :: chunk_size_prtl_gpu=gpu_prtl_arr_size/Nchunk_prtl_gpu
	integer :: np_gpu  !total number of particles on GPU
	integer, dimension(Nchunk_prtl_gpu) :: used_prtl_chunk
	integer            :: empty_prtl_chunk
	integer, parameter :: max_np_gpu=(Nchunk_prtl_gpu-1)*chunk_size_prtl_gpu  !one chunk must be empty for sorting 	
			 
	!parameters used in saving prtl using on CPU 
	integer, dimension(Nchunk_prtl_gpu) :: compact_chunk_offset_cpu		
	integer            :: buff_size_prtl_gpu=chunk_size_prtl_gpu ! size of the arrays that temporarily hold prtls coming/leaving the sub-domain
	
	integer, parameter :: pcount_gpu_size = 64
	integer, dimension(pcount_gpu_size), device :: pcount_gpu
	
    integer, dimension(:), allocatable, device:: flvp_gpu,tagp_gpu
    real(psn), dimension(:), allocatable, device:: qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu 
    integer, dimension(:), allocatable:: flvp_host,tagp_host
    real(psn), dimension(:), allocatable:: qp_host,xp_host,yp_host,zp_host,up_host,vp_host,wp_host,var1p_host
		
	integer, device :: np_send_gpu
	integer         :: np_recv_host	
    !-------------------------------------------------------------------
    !Curved BC related varaibles 
    !------------------------------------------------------------------- 
	real(psn), dimension(:,:,:), allocatable, device :: e_lx_gpu, e_ly_gpu , e_lz_gpu, b_arx_gpu, b_ary_gpu, b_arz_gpu ! used in conformal FDTD schemes 
	real(psn), dimension(:,:), allocatable, device :: usc_wtr_bx_gpu, usc_wtr_by_gpu, usc_wtr_bz_gpu ! wts used in the USC method, r (c)= row (column) of V matrix
	real(psn), dimension(:,:), allocatable, device :: usc_wtc_bx_gpu, usc_wtc_by_gpu, usc_wtc_bz_gpu
	integer, dimension(:,:,:), allocatable, device :: usc_fdtd_bx_gpu, usc_fdtd_by_gpu, usc_fdtd_bz_gpu ! 0: regular fdtd; positive integer, use the USC method, integer -> index in the wt. list 
	real(psn), dimension(:,:,:), allocatable, device :: usc_db1_gpu, usc_db2_gpu


end module var_gpu