module var_gpu
	use parameters
	use cudafor
	implicit none 
    
	integer, parameter :: gpu_test_prtl_arr_size=256*8 ! max. no. of test particles on GPU: Test Particles will be merged with prtls
	integer            :: buff_size_test_prtl_gpu=256 ! same as above but for the test particles :: To be removed 
	
	integer, parameter :: NthreadsGPU=256
    real(psn),dimension(:,:,:), allocatable, device :: Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: Jx_gpu,Jy_gpu,Jz_gpu
    real(psn),dimension(:,:,:), allocatable, device :: buffJx_gpu,buffJy_gpu,buffJz_gpu           
    real(psn),dimension(:,:,:), allocatable         :: buffJx_host,buffJy_host,buffJz_host
    real(psn), dimension(:,:,:), allocatable, device, target :: TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu
	integer :: AvailGPUSlots, AvailGPUSlotsTestPrtl
	
	integer, parameter :: Jwidth_gpu=16 
	real(psn),dimension(:,:,:,:), allocatable, device :: Jx_wide_gpu,Jy_wide_gpu,Jz_wide_gpu
	
	integer, dimension(:), allocatable, device :: cell_count_gpu
	integer :: Ncell_gpu
	
	!varaibles used to communicated Fld Data Between CPU and GPU
	integer :: mx_gpu_host,my_gpu_host,mz_gpu_host
	integer, device :: mx_gpu,my_gpu,mz_gpu 
	type(dim3)         :: grid, tBlock
	type(dim3)         :: tBlock_gpu_YZedge,tGrid_gpu_YZedge
    type(dim3)         :: tBlock_gpu_ZXedge,tGrid_gpu_ZXedge
	type(dim3)         :: tBlock_gpu_XYedge,tGrid_gpu_XYedge
	type(dim3)         :: tBlock_gpu_YZedge1,tGrid_gpu_YZedge1	
	
	
	real(psn),dimension(:,:,:), allocatable, device :: lEx_send_gpu,rEx_send_gpu,tEx_send_gpu,bEx_send_gpu,uEx_send_gpu,dEx_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lEx_recv_gpu,rEx_recv_gpu,tEx_recv_gpu,bEx_recv_gpu,uEx_recv_gpu,dEx_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lEx_send_host,rEx_send_host,tEx_send_host,bEx_send_host,uEx_send_host,dEx_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lEx_recv_host,rEx_recv_host,tEx_recv_host,bEx_recv_host,uEx_recv_host,dEx_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lEy_send_gpu,rEy_send_gpu,tEy_send_gpu,bEy_send_gpu,uEy_send_gpu,dEy_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lEy_recv_gpu,rEy_recv_gpu,tEy_recv_gpu,bEy_recv_gpu,uEy_recv_gpu,dEy_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lEy_send_host,rEy_send_host,tEy_send_host,bEy_send_host,uEy_send_host,dEy_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lEy_recv_host,rEy_recv_host,tEy_recv_host,bEy_recv_host,uEy_recv_host,dEy_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lEz_send_gpu,rEz_send_gpu,tEz_send_gpu,bEz_send_gpu,uEz_send_gpu,dEz_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lEz_recv_gpu,rEz_recv_gpu,tEz_recv_gpu,bEz_recv_gpu,uEz_recv_gpu,dEz_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lEz_send_host,rEz_send_host,tEz_send_host,bEz_send_host,uEz_send_host,dEz_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lEz_recv_host,rEz_recv_host,tEz_recv_host,bEz_recv_host,uEz_recv_host,dEz_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lBx_send_gpu,rBx_send_gpu,tBx_send_gpu,bBx_send_gpu,uBx_send_gpu,dBx_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lBx_recv_gpu,rBx_recv_gpu,tBx_recv_gpu,bBx_recv_gpu,uBx_recv_gpu,dBx_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lBx_send_host,rBx_send_host,tBx_send_host,bBx_send_host,uBx_send_host,dBx_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lBx_recv_host,rBx_recv_host,tBx_recv_host,bBx_recv_host,uBx_recv_host,dBx_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lBy_send_gpu,rBy_send_gpu,tBy_send_gpu,bBy_send_gpu,uBy_send_gpu,dBy_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lBy_recv_gpu,rBy_recv_gpu,tBy_recv_gpu,bBy_recv_gpu,uBy_recv_gpu,dBy_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lBy_send_host,rBy_send_host,tBy_send_host,bBy_send_host,uBy_send_host,dBy_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lBy_recv_host,rBy_recv_host,tBy_recv_host,bBy_recv_host,uBy_recv_host,dBy_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lBz_send_gpu,rBz_send_gpu,tBz_send_gpu,bBz_send_gpu,uBz_send_gpu,dBz_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lBz_recv_gpu,rBz_recv_gpu,tBz_recv_gpu,bBz_recv_gpu,uBz_recv_gpu,dBz_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lBz_send_host,rBz_send_host,tBz_send_host,bBz_send_host,uBz_send_host,dBz_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lBz_recv_host,rBz_recv_host,tBz_recv_host,bBz_recv_host,uBz_recv_host,dBz_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lJx_send_gpu,rJx_send_gpu,tJx_send_gpu,bJx_send_gpu,uJx_send_gpu,dJx_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lJx_recv_gpu,rJx_recv_gpu,tJx_recv_gpu,bJx_recv_gpu,uJx_recv_gpu,dJx_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lJx_send_host,rJx_send_host,tJx_send_host,bJx_send_host,uJx_send_host,dJx_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lJx_recv_host,rJx_recv_host,tJx_recv_host,bJx_recv_host,uJx_recv_host,dJx_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lJy_send_gpu,rJy_send_gpu,tJy_send_gpu,bJy_send_gpu,uJy_send_gpu,dJy_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lJy_recv_gpu,rJy_recv_gpu,tJy_recv_gpu,bJy_recv_gpu,uJy_recv_gpu,dJy_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lJy_send_host,rJy_send_host,tJy_send_host,bJy_send_host,uJy_send_host,dJy_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lJy_recv_host,rJy_recv_host,tJy_recv_host,bJy_recv_host,uJy_recv_host,dJy_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lJz_send_gpu,rJz_send_gpu,tJz_send_gpu,bJz_send_gpu,uJz_send_gpu,dJz_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lJz_recv_gpu,rJz_recv_gpu,tJz_recv_gpu,bJz_recv_gpu,uJz_recv_gpu,dJz_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lJz_send_host,rJz_send_host,tJz_send_host,bJz_send_host,uJz_send_host,dJz_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lJz_recv_host,rJz_recv_host,tJz_recv_host,bJz_recv_host,uJz_recv_host,dJz_recv_host
	
	real(psn),dimension(:,:,:), allocatable, device :: lCurr1_send_gpu,rCurr1_send_gpu,tCurr1_send_gpu,bCurr1_send_gpu,uCurr1_send_gpu,dCurr1_send_gpu
	real(psn),dimension(:,:,:), allocatable, device :: lCurr1_recv_gpu,rCurr1_recv_gpu,tCurr1_recv_gpu,bCurr1_recv_gpu,uCurr1_recv_gpu,dCurr1_recv_gpu
	real(psn),dimension(:,:,:), pinned, allocatable :: lCurr1_send_host,rCurr1_send_host,tCurr1_send_host,bCurr1_send_host,uCurr1_send_host,dCurr1_send_host
	real(psn),dimension(:,:,:), pinned, allocatable :: lCurr1_recv_host,rCurr1_recv_host,tCurr1_recv_host,bCurr1_recv_host,uCurr1_recv_host,dCurr1_recv_host
	
	
	
    real(psn):: xmin_host,xmax_host,ymax_host,ymin_host,zmax_host,zmin_host,xlen_host,ylen_host,zlen_host ! physical boundaries of particles at the gpu (local cordinate)
    integer  :: xmin1_host,xmax1_host,ymin1_host,ymax1_host, zmin1_host, zmax1_host
	
    real(psn), device :: xmin_gpu,xmax_gpu,ymax_gpu,ymin_gpu,zmax_gpu,zmin_gpu ! physical boundaries of particles at the gpu (local cordinate)
    integer, device  :: xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu, zmin1_gpu, zmax1_gpu
	real(psn) :: xmin0_host 
	
	!constants 
	real(psn), device :: Bx_ext0_gpu,By_ext0_gpu,Bz_ext0_gpu ! constant external magnetic field 
    real(psn), device :: fldc_gpu
    real(psn), device :: fld_halfc_gpu
    real(psn), device :: cinv_gpu
    real(psn), device :: sqc_gpu
    real(psn), device :: c_gpu
	real(psn), device :: wtm1_gpu,wt0_gpu,wtp1_gpu
	real(psn), dimension(:), allocatable, device :: flvrqm_gpu 
	
	
	!-----------------------------------------------------------------------------------------------------
	! Particles
	!-----------------------------------------------------------------------------------------------------
	!particle array 
    integer, dimension(:), allocatable, device:: flvp_gpu,tagp_gpu
    real(psn), dimension(:), allocatable, device:: qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu
    real(psn), dimension(:), allocatable, device:: qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu
	integer, dimension(:), allocatable, device :: flvp_send_gpu,tagp_send_gpu
    real(psn), dimension(:), allocatable, device:: qp_recv_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,var1p_recv_gpu
	integer, dimension(:), allocatable, device :: flvp_recv_gpu,tagp_recv_gpu
	integer, device :: np_send_gpu
	
	!prtl arr on host used for communication
    real(psn), dimension(:), pinned, allocatable:: qp_host_send,xp_host_send,yp_host_send,zp_host_send,up_host_send,vp_host_send,wp_host_send,var1p_host_send
	integer, dimension(:), pinned, allocatable  :: flvp_host_send,tagp_host_send
    real(psn), dimension(:), pinned, allocatable:: qp_host_recv,xp_host_recv,yp_host_recv,zp_host_recv,up_host_recv,vp_host_recv,wp_host_recv,var1p_host_recv
	integer, dimension(:), pinned, allocatable  :: flvp_host_recv,tagp_host_recv
	integer :: np_send_host
	integer :: np_recv_host	
	
	integer, dimension(:,:), allocatable :: lind_host,rind_host,bind_host,tind_host,dind_host,uind_host
	integer, dimension(:), allocatable  :: lc_host,rc_host,bc_host,tc_host,dc_host,uc_host
	

    !prtl arr size, host variables 
	integer, parameter :: Nchunk_prtl_gpu=64 ! total number of chunks
	integer, parameter :: chunk_size_prtl_gpu=gpu_prtl_arr_size/Nchunk_prtl_gpu
	!integer            :: buff_size_prtl_gpu=chunk_size_prtl_gpu
	integer :: np_gpu  !total number of particles on GPU
	integer, dimension(Nchunk_prtl_gpu) :: used_prtl_chunk
	integer            :: empty_prtl_chunk
		
    !derived parameters 	
	integer, parameter :: max_np_gpu=(Nchunk_prtl_gpu-1)*chunk_size_prtl_gpu 
			 
	!parameters used in saving prtl using on CPU 
	integer, dimension(Nchunk_prtl_gpu) :: compact_chunk_offset_cpu		
	integer            :: buff_size_prtl_gpu=chunk_size_prtl_gpu ! size of the arrays that temporarily hold prtls coming/leaving the sub-domain
	
	integer, parameter :: pcount_gpu_size = 64
	integer, dimension(pcount_gpu_size), device :: pcount_gpu
	
	!-----------------------------------------------------------------------------------------------------
	! Test Particles
	!-----------------------------------------------------------------------------------------------------
	
	!test particle array 
    integer, dimension(:), allocatable, device:: flvtp_gpu,tagtp_gpu
    real(psn), dimension(:), allocatable, device:: qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu
    real(psn), dimension(:), allocatable, device:: qtp_send_gpu,xtp_send_gpu,ytp_send_gpu,ztp_send_gpu,utp_send_gpu,vtp_send_gpu,wtp_send_gpu,var1tp_send_gpu
	integer, dimension(:), allocatable, device :: flvtp_send_gpu,tagtp_send_gpu
    real(psn), dimension(:), allocatable, device:: qtp_recv_gpu,xtp_recv_gpu,ytp_recv_gpu,ztp_recv_gpu,utp_recv_gpu,vtp_recv_gpu,wtp_recv_gpu,var1tp_recv_gpu
	integer, dimension(:), allocatable, device :: flvtp_recv_gpu,tagtp_recv_gpu
	integer, device :: ntp_send_gpu
	
	!test prtl arr on host used for communication
    real(psn), dimension(:), pinned, allocatable:: qtp_host_send,xtp_host_send,ytp_host_send,ztp_host_send,utp_host_send,vtp_host_send,wtp_host_send,var1tp_host_send
	integer, dimension(:),pinned, allocatable  :: flvtp_host_send,tagtp_host_send
    real(psn), dimension(:), pinned, allocatable:: qtp_host_recv,xtp_host_recv,ytp_host_recv,ztp_host_recv,utp_host_recv,vtp_host_recv,wtp_host_recv,var1tp_host_recv
	integer, dimension(:), pinned, allocatable  :: flvtp_host_recv,tagtp_host_recv
	integer :: ntp_send_host
	integer :: ntp_recv_host
	
	
	
    !test prtl arr size, host variables 
	integer, parameter :: Nchunk_test_prtl_gpu=8 ! total number of chunks
	integer, parameter :: chunk_size_test_prtl_gpu=gpu_test_prtl_arr_size/Nchunk_test_prtl_gpu
	!integer            :: buff_size_test_prtl_gpu=chunk_size_test_prtl_gpu
	integer :: ntp_gpu  !total number of test particles on GPU
	integer, dimension(Nchunk_test_prtl_gpu) :: used_test_prtl_chunk
	integer            :: empty_test_prtl_chunk
    !derived parameters for test particles
	integer, parameter :: max_ntp_gpu=(Nchunk_test_prtl_gpu-1)*chunk_size_test_prtl_gpu
	
	!parameters used in saving prtl using on CPU 
	integer, dimension(Nchunk_prtl_gpu) :: compact_chunk_offset_test_prtl_cpu
	
    !-------------------------------------------------------------------
    !Curved BC related varaibles 
    !------------------------------------------------------------------- 
	real(psn), dimension(:,:,:), allocatable, device :: e_lx_gpu, e_ly_gpu , e_lz_gpu, b_arx_gpu, b_ary_gpu, b_arz_gpu ! used in conformal FDTD schemes 
	real(psn), dimension(:,:), allocatable, device :: usc_wtr_bx_gpu, usc_wtr_by_gpu, usc_wtr_bz_gpu ! wts used in the USC method, r (c)= row (column) of V matrix
	real(psn), dimension(:,:), allocatable, device :: usc_wtc_bx_gpu, usc_wtc_by_gpu, usc_wtc_bz_gpu
	integer, dimension(:,:,:), allocatable, device :: usc_fdtd_bx_gpu, usc_fdtd_by_gpu, usc_fdtd_bz_gpu ! 0: regular fdtd; positive integer, use the USC method, integer -> index in the wt. list 
	real(psn), dimension(:,:,:), allocatable, device :: usc_db1_gpu, usc_db2_gpu
end module var_gpu