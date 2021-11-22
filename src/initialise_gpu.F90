module initialise_gpu
	use parameters 
	use vars
	use var_gpu 
	use particles_gpu
	use fields_gpu
	use cudafor
contains 
	subroutine InitAll_gpu 
		integer :: k1,k2
		real(psn), dimension(:), allocatable :: buff_temp_real
		integer, dimension(:), allocatable :: buff_temp_int
		
		!-----------------------------------------------------------------------------------------------------
		! Boundaries
		!-----------------------------------------------------------------------------------------------------		
		call SetDomainSizeGPU
		
		used_prtl_chunk=0
		empty_prtl_chunk=Nchunk_prtl_gpu
		used_test_prtl_chunk=0
		empty_test_prtl_chunk=Nchunk_test_prtl_gpu
!Define Boundaries of GPU domain 

		call SetThreadBlockGrid		
!-----------------------------------------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------------------------------------	 
!Copy parameters to device constants 		
		Bx_ext0_gpu=Bx_ext0
		By_ext0_gpu=By_ext0
		Bz_ext0_gpu=Bz_ext0
		c_gpu=c
	    fldc_gpu=fldc 
	    fld_halfc_gpu=fld_halfc
	    cinv_gpu=cinv
	    sqc_gpu=sqc
		allocate(flvrqm_gpu(Nflvr))
		flvrqm_gpu=flvrqm
		wtm1_gpu=wtm1
		wt0_gpu=wt0
		wtp1_gpu=wtp1
		!-----------------------------------------------------------------------------------------------------
		! Fields
		!-----------------------------------------------------------------------------------------------------			 		
	
        call AllocateFldGPU 
				
		!initialise the matrices 
		buffJx_host=0.0_psn
		buffJy_host=0.0_psn
		buffJz_host=0.0_psn
		buffJx_gpu=buffJx_host
		buffJy_gpu=buffJy_host
		buffJz_gpu=buffJz_host

	    FilteredEx_gpu=buffJx_host
		FilteredEy_gpu=buffJy_host
		FilteredEz_gpu=buffJz_host	
				
		
		! Initialise EM Fld and Current 
		call SendFullDomainEMFldtoGPU
		call ResetCurrentGPU
		
		!-----------------------------------------------------------------------------------------------------
		! Particle Arrays
		!-----------------------------------------------------------------------------------------------------			

		!Allocate the buffer memory on the Host to store prtldata recieved from GPU 
		allocate(qp_host_recv(buff_size_prtl_gpu))
		allocate(xp_host_recv(buff_size_prtl_gpu))
		allocate(yp_host_recv(buff_size_prtl_gpu))
		allocate(zp_host_recv(buff_size_prtl_gpu))
		allocate(up_host_recv(buff_size_prtl_gpu))
		allocate(vp_host_recv(buff_size_prtl_gpu))
		allocate(wp_host_recv(buff_size_prtl_gpu))
		allocate(var1p_host_recv(buff_size_prtl_gpu))
		allocate(flvp_host_recv(buff_size_prtl_gpu))
		allocate(tagp_host_recv(buff_size_prtl_gpu))
		qp_host_recv=0
		xp_host_recv=0
		yp_host_recv=0
		zp_host_recv=0
		up_host_recv=0
		vp_host_recv=0
		wp_host_recv=0
		var1p_host_recv=0
		flvp_host_recv=0
		tagp_host_recv=0
		!Allocate the buffer memory on the host for prtl arrays to be sent to GPU
		allocate(qp_host_send(buff_size_prtl_gpu))
		allocate(xp_host_send(buff_size_prtl_gpu))
		allocate(yp_host_send(buff_size_prtl_gpu))
		allocate(zp_host_send(buff_size_prtl_gpu))
		allocate(up_host_send(buff_size_prtl_gpu))
		allocate(vp_host_send(buff_size_prtl_gpu))
		allocate(wp_host_send(buff_size_prtl_gpu))
		allocate(var1p_host_send(buff_size_prtl_gpu))
		allocate(flvp_host_send(buff_size_prtl_gpu))
		allocate(tagp_host_send(buff_size_prtl_gpu))
		qp_host_send=0
		xp_host_send=0
		yp_host_send=0
		zp_host_send=0
		up_host_send=0
		vp_host_send=0
		wp_host_send=0
		var1p_host_send=0
		flvp_host_send=0
		tagp_host_send=0
 
        !--------GPU Arrays -----------------------!
		!Allocate prtl arrays on GPU  
		allocate(qp_gpu(gpu_prtl_arr_size))
		allocate(xp_gpu(gpu_prtl_arr_size))
		allocate(yp_gpu(gpu_prtl_arr_size))
		allocate(zp_gpu(gpu_prtl_arr_size))
		allocate(up_gpu(gpu_prtl_arr_size))
		allocate(vp_gpu(gpu_prtl_arr_size))
		allocate(wp_gpu(gpu_prtl_arr_size))
		allocate(var1p_gpu(gpu_prtl_arr_size))
		allocate(flvp_gpu(gpu_prtl_arr_size))
		allocate(tagp_gpu(gpu_prtl_arr_size))
		
		call SendFullDomainPrtltoGPU !copy all data from CPU to GPU	

		
		allocate(qp_send_gpu(buff_size_prtl_gpu)) !The buffer has smaller size, defined in var_gpu 
		allocate(xp_send_gpu(buff_size_prtl_gpu))
		allocate(yp_send_gpu(buff_size_prtl_gpu))
		allocate(zp_send_gpu(buff_size_prtl_gpu))
		allocate(up_send_gpu(buff_size_prtl_gpu))
		allocate(vp_send_gpu(buff_size_prtl_gpu))
		allocate(wp_send_gpu(buff_size_prtl_gpu))
		allocate(var1p_send_gpu(buff_size_prtl_gpu))
		allocate(flvp_send_gpu(buff_size_prtl_gpu))
		allocate(tagp_send_gpu(buff_size_prtl_gpu))
		!initialise arrays
		allocate(buff_temp_real(buff_size_prtl_gpu))
		allocate(buff_temp_int(buff_size_prtl_gpu))
		buff_temp_real=0
		buff_temp_int=0
		
		qp_send_gpu=buff_temp_real
		xp_send_gpu=buff_temp_real
		yp_send_gpu=buff_temp_real
		zp_send_gpu=buff_temp_real
		up_send_gpu=buff_temp_real
		vp_send_gpu=buff_temp_real
		wp_send_gpu=buff_temp_real
		var1p_send_gpu=buff_temp_real
		flvp_send_gpu=buff_temp_int
		tagp_send_gpu=buff_temp_int
		
		
		allocate(qp_recv_gpu(buff_size_prtl_gpu)) !The buffer has smaller size, defined in var_gpu 
		allocate(xp_recv_gpu(buff_size_prtl_gpu))
		allocate(yp_recv_gpu(buff_size_prtl_gpu))
		allocate(zp_recv_gpu(buff_size_prtl_gpu))
		allocate(up_recv_gpu(buff_size_prtl_gpu))
		allocate(vp_recv_gpu(buff_size_prtl_gpu))
		allocate(wp_recv_gpu(buff_size_prtl_gpu))
		allocate(var1p_recv_gpu(buff_size_prtl_gpu))
		allocate(flvp_recv_gpu(buff_size_prtl_gpu))
		allocate(tagp_recv_gpu(buff_size_prtl_gpu)) 
		qp_recv_gpu=buff_temp_real
		xp_recv_gpu=buff_temp_real
		yp_recv_gpu=buff_temp_real
		zp_recv_gpu=buff_temp_real
		up_recv_gpu=buff_temp_real
		vp_recv_gpu=buff_temp_real
		wp_recv_gpu=buff_temp_real
		var1p_recv_gpu=buff_temp_real
		flvp_recv_gpu=buff_temp_int
		tagp_recv_gpu=buff_temp_int
		
		deallocate(buff_temp_real)
		deallocate(buff_temp_int)
 

        np_recv_host=0
		np_send_host=0
		np_send_gpu=np_recv_host
		
		
		!-----------------------------------------------------------------------------------------------------
		! Test Particle Arrays
		!-----------------------------------------------------------------------------------------------------			

		!Allocate the buffer memory on the Host to store prtldata recieved from GPU 
		allocate(qtp_host_recv(buff_size_test_prtl_gpu))
		allocate(xtp_host_recv(buff_size_test_prtl_gpu))
		allocate(ytp_host_recv(buff_size_test_prtl_gpu))
		allocate(ztp_host_recv(buff_size_test_prtl_gpu))
		allocate(utp_host_recv(buff_size_test_prtl_gpu))
		allocate(vtp_host_recv(buff_size_test_prtl_gpu))
		allocate(wtp_host_recv(buff_size_test_prtl_gpu))
		allocate(var1tp_host_recv(buff_size_test_prtl_gpu))
		allocate(flvtp_host_recv(buff_size_test_prtl_gpu))
		allocate(tagtp_host_recv(buff_size_test_prtl_gpu))
		qtp_host_recv=0
		xtp_host_recv=0
		ytp_host_recv=0
		ztp_host_recv=0
		utp_host_recv=0
		vtp_host_recv=0
		wtp_host_recv=0
		var1tp_host_recv=0
		flvtp_host_recv=0
		tagtp_host_recv=0
		!Allocate the buffer memory on the host for prtl arrays to be sent to GPU
		allocate(qtp_host_send(buff_size_test_prtl_gpu))
		allocate(xtp_host_send(buff_size_test_prtl_gpu))
		allocate(ytp_host_send(buff_size_test_prtl_gpu))
		allocate(ztp_host_send(buff_size_test_prtl_gpu))
		allocate(utp_host_send(buff_size_test_prtl_gpu))
		allocate(vtp_host_send(buff_size_test_prtl_gpu))
		allocate(wtp_host_send(buff_size_test_prtl_gpu))
		allocate(var1tp_host_send(buff_size_test_prtl_gpu))
		allocate(flvtp_host_send(buff_size_test_prtl_gpu))
		allocate(tagtp_host_send(buff_size_test_prtl_gpu))
		qtp_host_send=0
		xtp_host_send=0
		ytp_host_send=0
		ztp_host_send=0
		utp_host_send=0
		vtp_host_send=0
		wtp_host_send=0
		var1tp_host_send=0
		flvtp_host_send=0
		tagtp_host_send=0
 
        !--------GPU Arrays -----------------------!
		!Allocate prtl arrays on GPU  
		allocate(qtp_gpu(gpu_test_prtl_arr_size))
		allocate(xtp_gpu(gpu_test_prtl_arr_size))
		allocate(ytp_gpu(gpu_test_prtl_arr_size))
		allocate(ztp_gpu(gpu_test_prtl_arr_size))
		allocate(utp_gpu(gpu_test_prtl_arr_size))
		allocate(vtp_gpu(gpu_test_prtl_arr_size))
		allocate(wtp_gpu(gpu_test_prtl_arr_size))
		allocate(var1tp_gpu(gpu_test_prtl_arr_size))
		allocate(flvtp_gpu(gpu_test_prtl_arr_size))
		allocate(tagtp_gpu(gpu_test_prtl_arr_size))
		qtp_gpu=qtp_host_send
		xtp_gpu=xtp_host_send
		ytp_gpu=ytp_host_send
		ztp_gpu=ztp_host_send
		utp_gpu=utp_host_send
		vtp_gpu=vtp_host_send
		wtp_gpu=wtp_host_send
		var1tp_gpu=var1tp_host_send
		flvtp_gpu=flvtp_host_send
		tagtp_gpu=tagtp_host_send
		
		!Now transfer particles and flds to GPU
		call SendFullDomainTestPrtltoGPU
		
		
		! allocate(qtp_send_gpu(buff_size_test_prtl_gpu)) !The buffer has smaller size, defined in var_gpu
		allocate(xtp_send_gpu(buff_size_test_prtl_gpu))
		allocate(ytp_send_gpu(buff_size_test_prtl_gpu))
		allocate(ztp_send_gpu(buff_size_test_prtl_gpu))
		allocate(utp_send_gpu(buff_size_test_prtl_gpu))
		allocate(vtp_send_gpu(buff_size_test_prtl_gpu))
		allocate(wtp_send_gpu(buff_size_test_prtl_gpu))
		allocate(var1tp_send_gpu(buff_size_test_prtl_gpu))
		allocate(flvtp_send_gpu(buff_size_test_prtl_gpu))
		allocate(tagtp_send_gpu(buff_size_test_prtl_gpu))
		!initialise arrays
		allocate(buff_temp_real(buff_size_test_prtl_gpu))
		allocate(buff_temp_int(buff_size_test_prtl_gpu))
		buff_temp_real=0
		buff_temp_int=0

		qtp_send_gpu=buff_temp_real
		xtp_send_gpu=buff_temp_real
		ytp_send_gpu=buff_temp_real
		ztp_send_gpu=buff_temp_real
		utp_send_gpu=buff_temp_real
		vtp_send_gpu=buff_temp_real
		wtp_send_gpu=buff_temp_real
		var1tp_send_gpu=buff_temp_real
		flvtp_send_gpu=buff_temp_int
		tagtp_send_gpu=buff_temp_int




		allocate(qtp_recv_gpu(buff_size_test_prtl_gpu)) !The buffer has smaller size, defined in var_gpu
		allocate(xtp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(ytp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(ztp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(utp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(vtp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(wtp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(var1tp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(flvtp_recv_gpu(buff_size_test_prtl_gpu))
		allocate(tagtp_recv_gpu(buff_size_test_prtl_gpu))
		qtp_recv_gpu=buff_temp_real
		xtp_recv_gpu=buff_temp_real
		ytp_recv_gpu=buff_temp_real
		ztp_recv_gpu=buff_temp_real
		utp_recv_gpu=buff_temp_real
		vtp_recv_gpu=buff_temp_real
		wtp_recv_gpu=buff_temp_real
		var1tp_recv_gpu=buff_temp_real
		flvtp_recv_gpu=buff_temp_int
		tagtp_recv_gpu=buff_temp_int

		deallocate(buff_temp_real)
		deallocate(buff_temp_int)


        ntp_recv_host=0
		ntp_send_host=0
		ntp_send_gpu=ntp_recv_host
		
		
		!varaibles used in loading particle outliers
		allocate(lind_host(buff_size_prtl_gpu,Nthreads))
	    allocate(rind_host(buff_size_prtl_gpu,Nthreads))
	    allocate(bind_host(buff_size_prtl_gpu,Nthreads))
		allocate(tind_host(buff_size_prtl_gpu,Nthreads))
		allocate(dind_host(buff_size_prtl_gpu,Nthreads))
		allocate(uind_host(buff_size_prtl_gpu,Nthreads))
		
		allocate(lc_host(Nthreads))
		allocate(rc_host(Nthreads))
		allocate(bc_host(Nthreads))
		allocate(tc_host(Nthreads))
		allocate(dc_host(Nthreads))
		allocate(uc_host(Nthreads))
		
				
	end subroutine InitAll_gpu 
	
	subroutine SetDomainSizeGPU 
		mx_gpu=mx
		my_gpu=my
		mz_gpu=mz

		xmin_host=xmin
		xmax_host=xmax 
		ymin_host=ymin 
		ymax_host=ymax 
		zmin_host=zmin
		zmax_host=zmax	
		
#ifdef twoD
		zmin_host=1.0_psn
		zmax_host=2.0_psn
#endif 		
		
		xmin_gpu=xmin_host !varaibles on GPU 
		xmax_gpu=xmax_host
		ymin_gpu=ymin_host
		ymax_gpu=ymax_host
		zmin_gpu=zmin_host
		zmax_gpu=zmax_host
		
		xmin1_host=xmin_host
		xmax1_host=xmax_host 
		ymin1_host=ymin_host
		ymax1_host=ymax_host 
		zmin1_host=zmin_host
		zmax1_host=zmax_host 
			
		xmin1_gpu=xmin1_host
		xmax1_gpu=xmax1_host 
		ymin1_gpu=ymin1_host
		ymax1_gpu=ymax1_host 
		zmin1_gpu=zmin1_host
		zmax1_gpu=zmax1_host 
        xmin0_host=xmin_host!because xmin_host is allowed to change in load balancing 	
		
	end subroutine SetDomainSizeGPU
	
	
	
	subroutine SetThreadBlockGrid
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(mx)/tBlock%x), ceiling(real(my)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(mx)/tBlock%x), ceiling(real(my)/tBlock%y), ceiling(real(mz)/tBlock%z)) 
#endif  

#ifdef twoD			 
		 tBlock_gpu_YZedge = dim3 (4 ,64 ,1)	
		 tGrid_gpu_YZedge = dim3(ceiling(real(3)/tBlock_gpu_YZedge%x), ceiling(real(my)/tBlock_gpu_YZedge%y), 1) 
		 tBlock_gpu_ZXedge = dim3 (64 ,4 ,1)	
		 tGrid_gpu_ZXedge = dim3(ceiling(real(mx)/tBlock_gpu_ZXedge%x), ceiling(real(3)/tBlock_gpu_ZXedge%y), 1) 
		 tBlock_gpu_XYedge = dim3 (1 , 1 ,1)	!not needed
		 tGrid_gpu_XYedge = dim3(ceiling(real(3)/tBlock_gpu_XYedge%x), ceiling(real(3)/tBlock_gpu_XYedge%y), 1) 
		 
		 tBlock_gpu_YZedge1 = dim3 (1 , 256 ,1)	
		 tGrid_gpu_YZedge1 = dim3(1, ceiling(real(my)/tBlock_gpu_YZedge1%y), 1) 
#else
         tBlock_gpu_YZedge = dim3 (4 , 1 ,16)
         tGrid_gpu_YZedge = dim3(ceiling(real(3)/tBlock_gpu_YZedge%x), ceiling(real(my)/tBlock_gpu_YZedge%y), ceiling(real(mz)/tBlock_gpu_YZedge%z)) 
         tBlock_gpu_ZXedge = dim3 (64 , 1 , 4)
         tGrid_gpu_ZXedge = dim3(ceiling(real(mx)/tBlock_gpu_ZXedge%x), ceiling(real(3)/tBlock_gpu_ZXedge%y), ceiling(real(mz)/tBlock_gpu_ZXedge%z)) 
         tBlock_gpu_XYedge = dim3 (64 , 4 , 1)
         tGrid_gpu_XYedge = dim3(ceiling(real(mx)/tBlock_gpu_XYedge%x), ceiling(real(my)/tBlock_gpu_XYedge%y), ceiling(real(3)/tBlock_gpu_XYedge%z)) 

		 tBlock_gpu_YZedge1 = dim3 (1 , 64 ,4)	
		 tGrid_gpu_YZedge1 = dim3(1, ceiling(real(my)/tBlock_gpu_YZedge1%y), ceiling(real(mz)/tBlock_gpu_YZedge1%z)) 
#endif 	

	end subroutine SetThreadBlockGrid 
	
    subroutine DeallocateFldGPU 
		deallocate(cell_count_gpu)
        deallocate(Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,Jx_gpu,Jy_gpu,Jz_gpu)
		deallocate(buffJx_gpu,buffJy_gpu,buffJz_gpu,buffJx_host,buffJy_host,buffJz_host)
		
		deallocate(Jx_wide_gpu,Jy_wide_gpu,Jz_wide_gpu)
		
        deallocate(TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu)
		deallocate(FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu)
	end subroutine DeallocateFldGPU
	subroutine AllocateFldGPU
		Ncell_gpu=mx*my*mz
		allocate(cell_count_gpu(mx*my*mz))
        allocate(Ex_gpu(mx,my,mz),Ey_gpu(mx,my,mz),Ez_gpu(mx,my,mz))
        allocate(Bx_gpu(mx,my,mz),By_gpu(mx,my,mz),Bz_gpu(mx,my,mz))
        allocate(Jx_gpu(mx,my,mz),Jy_gpu(mx,my,mz),Jz_gpu(mx,my,mz))
		
		allocate(Jx_wide_gpu(mx,my,mz,Jwidth_gpu),Jy_wide_gpu(mx,my,mz,Jwidth_gpu),Jz_wide_gpu(mx,my,mz,Jwidth_gpu))
		
		allocate(buffJx_gpu(mx,my,mz),buffJy_gpu(mx,my,mz),buffJz_gpu(mx,my,mz))
		allocate(buffJx_host(mx,my,mz),buffJy_host(mx,my,mz),buffJz_host(mx,my,mz))
		
        allocate(TexEx_gpu(mx,my,mz),TexEy_gpu(mx,my,mz),TexEz_gpu(mx,my,mz))
        allocate(TexBx_gpu(mx,my,mz),TexBy_gpu(mx,my,mz),TexBz_gpu(mx,my,mz))
		allocate(FilteredEx_gpu(mx,my,mz),FilteredEy_gpu(mx,my,mz),FilteredEz_gpu(mx,my,mz))
		
	end subroutine AllocateFldGPU
		
	
end module initialise_gpu 