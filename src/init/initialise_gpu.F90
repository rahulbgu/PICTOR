module initialise_gpu
	use parameters 
	use vars
	use var_gpu 
	use particles_gpu
	use fields_gpu
	use cudafor
contains 
	subroutine InitAll_gpu 
		integer, dimension(:), allocatable :: p_int
		real(psn), dimension(:), allocatable :: p_real		
		!-----------------------------------------------------------------------------------------------------
		! Boundaries and Threads
		!-----------------------------------------------------------------------------------------------------		
		call SetDomainSizeGPU

		call SetThreadBlockGrid				

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

		if(nMoverEMfilter.gt.0) then 
		    FilteredEx_gpu=buffJx_host
			FilteredEy_gpu=buffJy_host
			FilteredEz_gpu=buffJz_host	
		end if
		
		! Initialise EM Fld and Current 
		call SendFullDomainEMFldtoGPU
		call ResetCurrentGPU
		
		!-----------------------------------------------------------------------------------------------------
		! Particle Arrays
		!-----------------------------------------------------------------------------------------------------
		allocate(flvrqm_gpu(Nflvr))
		flvrqm_gpu=flvrqm
		
		used_prtl_chunk=0
		empty_prtl_chunk=Nchunk_prtl_gpu
		
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
		
		!initialise prtl arrays
		allocate(p_real(gpu_prtl_arr_size))
		p_real = 0 
		qp_gpu = p_real
		xp_gpu = p_real
		yp_gpu = p_real
		zp_gpu = p_real
		up_gpu = p_real
		vp_gpu = p_real
		wp_gpu = p_real
        var1p_gpu = p_real
		deallocate(p_real)
		allocate(p_int(gpu_prtl_arr_size))
		p_int = 0 
		flvp_gpu = p_int
		tagp_gpu = p_int
		deallocate(p_int)
		
		allocate(qp_host(chunk_size_prtl_gpu))
		allocate(xp_host(chunk_size_prtl_gpu))
		allocate(yp_host(chunk_size_prtl_gpu))
		allocate(zp_host(chunk_size_prtl_gpu))
		allocate(up_host(chunk_size_prtl_gpu))
		allocate(vp_host(chunk_size_prtl_gpu))
		allocate(wp_host(chunk_size_prtl_gpu))
		allocate(var1p_host(chunk_size_prtl_gpu))
		allocate(flvp_host(chunk_size_prtl_gpu))
		allocate(tagp_host(chunk_size_prtl_gpu))
		xp_host =0; yp_host=0; zp_host=0; 
		up_host =0; vp_host=0; wp_host=0;
		qp_host =0; flvp_host=0; tagp_host=0; var1p=0;  
		
		
		call SendFullDomainPrtltoGPU !copy all data from CPU to GPU	 

        np_recv_host=0
		np_send_host=0
		np_send_gpu=np_recv_host						
	end subroutine InitAll_gpu 
	
	subroutine SetDomainSizeGPU 
		mx_gpu=mx
		my_gpu=my
		mz_gpu=mz		
	end subroutine SetDomainSizeGPU
	
	subroutine ResetGridGPU
        call SetDomainSizeGPU
   	    call DeallocateFldGPU
   	    call AllocateFldGPU

        call SendFullDomainEMFldToGPU
   	    call SendFullDomainCurrToGPU
   	    call SetThreadBlockGrid
	end subroutine ResetGridGPU
	
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
		if(nMoverEMfilter.gt.0) deallocate(FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu)
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
		if(nMoverEMfilter.gt.0) allocate(FilteredEx_gpu(mx,my,mz),FilteredEy_gpu(mx,my,mz),FilteredEz_gpu(mx,my,mz))
		
	end subroutine AllocateFldGPU
		
	
end module initialise_gpu 