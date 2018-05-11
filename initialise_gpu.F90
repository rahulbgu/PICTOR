module initialise_gpu
	use parameters 
	use vars
	use var_gpu 
	use memory_gpu
	use cudafor
contains 
	subroutine InitAll_gpu 
		integer :: k1,k2
		real(psn), dimension(:), allocatable :: buff_temp_real
		integer, dimension(:), allocatable :: buff_temp_int
		

		if(mod(proc,16).eq.0) then 
			use_device=.true.
		else 
			use_device=.false.
		end if 
#ifdef GPU_EXCLUSIVE 
        use_device=.true.
#endif
        if(.not.use_device) return 
		
		!-----------------------------------------------------------------------------------------------------
		! Boundaries
		!-----------------------------------------------------------------------------------------------------		


	
		xmin_host=xmin+6 !These variables are avialable on host 
		xmax_host=xmax-6 
		ymin_host=ymin+6 
		ymax_host=ymax-6 
		zmin_host=zmin+6
		zmax_host=zmax-6
#ifdef GPU_EXCLUSIVE
		xmin_host=xmin
		xmax_host=xmax 
		ymin_host=ymin 
		ymax_host=ymax 
		zmin_host=zmin
		zmax_host=zmax
#endif 		
		
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
		
		used_prtl_chunk=0
		empty_prtl_chunk=Nchunk_prtl_gpu
		used_test_prtl_chunk=0
		empty_test_prtl_chunk=Nchunk_test_prtl_gpu
!Define Boundaries of GPU domain 


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
!Allocate Fld memory on gpu		
#ifdef twoD
     k1=1 
     k2=1 
#else
     k1=zmin1_host-2
	 k2=zmax1_host+2
#endif 		
		allocate(Ex_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Ey_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Ez_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Bx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(By_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Bz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Jx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Jy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(Jz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(buffJx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(buffJy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(buffJz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))

 
        allocate(TexEx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
        allocate(TexEy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
        allocate(TexEz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
        allocate(TexBx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
        allocate(TexBy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
        allocate(TexBz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))	
#ifdef twoD	
        allocate(VecJ_gpu(8,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
#else			
		allocate(VecJ_gpu(12,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
#endif 	
		
!filtered electric field is used to move test particles
		allocate(FilteredEx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(FilteredEy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(FilteredEz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
			!bffer current on host 
		allocate(buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		!array used in sorting particles 
		allocate(pcount_gpu(Nchunk_prtl_gpu,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
		allocate(pcount_host(Nchunk_prtl_gpu,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2))
				
		!initialise the buffer current memory 
		buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=0.0_psn
		buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=0.0_psn
		buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=0.0_psn
		buffJx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		buffJy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		buffJz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		pcount_host(Nchunk_prtl_gpu,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=0
		pcount_gpu(Nchunk_prtl_gpu,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=pcount_host(Nchunk_prtl_gpu,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
	    FilteredEx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		FilteredEy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		FilteredEz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)	
		
#ifdef twoD			 
		 tBlock_gpu_global = dim3 (16 ,16 ,1)	
		 tGrid_gpu_global = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock_gpu_global%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock_gpu_global%y), 1) 
#else
         tBlock_gpu_global = dim3 (8 ,8 ,4)
         tGrid_gpu_global = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock_gpu_global%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock_gpu_global%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock_gpu_global%z)) 
#endif  

!-----------------------------------------------------------------------------------------------------
! Fld Arrays for transferring data
!-----------------------------------------------------------------------------------------------------		
        mx_gpu_host=xmax1_host-xmin1_host+5 
		my_gpu_host=ymax1_host-ymin1_host+5
		mz_gpu_host=zmax1_host-zmin1_host+5 
#ifdef twoD
        mz_gpu_host=1 
#endif 		 
        mx_gpu=mx_gpu_host 
        my_gpu=my_gpu_host 
        mz_gpu=mz_gpu_host 
		Ncell_gpu=mx_gpu_host*my_gpu_host*mz_gpu_host
		
		
		allocate(lEx_send_gpu(3,my_gpu_host,mz_gpu_host),rEx_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bEx_send_gpu(mx_gpu_host,3,mz_gpu_host),tEx_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dEx_send_gpu(mx_gpu_host,my_gpu_host,3),uEx_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEy_send_gpu(3,my_gpu_host,mz_gpu_host),rEy_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bEy_send_gpu(mx_gpu_host,3,mz_gpu_host),tEy_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dEy_send_gpu(mx_gpu_host,my_gpu_host,3),uEy_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEz_send_gpu(3,my_gpu_host,mz_gpu_host),rEz_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bEz_send_gpu(mx_gpu_host,3,mz_gpu_host),tEz_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dEz_send_gpu(mx_gpu_host,my_gpu_host,3),uEz_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBx_send_gpu(3,my_gpu_host,mz_gpu_host),rBx_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bBx_send_gpu(mx_gpu_host,3,mz_gpu_host),tBx_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dBx_send_gpu(mx_gpu_host,my_gpu_host,3),uBx_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBy_send_gpu(3,my_gpu_host,mz_gpu_host),rBy_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bBy_send_gpu(mx_gpu_host,3,mz_gpu_host),tBy_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dBy_send_gpu(mx_gpu_host,my_gpu_host,3),uBy_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBz_send_gpu(3,my_gpu_host,mz_gpu_host),rBz_send_gpu(2,my_gpu_host,mz_gpu_host))
		allocate(bBz_send_gpu(mx_gpu_host,3,mz_gpu_host),tBz_send_gpu(mx_gpu_host,2,mz_gpu_host))
		allocate(dBz_send_gpu(mx_gpu_host,my_gpu_host,3),uBz_send_gpu(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEx_recv_gpu(2,my_gpu_host,mz_gpu_host),rEx_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bEx_recv_gpu(mx_gpu_host,2,mz_gpu_host),tEx_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dEx_recv_gpu(mx_gpu_host,my_gpu_host,2),uEx_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lEy_recv_gpu(2,my_gpu_host,mz_gpu_host),rEy_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bEy_recv_gpu(mx_gpu_host,2,mz_gpu_host),tEy_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dEy_recv_gpu(mx_gpu_host,my_gpu_host,2),uEy_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lEz_recv_gpu(2,my_gpu_host,mz_gpu_host),rEz_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bEz_recv_gpu(mx_gpu_host,2,mz_gpu_host),tEz_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dEz_recv_gpu(mx_gpu_host,my_gpu_host,2),uEz_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBx_recv_gpu(2,my_gpu_host,mz_gpu_host),rBx_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bBx_recv_gpu(mx_gpu_host,2,mz_gpu_host),tBx_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dBx_recv_gpu(mx_gpu_host,my_gpu_host,2),uBx_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBy_recv_gpu(2,my_gpu_host,mz_gpu_host),rBy_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bBy_recv_gpu(mx_gpu_host,2,mz_gpu_host),tBy_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dBy_recv_gpu(mx_gpu_host,my_gpu_host,2),uBy_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBz_recv_gpu(2,my_gpu_host,mz_gpu_host),rBz_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bBz_recv_gpu(mx_gpu_host,2,mz_gpu_host),tBz_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dBz_recv_gpu(mx_gpu_host,my_gpu_host,2),uBz_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		
		allocate(lEx_send_host(3,my_gpu_host,mz_gpu_host),rEx_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bEx_send_host(mx_gpu_host,3,mz_gpu_host),tEx_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dEx_send_host(mx_gpu_host,my_gpu_host,3),uEx_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEy_send_host(3,my_gpu_host,mz_gpu_host),rEy_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bEy_send_host(mx_gpu_host,3,mz_gpu_host),tEy_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dEy_send_host(mx_gpu_host,my_gpu_host,3),uEy_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEz_send_host(3,my_gpu_host,mz_gpu_host),rEz_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bEz_send_host(mx_gpu_host,3,mz_gpu_host),tEz_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dEz_send_host(mx_gpu_host,my_gpu_host,3),uEz_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBx_send_host(3,my_gpu_host,mz_gpu_host),rBx_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bBx_send_host(mx_gpu_host,3,mz_gpu_host),tBx_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dBx_send_host(mx_gpu_host,my_gpu_host,3),uBx_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBy_send_host(3,my_gpu_host,mz_gpu_host),rBy_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bBy_send_host(mx_gpu_host,3,mz_gpu_host),tBy_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dBy_send_host(mx_gpu_host,my_gpu_host,3),uBy_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lBz_send_host(3,my_gpu_host,mz_gpu_host),rBz_send_host(2,my_gpu_host,mz_gpu_host))
		allocate(bBz_send_host(mx_gpu_host,3,mz_gpu_host),tBz_send_host(mx_gpu_host,2,mz_gpu_host))
		allocate(dBz_send_host(mx_gpu_host,my_gpu_host,3),uBz_send_host(mx_gpu_host,my_gpu_host,2))
		
		allocate(lEx_recv_host(2,my_gpu_host,mz_gpu_host),rEx_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bEx_recv_host(mx_gpu_host,2,mz_gpu_host),tEx_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dEx_recv_host(mx_gpu_host,my_gpu_host,2),uEx_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lEy_recv_host(2,my_gpu_host,mz_gpu_host),rEy_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bEy_recv_host(mx_gpu_host,2,mz_gpu_host),tEy_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dEy_recv_host(mx_gpu_host,my_gpu_host,2),uEy_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lEz_recv_host(2,my_gpu_host,mz_gpu_host),rEz_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bEz_recv_host(mx_gpu_host,2,mz_gpu_host),tEz_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dEz_recv_host(mx_gpu_host,my_gpu_host,2),uEz_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBx_recv_host(2,my_gpu_host,mz_gpu_host),rBx_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bBx_recv_host(mx_gpu_host,2,mz_gpu_host),tBx_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dBx_recv_host(mx_gpu_host,my_gpu_host,2),uBx_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBy_recv_host(2,my_gpu_host,mz_gpu_host),rBy_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bBy_recv_host(mx_gpu_host,2,mz_gpu_host),tBy_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dBy_recv_host(mx_gpu_host,my_gpu_host,2),uBy_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lBz_recv_host(2,my_gpu_host,mz_gpu_host),rBz_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bBz_recv_host(mx_gpu_host,2,mz_gpu_host),tBz_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dBz_recv_host(mx_gpu_host,my_gpu_host,2),uBz_recv_host(mx_gpu_host,my_gpu_host,3))
		
		!Current 
		allocate(lJx_send_gpu(3,my_gpu_host,mz_gpu_host),rJx_send_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJx_send_gpu(mx_gpu_host,3,mz_gpu_host),tJx_send_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJx_send_gpu(mx_gpu_host,my_gpu_host,3),uJx_send_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJy_send_gpu(3,my_gpu_host,mz_gpu_host),rJy_send_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJy_send_gpu(mx_gpu_host,3,mz_gpu_host),tJy_send_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJy_send_gpu(mx_gpu_host,my_gpu_host,3),uJy_send_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJz_send_gpu(3,my_gpu_host,mz_gpu_host),rJz_send_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJz_send_gpu(mx_gpu_host,3,mz_gpu_host),tJz_send_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJz_send_gpu(mx_gpu_host,my_gpu_host,3),uJz_send_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJx_recv_gpu(3,my_gpu_host,mz_gpu_host),rJx_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJx_recv_gpu(mx_gpu_host,3,mz_gpu_host),tJx_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJx_recv_gpu(mx_gpu_host,my_gpu_host,3),uJx_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJy_recv_gpu(3,my_gpu_host,mz_gpu_host),rJy_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJy_recv_gpu(mx_gpu_host,3,mz_gpu_host),tJy_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJy_recv_gpu(mx_gpu_host,my_gpu_host,3),uJy_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJz_recv_gpu(3,my_gpu_host,mz_gpu_host),rJz_recv_gpu(3,my_gpu_host,mz_gpu_host))
		allocate(bJz_recv_gpu(mx_gpu_host,3,mz_gpu_host),tJz_recv_gpu(mx_gpu_host,3,mz_gpu_host))
		allocate(dJz_recv_gpu(mx_gpu_host,my_gpu_host,3),uJz_recv_gpu(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJx_send_host(3,my_gpu_host,mz_gpu_host),rJx_send_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJx_send_host(mx_gpu_host,3,mz_gpu_host),tJx_send_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJx_send_host(mx_gpu_host,my_gpu_host,3),uJx_send_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJy_send_host(3,my_gpu_host,mz_gpu_host),rJy_send_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJy_send_host(mx_gpu_host,3,mz_gpu_host),tJy_send_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJy_send_host(mx_gpu_host,my_gpu_host,3),uJy_send_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJz_send_host(3,my_gpu_host,mz_gpu_host),rJz_send_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJz_send_host(mx_gpu_host,3,mz_gpu_host),tJz_send_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJz_send_host(mx_gpu_host,my_gpu_host,3),uJz_send_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJx_recv_host(3,my_gpu_host,mz_gpu_host),rJx_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJx_recv_host(mx_gpu_host,3,mz_gpu_host),tJx_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJx_recv_host(mx_gpu_host,my_gpu_host,3),uJx_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJy_recv_host(3,my_gpu_host,mz_gpu_host),rJy_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJy_recv_host(mx_gpu_host,3,mz_gpu_host),tJy_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJy_recv_host(mx_gpu_host,my_gpu_host,3),uJy_recv_host(mx_gpu_host,my_gpu_host,3))
		
		allocate(lJz_recv_host(3,my_gpu_host,mz_gpu_host),rJz_recv_host(3,my_gpu_host,mz_gpu_host))
		allocate(bJz_recv_host(mx_gpu_host,3,mz_gpu_host),tJz_recv_host(mx_gpu_host,3,mz_gpu_host))
		allocate(dJz_recv_host(mx_gpu_host,my_gpu_host,3),uJz_recv_host(mx_gpu_host,my_gpu_host,3))	
		
		!used in filter
		allocate(lCurr1_send_gpu(1,my_gpu_host,mz_gpu_host),rCurr1_send_gpu(1,my_gpu_host,mz_gpu_host))
		allocate(bCurr1_send_gpu(mx_gpu_host,1,mz_gpu_host),tCurr1_send_gpu(mx_gpu_host,1,mz_gpu_host))
		allocate(dCurr1_send_gpu(mx_gpu_host,my_gpu_host,1),uCurr1_send_gpu(mx_gpu_host,my_gpu_host,1))
		
		allocate(lCurr1_recv_gpu(1,my_gpu_host,mz_gpu_host),rCurr1_recv_gpu(1,my_gpu_host,mz_gpu_host))
		allocate(bCurr1_recv_gpu(mx_gpu_host,1,mz_gpu_host),tCurr1_recv_gpu(mx_gpu_host,1,mz_gpu_host))
		allocate(dCurr1_recv_gpu(mx_gpu_host,my_gpu_host,1),uCurr1_recv_gpu(mx_gpu_host,my_gpu_host,1))
		
		allocate(lCurr1_send_host(1,my_gpu_host,mz_gpu_host),rCurr1_send_host(1,my_gpu_host,mz_gpu_host))
		allocate(bCurr1_send_host(mx_gpu_host,1,mz_gpu_host),tCurr1_send_host(mx_gpu_host,1,mz_gpu_host))
		allocate(dCurr1_send_host(mx_gpu_host,my_gpu_host,1),uCurr1_send_host(mx_gpu_host,my_gpu_host,1))
		
		allocate(lCurr1_recv_host(1,my_gpu_host,mz_gpu_host),rCurr1_recv_host(1,my_gpu_host,mz_gpu_host))
		allocate(bCurr1_recv_host(mx_gpu_host,1,mz_gpu_host),tCurr1_recv_host(mx_gpu_host,1,mz_gpu_host))
		allocate(dCurr1_recv_host(mx_gpu_host,my_gpu_host,1),uCurr1_recv_host(mx_gpu_host,my_gpu_host,1))
		
		
		
#ifdef twoD			 
		 tBlock_gpu_YZedge = dim3 (4 ,64 ,1)	
		 tGrid_gpu_YZedge = dim3(ceiling(real(3)/tBlock_gpu_YZedge%x), ceiling(real(my_gpu_host)/tBlock_gpu_YZedge%y), 1) 
		 tBlock_gpu_ZXedge = dim3 (64 ,4 ,1)	
		 tGrid_gpu_ZXedge = dim3(ceiling(real(mx_gpu_host)/tBlock_gpu_ZXedge%x), ceiling(real(3)/tBlock_gpu_ZXedge%y), 1) 
		 tBlock_gpu_XYedge = dim3 (1 , 1 ,1)	!not needed
		 tGrid_gpu_XYedge = dim3(ceiling(real(3)/tBlock_gpu_XYedge%x), ceiling(real(3)/tBlock_gpu_XYedge%y), 1) 
#else
         tBlock_gpu_YZedge = dim3 (4 , 1 ,16)
         tGrid_gpu_YZedge = dim3(ceiling(real(3)/tBlock_gpu_YZedge%x), ceiling(real(my_gpu_host)/tBlock_gpu_YZedge%y), ceiling(real(mz_gpu_host)/tBlock_gpu_YZedge%z)) 
         tBlock_gpu_ZXedge = dim3 (64 , 1 , 4)
         tGrid_gpu_ZXedge = dim3(ceiling(real(mx_gpu_host)/tBlock_gpu_ZXedge%x), ceiling(real(3)/tBlock_gpu_ZXedge%y), ceiling(real(mz_gpu_host)/tBlock_gpu_ZXedge%z)) 
         tBlock_gpu_XYedge = dim3 (64 , 4 , 1)
         tGrid_gpu_XYedge = dim3(ceiling(real(mx_gpu_host)/tBlock_gpu_XYedge%x), ceiling(real(my_gpu_host)/tBlock_gpu_XYedge%y), ceiling(real(3)/tBlock_gpu_XYedge%z)) 
#endif 		
		
         allocate(cell_count_gpu(mx_gpu_host*my_gpu_host*mz_gpu_host))
		
		!-----------------------------------------------------------------------------------------------------
		! Particle Arrays
		!-----------------------------------------------------------------------------------------------------			

		!Allocate the buffer memory on the Host to store prtldata recieved from GPU 
		allocate(qp_host_recv(gpu_prtl_arr_size))
		allocate(xp_host_recv(gpu_prtl_arr_size))
		allocate(yp_host_recv(gpu_prtl_arr_size))
		allocate(zp_host_recv(gpu_prtl_arr_size))
		allocate(up_host_recv(gpu_prtl_arr_size))
		allocate(vp_host_recv(gpu_prtl_arr_size))
		allocate(wp_host_recv(gpu_prtl_arr_size))
		allocate(var1p_host_recv(gpu_prtl_arr_size))
		allocate(flvp_host_recv(gpu_prtl_arr_size))
		allocate(tagp_host_recv(gpu_prtl_arr_size))
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
		allocate(qp_host_send(gpu_prtl_arr_size))
		allocate(xp_host_send(gpu_prtl_arr_size))
		allocate(yp_host_send(gpu_prtl_arr_size))
		allocate(zp_host_send(gpu_prtl_arr_size))
		allocate(up_host_send(gpu_prtl_arr_size))
		allocate(vp_host_send(gpu_prtl_arr_size))
		allocate(wp_host_send(gpu_prtl_arr_size))
		allocate(var1p_host_send(gpu_prtl_arr_size))
		allocate(flvp_host_send(gpu_prtl_arr_size))
		allocate(tagp_host_send(gpu_prtl_arr_size))
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
		qp_gpu=qp_host_send
		xp_gpu=xp_host_send
		yp_gpu=yp_host_send
		zp_gpu=zp_host_send
		up_gpu=up_host_send
		vp_gpu=vp_host_send
		wp_gpu=wp_host_send
		var1p_gpu=var1p_host_send
		flvp_gpu=flvp_host_send
		tagp_gpu=tagp_host_send
		
		np_gpu=0 !Host variable		
		
		
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
		buff_tmep_int=0
		
		qp_send_gpu=buff_temp_real
		xp_send_gpu=buff_temp_real
		yp_send_gpu=buff_temp_real
		zp_send_gpu=buff_temp_real
		up_send_gpu=buff_temp_real
		vp_send_gpu=buff_temp_real
		wp_send_gpu=buff_temp_real
		var1p_send_gpu=buff_temp_real
		flvp_send_gpu=buff_tmep_int
		tagp_send_gpu=buff_tmep_int
		
		
		
		
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
		flvp_recv_gpu=buff_tmep_int
		tagp_recv_gpu=buff_tmep_int
		
		deallocate(buff_temp_real)
		deallocate(buff_temp_int)
 

        np_recv_host=0
		np_send_host=0
		np_send_gpu=np_recv_host
		
		!Now transfer particles and flds to GPU
		call UpdateAvailGPUSlots	
        FullCurrentExchange=.true.
		SpillOverSendHost=.true.
		call LoadHostOutlier 
		call SendFullDomainPrtltoGPU
		
		
		!-----------------------------------------------------------------------------------------------------
		! Test Particle Arrays
		!-----------------------------------------------------------------------------------------------------			

		!Allocate the buffer memory on the Host to store prtldata recieved from GPU 
		allocate(qtp_host_recv(gpu_test_prtl_arr_size))
		allocate(xtp_host_recv(gpu_test_prtl_arr_size))
		allocate(ytp_host_recv(gpu_test_prtl_arr_size))
		allocate(ztp_host_recv(gpu_test_prtl_arr_size))
		allocate(utp_host_recv(gpu_test_prtl_arr_size))
		allocate(vtp_host_recv(gpu_test_prtl_arr_size))
		allocate(wtp_host_recv(gpu_test_prtl_arr_size))
		allocate(var1tp_host_recv(gpu_test_prtl_arr_size))
		allocate(flvtp_host_recv(gpu_test_prtl_arr_size))
		allocate(tagtp_host_recv(gpu_test_prtl_arr_size))
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
		allocate(qtp_host_send(gpu_test_prtl_arr_size))
		allocate(xtp_host_send(gpu_test_prtl_arr_size))
		allocate(ytp_host_send(gpu_test_prtl_arr_size))
		allocate(ztp_host_send(gpu_test_prtl_arr_size))
		allocate(utp_host_send(gpu_test_prtl_arr_size))
		allocate(vtp_host_send(gpu_test_prtl_arr_size))
		allocate(wtp_host_send(gpu_test_prtl_arr_size))
		allocate(var1tp_host_send(gpu_test_prtl_arr_size))
		allocate(flvtp_host_send(gpu_test_prtl_arr_size))
		allocate(tagtp_host_send(gpu_test_prtl_arr_size))
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
		
		ntp_gpu=0 !Host variable		
		
		
		allocate(qtp_send_gpu(buff_size_test_prtl_gpu)) !The buffer has smaller size, defined in var_gpu 
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
		buff_tmep_int=0
		
		qtp_send_gpu=buff_temp_real
		xtp_send_gpu=buff_temp_real
		ytp_send_gpu=buff_temp_real
		ztp_send_gpu=buff_temp_real
		utp_send_gpu=buff_temp_real
		vtp_send_gpu=buff_temp_real
		wtp_send_gpu=buff_temp_real
		var1tp_send_gpu=buff_temp_real
		flvtp_send_gpu=buff_tmep_int
		tagtp_send_gpu=buff_tmep_int
		
		
		
		
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
		flvtp_recv_gpu=buff_tmep_int
		tagtp_recv_gpu=buff_tmep_int
		
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
		

		
		
		
		
		!Now transfer particles and flds to GPU
		call UpdateAvailGPUSlotsTestPrtl	
		call LoadHostOutlierTestPrtl 
		call SendFullDomainTestPrtltoGPU
		
		
		!-----------------------------------------------------------------------------------------------------
		! Initialise EM Fld and Current 
		!-----------------------------------------------------------------------------------------------------				
		
		call SendFullDomainEMFldtoGPU
		call ResetCurrentGPU	
				
	end subroutine InitAll_gpu 
	
	
end module initialise_gpu 