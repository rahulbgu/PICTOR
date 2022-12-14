module particles_gpu 
	use parameters
	use vars 
	use var_gpu 
	use mem_prtl
	use cudafor
	use communication
	use particles
	!use thrust_module
	use cudadevice
#ifdef cyl
    use cyl_comm_fldprtl
#endif 	
    implicit none
contains 
	
	subroutine SendFullDomainPrtltoGPU
		!copy all particles, use all chunks
		call SendPrtlToGPU(Nchunk_prtl_gpu)
	
		if(used_prtl_arr_size.ne.0) call Abort(12) ! some particles could not be copied; not enough space on GPU
	end subroutine SendFullDomainPrtltoGPU
	
	
	!----------------------------------------------------------------------------
	!  Copy particle data to/from Host/GPU  
	!----------------------------------------------------------------------------
	subroutine CopyPrtlToGPU(np_copy,off_host,off_gpu)! from CPU main prtl array to GPU 	
		integer :: np_copy, off_host, off_gpu
		
		if(np_copy.eq.0) return
		flvp_gpu(off_gpu+1:off_gpu+np_copy) = flvp(off_host+1:off_host+np_copy)
		tagp_gpu(off_gpu+1:off_gpu+np_copy) = tagp(off_host+1:off_host+np_copy)
	 	qp_gpu(off_gpu+1:off_gpu+np_copy) = qp(off_host+1:off_host+np_copy)
		xp_gpu(off_gpu+1:off_gpu+np_copy) = xp(off_host+1:off_host+np_copy)
		yp_gpu(off_gpu+1:off_gpu+np_copy) = yp(off_host+1:off_host+np_copy)
		zp_gpu(off_gpu+1:off_gpu+np_copy) = zp(off_host+1:off_host+np_copy)
		up_gpu(off_gpu+1:off_gpu+np_copy) = up(off_host+1:off_host+np_copy)
		vp_gpu(off_gpu+1:off_gpu+np_copy) = vp(off_host+1:off_host+np_copy)
		wp_gpu(off_gpu+1:off_gpu+np_copy) = wp(off_host+1:off_host+np_copy)
		var1p_gpu(off_gpu+1:off_gpu+np_copy) = var1p(off_host+1:off_host+np_copy)
	end subroutine CopyPrtlToGPU
	
	subroutine CopyPrtlToHost(np_copy,off_host,off_gpu)! from GPU prtl array to the main CPU prtl array
		integer :: np_copy, off_host, off_gpu
		if(np_copy.eq.0) return
		flvp(off_host+1:off_host+np_copy) = flvp_gpu(off_gpu+1:off_gpu+np_copy)
		tagp(off_host+1:off_host+np_copy) = tagp_gpu(off_gpu+1:off_gpu+np_copy)
	 	qp(off_host+1:off_host+np_copy) = qp_gpu(off_gpu+1:off_gpu+np_copy)
		xp(off_host+1:off_host+np_copy) = xp_gpu(off_gpu+1:off_gpu+np_copy)
		yp(off_host+1:off_host+np_copy) = yp_gpu(off_gpu+1:off_gpu+np_copy)
		zp(off_host+1:off_host+np_copy) = zp_gpu(off_gpu+1:off_gpu+np_copy)
		up(off_host+1:off_host+np_copy) = up_gpu(off_gpu+1:off_gpu+np_copy)
		vp(off_host+1:off_host+np_copy) = vp_gpu(off_gpu+1:off_gpu+np_copy)
		wp(off_host+1:off_host+np_copy) = wp_gpu(off_gpu+1:off_gpu+np_copy)
		var1p(off_host+1:off_host+np_copy) = var1p_gpu(off_gpu+1:off_gpu+np_copy)
	end subroutine CopyPrtlToHost
	
	subroutine CopyPrtlToHostBuff(np_copy,off_host,off_gpu)! from GPU prtl array to the main CPU prtl array
		integer :: np_copy, off_host, off_gpu
		if(np_copy.eq.0) return
		flvp_host(off_host+1:off_host+np_copy) = flvp_gpu(off_gpu+1:off_gpu+np_copy)
		tagp_host(off_host+1:off_host+np_copy) = tagp_gpu(off_gpu+1:off_gpu+np_copy)
	 	qp_host(off_host+1:off_host+np_copy) = qp_gpu(off_gpu+1:off_gpu+np_copy)
		xp_host(off_host+1:off_host+np_copy) = xp_gpu(off_gpu+1:off_gpu+np_copy)
		yp_host(off_host+1:off_host+np_copy) = yp_gpu(off_gpu+1:off_gpu+np_copy)
		zp_host(off_host+1:off_host+np_copy) = zp_gpu(off_gpu+1:off_gpu+np_copy)
		up_host(off_host+1:off_host+np_copy) = up_gpu(off_gpu+1:off_gpu+np_copy)
		vp_host(off_host+1:off_host+np_copy) = vp_gpu(off_gpu+1:off_gpu+np_copy)
		wp_host(off_host+1:off_host+np_copy) = wp_gpu(off_gpu+1:off_gpu+np_copy)
		var1p_host(off_host+1:off_host+np_copy) = var1p_gpu(off_gpu+1:off_gpu+np_copy)
	end subroutine CopyPrtlToHostBuff
	
	!----------------------------------------------------------------------------
	!  Send particle data from CPU (main array) to GPU (main array)  
	!----------------------------------------------------------------------------
	subroutine SendPrtlToGPU(kc_max)
		integer :: kc_max!max index of the chunck allowed to be filled in with incoming particles
		integer :: kc, np_copy
		do kc = 1, kc_max
			if(kc.eq.empty_prtl_chunk) cycle !don't insert any particle in the empty chunk
			
			np_copy = min(chunk_size_prtl_gpu-used_prtl_chunk(kc),used_prtl_arr_size)
			
			if(np_copy.ne.0) then 
				
				call CopyPrtlToGPU(np_copy,used_prtl_arr_size-np_copy,used_prtl_chunk(kc))
				
				used_prtl_arr_size = used_prtl_arr_size - np_copy
				qp(used_prtl_arr_size+1:used_prtl_arr_size+np_copy) = 0 
			
			    used_prtl_chunk(kc) = used_prtl_chunk(kc) + np_copy
			    np_gpu = np_gpu + np_copy
		    end if 
			
			if(used_prtl_arr_size.eq.0) exit
		end do 
	end subroutine SendPrtlToGPU	
	
	!----------------------------------------------------------------------------
	!  a short array of parallel counters used on gpus 
	!----------------------------------------------------------------------------
	subroutine Reset_pcount_gpu
		integer, dimension(pcount_gpu_size) :: count
		count = 0 
		pcount_gpu = count
	end subroutine Reset_pcount_gpu
	
	integer function gather_pcount_gpu()
		integer, dimension(pcount_gpu_size) :: count_cpu 
		integer :: n
		count_cpu = pcount_gpu 
		gather_pcount_gpu = 0 
		do n = 1, pcount_gpu_size
			gather_pcount_gpu = gather_pcount_gpu + count_cpu(n)
		end do
	end function gather_pcount_gpu

	!----------------------------------------------------------------------------
	!  Exchange particles with neighbouring subdomains
	!----------------------------------------------------------------------------
	
	subroutine ExchangePrtlGPU
		integer :: n !used in temporary fix 
		call StartTimer(32)
		call LoadGPUOutlier
		call StopTimer(32)
		if(proc.eq.0) print*,'total time in loading the outliers', real(exec_time(32))
		
   		!GPU To CPU 
		call StartTimer(33)  
   		if(np_recv_host.gt.0) then 
			call CopyPrtlToHostBuff(np_recv_host,0,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
   	    end if 
   		np_gpu=np_gpu-np_recv_host
		
		call StopTimer(33)

		call ExchangePrtl(qp_host,xp_host,yp_host,zp_host,up_host,vp_host,wp_host,var1p_host,flvp_host,tagp_host,np_recv_host)
		call AppendParticles   

 		call StopTimer(35)
 		 ! print*,'total time MPI comm of Prtl data', real(exec_time(35))	  
		  
        call StartTimer(36)
		  
		call SendPrtlToGPU(Nchunk_prtl_gpu)
		    
		call StopTimer(36) 

	end subroutine ExchangePrtlGPU
	


!     subroutine SendPrtlToGPU_OLD(size,x,y,z,u,v,w,q,tag,flv,var,buff_size,ind_max)
! 		integer :: size
! 		real(psn), dimension(size) :: x,y,z,u,v,w,var,q
! 		integer, dimension(size) :: tag,flv
! 		integer :: ind_max, buff_size
! 		integer :: off, ind
!
! 		do off=0,ind_max-1,buff_size
! 			ind=min(buff_size,ind_max-off)
!
! 			qp_recv_gpu(1:ind)=q(off+1:off+ind)
!   			xp_recv_gpu(1:ind)=x(off+1:off+ind)
!   			yp_recv_gpu(1:ind)=y(off+1:off+ind)
!   			zp_recv_gpu(1:ind)=z(off+1:off+ind)
!   			up_recv_gpu(1:ind)=u(off+1:off+ind)
!   			vp_recv_gpu(1:ind)=v(off+1:off+ind)
!   			wp_recv_gpu(1:ind)=w(off+1:off+ind)
!   			var1p_recv_gpu(1:ind)=var(off+1:off+ind)
!   			flvp_recv_gpu(1:ind)=flv(off+1:off+ind)
!   			tagp_recv_gpu(1:ind)=tag(off+1:off+ind)
!
! 			call AppendGPUPrtlArr(ind)
! 			np_gpu=np_gpu+ind
! 		end do
! 	end subroutine SendPrtlToGPU_OLD
!
!
! 	subroutine AppendGPUPrtlArr(np_append)
! 		integer :: np_append
! 		integer :: kc !Kernel Call
! 		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted
!
!         if(np_append.eq.0) return
! 		Inserted=0
! 		do kc=1,Nchunk_prtl_gpu
! 			if(kc.eq.empty_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
! 			indi_prtl=used_prtl_chunk(kc)+1+(kc-1)*chunk_size_prtl_gpu
! 			indf_prtl=chunk_size_prtl_gpu+(kc-1)*chunk_size_prtl_gpu
! 			indi_buff=Inserted+1
! 			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),np_append)
! 				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(gpu_prtl_arr_size,buff_size_prtl_gpu,qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
! 				qp_recv_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,var1p_recv_gpu,flvp_recv_gpu,tagp_recv_gpu,&
! 				indi_prtl,indf_prtl,indi_buff,indf_buff)
!
! 		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles
! 		    used_prtl_chunk(kc)=used_prtl_chunk(kc)+(indf_buff-indi_buff)+1
! 			if(Inserted.ge.np_append) exit ! no need to continue
!  		end do
! 	end subroutine AppendGPUPrtlArr
!
! 	attributes(global) subroutine  AppendGPUPrtlArrKernel(size1,size2,q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,i1,j1,i2,j2)
! 	     integer, value :: size1, size2
! 		 real, dimension(size1) ::q1,x1,y1,z1,u1,v1,w1,var1
! 		 integer, dimension(size1) :: flv1,tag1
! 	     real, dimension(size2) ::q2,x2,y2,z2,u2,v2,w2,var2
! 		 integer, dimension(size2) :: flv2,tag2
! 		 integer, value :: i1,i2,j1,j2
! 		 integer :: n1,n2
!
! 		 n1 = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
! 		 n2 = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i2-1)
!
! 		 if((n1.le.j1).and.(n2.le.j2)) then
! 			 q1(n1)=q2(n2)
! 		     x1(n1)=x2(n2)
! 			 y1(n1)=y2(n2)
! 			 z1(n1)=z2(n2)
! 			 u1(n1)=u2(n2)
! 			 v1(n1)=v2(n2)
! 			 w1(n1)=w2(n2)
! 			 var1(n1)=var2(n2)
! 			 flv1(n1)=flv2(n2)
! 			 tag1(n1)=tag2(n2)
! 		 end if
!
! 	end subroutine  AppendGPUPrtlArrKernel
	
	
	!-----------------------------------------------------------------------------------
	!
	!   Load prtl data into the buffer memory if the particle has moved outside all boundaries
	!
	!-----------------------------------------------------------------------------------

	subroutine LoadGPUOutlier
		integer :: kc, indi, indf
	
		np_recv_host=0 
		np_send_gpu=np_recv_host

		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call LoadGPUOutlierKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(chunk_size_prtl_gpu,qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		                                                                                               xmin,xmax,ymin,ymax,zmin,zmax,np_send_gpu,indi,indf,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
	    end do 
		
		np_recv_host=np_send_gpu !update the incoming prtl count on the host 	
	end subroutine LoadGPUOutlier	


	attributes(global) subroutine LoadGPUOutlierKernel(size,q,x,y,z,u,v,w,var,flv,tag,xmin,xmax,ymin,ymax,zmin,zmax,np_send_gpu,i1,i2,off)
	     integer, value :: size
		 real, dimension(size) ::q,x,y,z,u,v,w,var
		 integer, dimension(size) :: flv,tag   
		 real, value  :: xmin,xmax,ymin,ymax,zmin,zmax 
		 integer :: np_send_gpu
		 integer, value :: off
		 integer, value :: i1,i2
		 integer :: InsertAt
		 integer :: n
		 
		 
		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
		 
		 if(n.gt.i2) return !to avoid out of bound memory access 
		 if((q(n).ne.0)) then	

#ifndef twoD			 
         if((x(n).lt.xmin).or.(x(n).gt.xmax).or.(y(n).lt.ymin).or.(y(n).gt.ymax).or.(z(n).lt.zmin).or.(z(n).gt.zmax)) then
#else  					 
         if((x(n).lt.xmin).or.(x(n).gt.xmax).or.(y(n).lt.ymin).or.(y(n).gt.ymax)) then
#endif 						 
						 InsertAt=atomicinc(np_send_gpu,400000000)
						 InsertAt=InsertAt+1 + off
												  
						 q(InsertAt)=q(n)		 
					     x(InsertAt)=x(n)
						 y(InsertAt)=y(n)
						 z(InsertAt)=z(n)
						 u(InsertAt)=u(n)
						 v(InsertAt)=v(n)
						 w(InsertAt)=w(n)
						 var(InsertAt)=var(n)
						 flv(InsertAt)=flv(n)
						 tag(InsertAt)=tag(n)
						 
						 ! delete the particle
						 q(n)=0
						 u(n)=0
						 v(n)=0
						 w(n)=0
						 tag(n)=0
						 flv(n)=0

		 end if 			  
		 end if 
		 	
	end subroutine LoadGPUOutlierKernel

	
! #ifdef GPU_USE_INTRINSICS
! 	attributes(global) subroutine LoadGPUOutlierKernel(q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,xmin1,xmax1,ymin1,ymax1,zmin1,zmax1,np_send_gpu,Nmax,i1,i2)
! 	     real, dimension(:) ::q1,x1,y1,z1,u1,v1,w1,var1
! 		 integer, dimension(:) :: flv1,tag1
! 	     real, dimension(:) ::q2,x2,y2,z2,u2,v2,w2,var2
! 		 integer, dimension(:) :: flv2,tag2
! 		 real, value  :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1
! 		 integer :: np_send_gpu
! 		 integer, value  :: Nmax
! 		 integer, value :: i1,i2
! 		 integer :: InsertAt
! 		 integer :: n
! 		 integer :: res
! 		 integer :: res_pred
! 		 integer :: WarpID
! 		 integer :: warp_count
! 		 integer :: warp_pred
!
! 		 !integer :: ALL_TRUE
! 		 !DATA ALL_TRUE /B'11111111111111111111111111111111'/
!
!
! 		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
!
! 		 if(n.gt.i2) return !to avoid out of bound memory access
!
!
! #ifndef twoD
!          if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1).or.(z1(n).lt.zmin1).or.(z1(n).gt.zmax1)) then
! #else
!          if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1)) then
! #endif
!              res=1
! 		 else
! 			 res=0
! 		 end if
! 		 if((q1(n).eq.0)) res=0
!
! 		 warp_pred=ballot(res)
! 		 !count=syncthreads_count(res)
! 		 warp_count=__popc(warp_pred)
!
! 		 if(warp_count.ne.0) then
! 			 WarpID=mod((threadIdx%x-1), 32) + 1 !WARP_SIZE is assumed to be 32
! 			 if(WarpID.eq.1) InsertAt=atomicadd(np_send_gpu,warp_count)
! 			 InsertAt = __shfl(InsertAt, 1)
!
! 			 res_pred=IAND(warp_pred,IBITS(B'11111111111111111111111111111111',warpID,warpID))
! 			 warp_count=__popc(res_pred)
! 			 InsertAt=InsertAt+warp_count
! 			         if(res.eq.1) then
! 						 q2(InsertAt)=q1(n)
! 					     x2(InsertAt)=x1(n)
! 						 y2(InsertAt)=y1(n)
! 						 z2(InsertAt)=z1(n)
! 						 u2(InsertAt)=u1(n)
! 						 v2(InsertAt)=v1(n)
! 						 w2(InsertAt)=w1(n)
! 						 var2(InsertAt)=var1(n)
! 						 flv2(InsertAt)=flv1(n)
! 						 tag2(InsertAt)=tag1(n)
!
! 						 ! delete the particle
! 						 q1(n)=0
! 						 x1(n)=xmin1+0.7
! 						 y1(n)=ymin1+0.7
! 						 z1(n)=zmin1+0.7
! 						 u1(n)=0
! 						 v1(n)=0
! 						 w1(n)=0
! 						 var1(n)=0
! 						 tag1(n)=0
! 						 flv1(n)=0
! 					 end if
! 		    end if
! 	end subroutine LoadGPUOutlierKernel
! #endif	

!-----------------------------------------------------------------------------------
!
!   Load prtl data into the buffer memory if the particle x-position is larger than 'xb' 
!
!-----------------------------------------------------------------------------------

	subroutine LoadGPUOutlierRight(xb)
		integer :: kc, indi, indf
		real(psn) :: xb 
		real(psn) :: inf = 2147483647.0
	
		np_recv_host=0 
		np_send_gpu=np_recv_host

		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call LoadGPUOutlierKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(chunk_size_prtl_gpu,qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		                                                                                               -inf,xb,-inf,inf,-inf,inf,np_send_gpu,indi,indf,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
	    end do 
	
		np_recv_host=np_send_gpu !update the incoming prtl count on the host 	
		call CopyPrtlToHost(np_recv_host,used_prtl_arr_size,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
		used_prtl_arr_size = used_prtl_arr_size + np_recv_host
	end subroutine LoadGPUOutlierRight	

	
!--------------------------------------------------------------------------------------------
!
!   Particle Sorting 
!
!--------------------------------------------------------------------------------------------	
subroutine ReorderPrtlGPU
	 if(modulo(t,prtl_reorder_period).ne.0) return 
	 !call ReorderPrtlArrGPU
end subroutine ReorderPrtlGPU	

subroutine ReorderPrtlArrGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: total_pcount_chunk
	integer, dimension(Ncell_gpu) :: count_temp
	
	
	kc=empty_prtl_chunk
	
	do nn=1,Nchunk_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
		
		   indi=(kc-1)*chunk_size_prtl_gpu+1
		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(cell_count_gpu,Ncell_gpu)
		   
		   call BucketCountPrtl<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,mx,my)
		 	   
		   !call thrustscan(cell_count_gpu,Ncell_gpu)
		   call cpu_scan_temp(Ncell_gpu)
		   total_pcount_chunk=cell_count_gpu(Ncell_gpu)
		   
		   call BucketSortCopyPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,empty_prtl_chunk,chunk_size_prtl_gpu,mx,my)
		   
		used_prtl_chunk(empty_prtl_chunk)=total_pcount_chunk
		
		used_prtl_chunk(kc)=0
	    empty_prtl_chunk=kc
	end do 
	!print*,'After Sorting used prtl chunk',used_prtl_chunk(:)
	!if(proc.eq.0) print*,'Total Prtl After Sorting',sum(used_prtl_chunk)
	
end subroutine ReorderPrtlArrGPU  

	
subroutine cpu_scan_temp(ncell)
	integer :: ncell
	integer, dimension(ncell) :: count_cpu
	integer :: i
	count_cpu=cell_count_gpu
	do i=2,ncell
		count_cpu(i)=count_cpu(i)+count_cpu(i-1)
	end do 
	cell_count_gpu=count_cpu
end subroutine cpu_scan_temp


attributes(global) subroutine ResetCellCount(Arr,N)
         integer, value :: N 
		 integer, dimension(N) :: Arr
		 integer :: i
		 i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
		 if(i.le.N) Arr(i)=0
end subroutine ResetCellCount  

attributes(global) subroutine BucketCountPrtl(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,mx_gpu,my_gpu)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
    integer, value :: i1,i2 
    integer, value :: arr_size 
    integer, dimension(arr_size) :: Arr 
	integer, value :: mx_gpu,my_gpu
	
	integer :: i,j,k		
	integer :: stat
	integer :: ind
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(flv(n).ne.0) then 
			i=floor(x(n))
			j=floor(y(n))
			k=floor(z(n))
			
! 			i=floor((x(n)-3)/8)+3
! 			j=floor((y(n)-3)/8)+3
! 			k=floor((z(n)-3)/8)+3
#ifdef twoD
            k=1  
#endif			
            ind=my_gpu*mx_gpu*(k-1)+mx_gpu*(j-1)+i 
			stat=atomicadd(Arr(ind),1)
		end if 
	end if  
end subroutine BucketCountPrtl

attributes(global) subroutine BucketSortCopyPrtlKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,chunk_dest,chunk_size,mx_gpu,my_gpu)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
	integer, value :: chunk_dest,chunk_size
 	integer, value :: i1,i2 
	integer, value :: arr_size
	integer, dimension(arr_size) :: Arr
	integer, value :: mx_gpu,my_gpu
	integer :: i,j,k,ind,stat,ind1
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(flv(n).ne.0) then 
			i=floor(x(n))
			j=floor(y(n))
			k=floor(z(n))

! 			i=floor((x(n)-3)/8)+3
! 			j=floor((y(n)-3)/8)+3
! 			k=floor((z(n)-3)/8)+3
#ifdef twoD 
            k=1
#endif 		
            ind1=my_gpu*mx_gpu*(k-1)+mx_gpu*(j-1)+i -1
			stat=atomicadd(Arr(ind1),1)
			ind=1+stat+(chunk_dest-1)*chunk_size
			!copy particle
			q(ind)=q(n)
			x(ind)=x(n)
			y(ind)=y(n)
			z(ind)=z(n)
			u(ind)=u(n)
			v(ind)=v(n)
			w(ind)=w(n)
			var1(ind)=var1(n)
			flv(ind)=flv(n)
			tag(ind)=tag(n)

			!delete the old particle
			q(n)=0
			flv(n)=0
			u(n)=0
			v(n)=0
			w(n)=0
			tag(n)=0
			var1(n)=0
		end if 
	end if  
end subroutine BucketSortCopyPrtlKernel


!--------------------------------------------------------------------------------------------
!
!   Data Transfer for saving output data 
!
!--------------------------------------------------------------------------------------------	

subroutine CopyGPUDataToHostForOutput
	if((modulo(t,prtl_save_period).eq.0).or.(modulo(t,restart_save_period).eq.0)) then
		call RecvFullDomainPrtlFromGPU
	end if 	
	
	if((modulo(t,fld_save_period).eq.0).or.(modulo(t,restart_save_period).eq.0).or.(modulo(t,spec_save_period).eq.0)) then 
		Ex=Ex_gpu; Ey=Ey_gpu; Ez=Ez_gpu;
		Bx=Bx_gpu; By=By_gpu; Bz=Bz_gpu;
		Jx=Jx_gpu
		Jy=Jy_gpu
		Jz=Jz_gpu
	end if 
end subroutine CopyGPUDataToHostForOutput

subroutine FinishSaveDataGPU
	if((modulo(t,prtl_save_period).eq.0).or.(modulo(t,restart_save_period).eq.0)) then 		
		call ClearPrtlArrCPU
	end if
end subroutine FinishSaveDataGPU

subroutine ClearPrtlArrCPU
	qp(1:used_prtl_arr_size)=0
	used_prtl_arr_size=0
end subroutine ClearPrtlArrCPU


subroutine RecvFullDomainPrtlFromGPU
	integer :: kc
	
	used_prtl_arr_size=0
	do kc=1,Nchunk_prtl_gpu
	    if(used_prtl_chunk(kc).ne.0) then 
			call CopyPrtlToHost(used_prtl_chunk(kc),used_prtl_arr_size,(kc-1)*chunk_size_prtl_gpu)
			used_prtl_arr_size=used_prtl_arr_size+used_prtl_chunk(kc)
	    end if 		
	end do

end subroutine RecvFullDomainPrtlFromGPU

			
end module particles_gpu 

