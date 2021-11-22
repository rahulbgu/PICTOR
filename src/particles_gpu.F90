module particles_gpu 
	use parameters
	use vars 
	use var_gpu 
	use memory
	use cudafor
	use communication
	use comm_fldprtl
	use comm_prtl_gpu
	!use thrust_module
	use cudadevice
#ifdef cyl
    use cyl_comm_fldprtl
#endif 	
contains 
	
	subroutine SendFullDomainPrtltoGPU
		integer :: n,kc,off
		integer, dimension(Nchunk_prtl_gpu) :: np_chunk
		!This is a direct copy of all particles
		qp_gpu=qp
		xp_gpu=xp
		yp_gpu=yp
		zp_gpu=zp
		up_gpu=up
		vp_gpu=vp
		wp_gpu=wp
		var1p_gpu=var1p
		flvp_gpu=flvp
		tagp_gpu=tagp
		!now update the used_index (on Host) for further prtl sends
		used_prtl_chunk=0
		np_chunk=0
		do kc=1,Nchunk_prtl_gpu !OpenMP candidate
			off=(kc-1)*chunk_size_prtl_gpu
			do n=1,chunk_size_prtl_gpu
				if(qp(n+off).ne.0) then 
					used_prtl_chunk(kc)=n
				    np_chunk(kc)=np_chunk(kc)+1
				end if
			end do 
		end do 	
		
		np_gpu=sum(np_chunk)	
		if(used_prtl_chunk(empty_prtl_chunk).ne.0) call Abort(12)
		qp=0 ! mark all CPU particles as deleted
		used_prtl_arr_size=0
	end subroutine SendFullDomainPrtltoGPU
	

	
	subroutine ExchangePrtlGPU
		integer :: n !used in temporary fix 
		call StartTimer(32)
		call LoadGPUOutlier
	    call LoadGPUOutlierTestPrtl
		call StopTimer(32)
		if(proc.eq.0) print*,'total time in loading the outliers', real(exec_time(32))
		
   		!GPU To CPU 
		call StartTimer(33)  
   		if(np_recv_host.gt.0) then 
   			qp_host_recv(1:np_recv_host)=qp_send_gpu(1:np_recv_host)
   			xp_host_recv(1:np_recv_host)=xp_send_gpu(1:np_recv_host)
   			yp_host_recv(1:np_recv_host)=yp_send_gpu(1:np_recv_host)
   			zp_host_recv(1:np_recv_host)=zp_send_gpu(1:np_recv_host)
   			up_host_recv(1:np_recv_host)=up_send_gpu(1:np_recv_host)
   			vp_host_recv(1:np_recv_host)=vp_send_gpu(1:np_recv_host)
   			wp_host_recv(1:np_recv_host)=wp_send_gpu(1:np_recv_host)
   			var1p_host_recv(1:np_recv_host)=var1p_send_gpu(1:np_recv_host)
   			flvp_host_recv(1:np_recv_host)=flvp_send_gpu(1:np_recv_host)
   			tagp_host_recv(1:np_recv_host)=tagp_send_gpu(1:np_recv_host)
   	    end if 
   		np_gpu=np_gpu-np_recv_host
		
		if(ntp_recv_host.gt.0) then
			qtp_host_recv(1:ntp_recv_host)=qtp_send_gpu(1:ntp_recv_host)
			xtp_host_recv(1:ntp_recv_host)=xtp_send_gpu(1:ntp_recv_host)
			ytp_host_recv(1:ntp_recv_host)=ytp_send_gpu(1:ntp_recv_host)
			ztp_host_recv(1:ntp_recv_host)=ztp_send_gpu(1:ntp_recv_host)
			utp_host_recv(1:ntp_recv_host)=utp_send_gpu(1:ntp_recv_host)
			vtp_host_recv(1:ntp_recv_host)=vtp_send_gpu(1:ntp_recv_host)
			wtp_host_recv(1:ntp_recv_host)=wtp_send_gpu(1:ntp_recv_host)
			var1tp_host_recv(1:ntp_recv_host)=var1tp_send_gpu(1:ntp_recv_host)
			flvtp_host_recv(1:ntp_recv_host)=flvtp_send_gpu(1:ntp_recv_host)
			tagtp_host_recv(1:ntp_recv_host)=tagtp_send_gpu(1:ntp_recv_host)
	    end if
		call StopTimer(33)
		!print*,'total time in GPU to CPU data transfer:', real(exec_time(33))
		ntp_gpu=ntp_gpu-ntp_recv_host
        lcross=0
        rcross=0
        tcross=0
        bcross=0
        ucross=0
        dcross=0
 		call StartTimer(34)
 		call  LoadPrtlOutliersGPU(qp_host_recv,xp_host_recv,yp_host_recv,zp_host_recv,up_host_recv,vp_host_recv,wp_host_recv,flvp_host_recv,var1p_host_recv,tagp_host_recv,gpu_prtl_arr_size,np_recv_host,&
 		      lind_host,rind_host,bind_host,tind_host,dind_host,uind_host,lc_host,rc_host,bc_host,tc_host,dc_host,uc_host,buff_size_prtl_gpu)
 	      lpcross=lcross
 	      rpcross=rcross
 	      tpcross=tcross
 	      bpcross=bcross
 	      upcross=ucross
 	      dpcross=dcross
 		  np=np-(lpcross+rpcross+tpcross+bpcross+upcross+dpcross)
   		call  LoadPrtlOutliersGPU(qtp_host_recv,xtp_host_recv,ytp_host_recv,ztp_host_recv,utp_host_recv,vtp_host_recv,wtp_host_recv,flvtp_host_recv,var1tp_host_recv,tagtp_host_recv,gpu_test_prtl_arr_size,ntp_recv_host,&
   		      lind_host,rind_host,bind_host,tind_host,dind_host,uind_host,lc_host,rc_host,bc_host,tc_host,dc_host,uc_host,buff_size_prtl_gpu)
	  	call StopTimer(34)
 	  	!print*,'total time in Loading outliers on CPU:', real(exec_time(34))
          ltpcross=lcross-lpcross
          rtpcross=rcross-rpcross
          ttpcross=tcross-tpcross
          btpcross=bcross-bpcross
          utpcross=ucross-upcross
          dtpcross=dcross-dpcross
 		  ntp=ntp-(ltpcross+rtpcross+ttpcross+btpcross+utpcross+dtpcross)
 		  call StartTimer(35)
          call SendRecvPrtlSize !test particles are transferred in the same array
#ifdef cyl
          if(inc_axis) call ExchangePrtlAxis
#endif
          call UpdateTransferInSize
          call SendRecvPrtl
		   

 		  call StopTimer(35)
 		 ! print*,'total time MPI comm of Prtl data', real(exec_time(35))	  
		  
            call StartTimer(36)
#ifdef twoD
 	   np_send_host=linp_count+rinp_count+tinp_count+binp_count
 	   ntp_send_host=lintp_count+rintp_count+tintp_count+bintp_count
#else
 	   np_send_host=linp_count+rinp_count+tinp_count+binp_count+uinp_count+dinp_count
 	   ntp_send_host=lintp_count+rintp_count+tintp_count+bintp_count+uintp_count+dintp_count
#endif
        call  UnloadPrtlOutliersGPU(qp_host_send,xp_host_send,yp_host_send,zp_host_send,up_host_send,vp_host_send,wp_host_send,flvp_host_send,var1p_host_send,tagp_host_send,gpu_prtl_arr_size,&
        1,linp_count,1,rinp_count,1,binp_count,1,tinp_count,1,dinp_count,1,uinp_count)
	   
      call  UnloadPrtlOutliersGPU(qtp_host_send,xtp_host_send,ytp_host_send,ztp_host_send,utp_host_send,vtp_host_send,wtp_host_send,flvtp_host_send,var1tp_host_send,tagtp_host_send,gpu_test_prtl_arr_size,&
      linp_count+1,linp_count+lintp_count,rinp_count+1,rinp_count+rintp_count,binp_count+1,binp_count+bintp_count,tinp_count+1,tinp_count+tintp_count,dinp_count+1,dinp_count+dintp_count,uinp_count+1,uinp_count+uintp_count)   	 
	   		 
	
		  ! print*,'total time in unloading data', real(exec_time(36))	
		
		
		call SendPrtlToGPU(buff_size_prtl_gpu,xp_host_send,yp_host_send,zp_host_send,up_host_send,vp_host_send,wp_host_send,qp_host_send,tagp_host_send,flvp_host_send,var1p_host_send,buff_size_prtl_gpu,np_send_host)
		np=np+np_send_host
		
		
        call SendTestPrtlToGPU(buff_size_test_prtl_gpu,xtp_host_send,ytp_host_send,ztp_host_send,utp_host_send,vtp_host_send,wtp_host_send,qtp_host_send,tagtp_host_send,flvtp_host_send,var1tp_host_send,buff_size_test_prtl_gpu,ntp_send_host)
		ntp=ntp+ntp_send_host
		  
		call StopTimer(36) 

	end subroutine ExchangePrtlGPU
	


	subroutine LoadGPUOutlier
		integer :: kc,indi,indf
		integer :: np_send_gpu_temp
		np_send_gpu_temp=0 
		np_send_gpu=np_send_gpu_temp

		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call LoadGPUOutlierKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(gpu_prtl_arr_size,buff_size_prtl_gpu,qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		    qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu,flvp_send_gpu,tagp_send_gpu,&
		    xmin,xmax,ymin,ymax,zmin,zmax,np_send_gpu,buff_size_prtl_gpu,indi,indf)
	    end do 
		
		np_recv_host=np_send_gpu !update the incoming prtl count on the host 	
	end subroutine LoadGPUOutlier

    subroutine SendPrtlToGPU(size,x,y,z,u,v,w,q,tag,flv,var,buff_size,ind_max)
		integer :: size
		real(psn), dimension(size) :: x,y,z,u,v,w,var,q
		integer, dimension(size) :: tag,flv
		integer :: ind_max, buff_size 
		integer :: off, ind 
		
		do off=0,ind_max-1,buff_size 
			ind=min(buff_size,ind_max-off)
  			
			qp_recv_gpu(1:ind)=q(off+1:off+ind)
  			xp_recv_gpu(1:ind)=x(off+1:off+ind)
  			yp_recv_gpu(1:ind)=y(off+1:off+ind)
  			zp_recv_gpu(1:ind)=z(off+1:off+ind)
  			up_recv_gpu(1:ind)=u(off+1:off+ind)
  			vp_recv_gpu(1:ind)=v(off+1:off+ind)
  			wp_recv_gpu(1:ind)=w(off+1:off+ind)
  			var1p_recv_gpu(1:ind)=var(off+1:off+ind)
  			flvp_recv_gpu(1:ind)=flv(off+1:off+ind)
  			tagp_recv_gpu(1:ind)=tag(off+1:off+ind)
			
			call AppendGPUPrtlArr(ind)
			np_gpu=np_gpu+ind
		end do 
	end subroutine SendPrtlToGPU
	

	subroutine AppendGPUPrtlArr(np_append)
		integer :: np_append
		integer :: kc !Kernel Call
		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted

        if(np_append.eq.0) return
		Inserted=0
		do kc=1,Nchunk_prtl_gpu
			if(kc.eq.empty_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
			indi_prtl=used_prtl_chunk(kc)+1+(kc-1)*chunk_size_prtl_gpu
			indf_prtl=chunk_size_prtl_gpu+(kc-1)*chunk_size_prtl_gpu
			indi_buff=Inserted+1
			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),np_append)
				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(gpu_prtl_arr_size,buff_size_prtl_gpu,qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
				qp_recv_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,var1p_recv_gpu,flvp_recv_gpu,tagp_recv_gpu,&
				indi_prtl,indf_prtl,indi_buff,indf_buff)

		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles
		    used_prtl_chunk(kc)=used_prtl_chunk(kc)+(indf_buff-indi_buff)+1
			if(Inserted.ge.np_append) exit ! no need to coninue
 		end do
	end subroutine AppendGPUPrtlArr
	
	attributes(global) subroutine  AppendGPUPrtlArrKernel(size1,size2,q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,i1,j1,i2,j2)
	     integer, value :: size1, size2 
		 real, dimension(size1) ::q1,x1,y1,z1,u1,v1,w1,var1
		 integer, dimension(size1) :: flv1,tag1  
	     real, dimension(size2) ::q2,x2,y2,z2,u2,v2,w2,var2
		 integer, dimension(size2) :: flv2,tag2  
		 integer, value :: i1,i2,j1,j2 
		 integer :: n1,n2
	     
		 n1 = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
		 n2 = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i2-1)
		 
		 if((n1.le.j1).and.(n2.le.j2)) then
			 q1(n1)=q2(n2)
		     x1(n1)=x2(n2)
			 y1(n1)=y2(n2)
			 z1(n1)=z2(n2)
			 u1(n1)=u2(n2)
			 v1(n1)=v2(n2)
			 w1(n1)=w2(n2)
			 var1(n1)=var2(n2)
			 flv1(n1)=flv2(n2)
			 tag1(n1)=tag2(n2)
		 end if 		 
		 
	end subroutine  AppendGPUPrtlArrKernel

#ifndef GPU_USE_INTRINSICS	
	attributes(global) subroutine LoadGPUOutlierKernel(size1,size2,q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,xmin1,xmax1,ymin1,ymax1,zmin1,zmax1,np_send_gpu,Nmax,i1,i2)
	     integer, value :: size1, size2 
		 real, dimension(size1) ::q1,x1,y1,z1,u1,v1,w1,var1
		 integer, dimension(size1) :: flv1,tag1  
	     real, dimension(size2) ::q2,x2,y2,z2,u2,v2,w2,var2
		 integer, dimension(size2) :: flv2,tag2  
		 real, value  :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1 
		 integer :: np_send_gpu
		 integer, value  :: Nmax
		 integer, value :: i1,i2
		 integer :: InsertAt
		 integer :: n
		 
		 
		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
		 
		 if(n.gt.i2) return !to avoid out of bound memory access 
		 if((q1(n).ne.0)) then	

#ifndef twoD			 
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1).or.(z1(n).lt.zmin1).or.(z1(n).gt.zmax1)) then
#else  					 
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1)) then
#endif 						 
						 InsertAt=atomicinc(np_send_gpu,400000000)
						 InsertAt=InsertAt+1
												  
						 q2(InsertAt)=q1(n)		 
					     x2(InsertAt)=x1(n)
						 y2(InsertAt)=y1(n)
						 z2(InsertAt)=z1(n)
						 u2(InsertAt)=u1(n)
						 v2(InsertAt)=v1(n)
						 w2(InsertAt)=w1(n)
						 var2(InsertAt)=var1(n)
						 flv2(InsertAt)=flv1(n)
						 tag2(InsertAt)=tag1(n)
						 
						 ! delete the particle
						 q1(n)=0
						 x1(n)=xmin1+0.7
						 y1(n)=ymin1+0.7
						 z1(n)=zmin1+0.7
						 u1(n)=0
						 v1(n)=0
						 w1(n)=0
						 var1(n)=0
						 tag1(n)=0
						 flv1(n)=0

		 end if 			  
		 end if 
		 	
	end subroutine LoadGPUOutlierKernel
#endif	
	
#ifdef GPU_USE_INTRINSICS	
	attributes(global) subroutine LoadGPUOutlierKernel(q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,xmin1,xmax1,ymin1,ymax1,zmin1,zmax1,np_send_gpu,Nmax,i1,i2)
	     real, dimension(:) ::q1,x1,y1,z1,u1,v1,w1,var1
		 integer, dimension(:) :: flv1,tag1
	     real, dimension(:) ::q2,x2,y2,z2,u2,v2,w2,var2
		 integer, dimension(:) :: flv2,tag2
		 real, value  :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1
		 integer :: np_send_gpu
		 integer, value  :: Nmax
		 integer, value :: i1,i2
		 integer :: InsertAt
		 integer :: n
		 integer :: res
		 integer :: res_pred
		 integer :: WarpID
		 integer :: warp_count
		 integer :: warp_pred

		 !integer :: ALL_TRUE
		 !DATA ALL_TRUE /B'11111111111111111111111111111111'/


		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)

		 if(n.gt.i2) return !to avoid out of bound memory access


#ifndef twoD
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1).or.(z1(n).lt.zmin1).or.(z1(n).gt.zmax1)) then
#else
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1)) then
#endif
             res=1
		 else
			 res=0
		 end if
		 if((q1(n).eq.0)) res=0

		 warp_pred=ballot(res)
		 !count=syncthreads_count(res)
		 warp_count=__popc(warp_pred)

		 if(warp_count.ne.0) then
			 WarpID=mod((threadIdx%x-1), 32) + 1 !WARP_SIZE is assumed to be 32
			 if(WarpID.eq.1) InsertAt=atomicadd(np_send_gpu,warp_count)
			 InsertAt = __shfl(InsertAt, 1)

			 res_pred=IAND(warp_pred,IBITS(B'11111111111111111111111111111111',warpID,warpID))
			 warp_count=__popc(res_pred)
			 InsertAt=InsertAt+warp_count
			         if(res.eq.1) then
						 q2(InsertAt)=q1(n)
					     x2(InsertAt)=x1(n)
						 y2(InsertAt)=y1(n)
						 z2(InsertAt)=z1(n)
						 u2(InsertAt)=u1(n)
						 v2(InsertAt)=v1(n)
						 w2(InsertAt)=w1(n)
						 var2(InsertAt)=var1(n)
						 flv2(InsertAt)=flv1(n)
						 tag2(InsertAt)=tag1(n)

						 ! delete the particle
						 q1(n)=0
						 x1(n)=xmin1+0.7
						 y1(n)=ymin1+0.7
						 z1(n)=zmin1+0.7
						 u1(n)=0
						 v1(n)=0
						 w1(n)=0
						 var1(n)=0
						 tag1(n)=0
						 flv1(n)=0
					 end if 
		    end if
	end subroutine LoadGPUOutlierKernel
#endif	
	
	
	!------------------------------------------------------------
	! Test Particle Subroutines  
	!------------------------------------------------------------
	subroutine SendFullDomainTestPrtltoGPU
		integer :: n,kc,off
		!This is simply a direct copy of all the particles from buffer to GPU prtl Arr
		qtp_gpu=qtp
		xtp_gpu=xtp
		ytp_gpu=ytp
		ztp_gpu=ztp
		utp_gpu=utp
		vtp_gpu=vtp
		wtp_gpu=wtp
		var1tp_gpu=var1tp
		flvtp_gpu=flvtp
		tagtp_gpu=tagtp
		!now update the used_index (on Host) for further prtl sends
		used_test_prtl_chunk=0
		do kc=1,Nchunk_test_prtl_gpu !OpenMP candidate
			off=(kc-1)*chunk_size_test_prtl_gpu
			do n=1,chunk_size_test_prtl_gpu
				if(qtp_host_send(n+off).ne.0) used_test_prtl_chunk(kc)=used_test_prtl_chunk(kc)+1
			end do 
		end do 	
		ntp_gpu=sum(used_test_prtl_chunk)
		qtp=0
		used_test_prtl_arr_size=0

	end subroutine SendFullDomainTestPrtltoGPU
	
    subroutine SendTestPrtlToGPU(size,x,y,z,u,v,w,q,tag,flv,var,buff_size,ind_max)
		integer :: size
		real(psn), dimension(size) :: x,y,z,u,v,w,var,q
		integer, dimension(size) :: tag,flv
		integer :: ind_max, buff_size 
		integer :: off, ind 
		
		do off=0,ind_max-1,buff_size 
			ind=min(buff_size,ind_max-off)
  			
			qtp_recv_gpu(1:ind)=q(off+1:off+ind)
  			xtp_recv_gpu(1:ind)=x(off+1:off+ind)
  			ytp_recv_gpu(1:ind)=y(off+1:off+ind)
  			ztp_recv_gpu(1:ind)=z(off+1:off+ind)
  			utp_recv_gpu(1:ind)=u(off+1:off+ind)
  			vtp_recv_gpu(1:ind)=v(off+1:off+ind)
  			wtp_recv_gpu(1:ind)=w(off+1:off+ind)
  			var1tp_recv_gpu(1:ind)=var(off+1:off+ind)
  			flvtp_recv_gpu(1:ind)=flv(off+1:off+ind)
  			tagtp_recv_gpu(1:ind)=tag(off+1:off+ind)
			
			call AppendGPUTestPrtlArr(ind)
			ntp_gpu=ntp_gpu+ind
		end do 
	end subroutine SendTestPrtlToGPU
	
	
	subroutine AppendGPUTestPrtlArr(np_append)
		integer :: np_append
		integer :: kc !Kernel Call
		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted  
		
		if(np_append.eq.0) return
		Inserted=0
		do kc=1,Nchunk_test_prtl_gpu
			if(kc.eq.empty_test_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
			indi_prtl=used_test_prtl_chunk(kc)+1+(kc-1)*chunk_size_test_prtl_gpu
			indf_prtl=chunk_size_test_prtl_gpu+(kc-1)*chunk_size_test_prtl_gpu
			indi_buff=Inserted+1
			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),np_append)	
				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(gpu_test_prtl_arr_size,buff_size_test_prtl_gpu,qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
				qtp_recv_gpu,xtp_recv_gpu,ytp_recv_gpu,ztp_recv_gpu,utp_recv_gpu,vtp_recv_gpu,wtp_recv_gpu,var1tp_recv_gpu,flvtp_recv_gpu,tagtp_recv_gpu,&
				indi_prtl,indf_prtl,indi_buff,indf_buff)
				
		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles 
		    used_test_prtl_chunk(kc)=used_test_prtl_chunk(kc)+(indf_buff-indi_buff)+1
			if(Inserted.ge.np_append) exit ! no need to coninue 
 		end do 
			
	end subroutine AppendGPUTestPrtlArr
	
	
	subroutine LoadGPUOutlierTestPrtl
		integer :: kc,indi,indf
		integer :: ntp_send_gpu_temp
		ntp_send_gpu_temp=0 
		ntp_send_gpu=ntp_send_gpu_temp

		do kc=1,Nchunk_test_prtl_gpu
		    indi=(kc-1)*chunk_size_test_prtl_gpu+1
		    indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
			!prtl kernel is used for test prtl as well
		    call LoadGPUOutlierKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(gpu_test_prtl_arr_size, buff_size_test_prtl_gpu, qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		    qtp_send_gpu,xtp_send_gpu,ytp_send_gpu,ztp_send_gpu,utp_send_gpu,vtp_send_gpu,wtp_send_gpu,var1tp_send_gpu,flvtp_send_gpu,tagtp_send_gpu,&
		    xmin,xmax,ymin,ymax,zmin,zmax,ntp_send_gpu,buff_size_test_prtl_gpu,indi,indf)		
		end do 
		ntp_recv_host=ntp_send_gpu !update the incoming prtl count on the host 	
	end subroutine LoadGPUOutlierTestPrtl
	
!--------------------------------------------------------------------------------------------
!
!   Particle Sorting 
!
!--------------------------------------------------------------------------------------------	
subroutine ReorderPrtlGPU
	 if(modulo(t,prtl_reorder_period).ne.0) return 
	 call ReorderPrtlArrGPU
	 !call RegroupPrtlArrGPU
	 !call CompactPrtlArrGPU
	 call ReorderTestPrtlArrGPU
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

	
subroutine ReorderTestPrtlArrGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: total_pcount_chunk
	
	kc=empty_test_prtl_chunk
	
	do nn=1,Nchunk_test_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_test_prtl_gpu) kc=kc-Nchunk_test_prtl_gpu
		
		   indi=(kc-1)*chunk_size_test_prtl_gpu+1
		   indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
		   !call Reset_pcount_gpu(kc)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(cell_count_gpu,Ncell_gpu)

		   call BucketCountPrtl<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,mx,my)
		   
		   !call thrustscan(cell_count_gpu,Ncell_gpu)
		   call cpu_scan_temp(Ncell_gpu)
		   total_pcount_chunk=cell_count_gpu(Ncell_gpu)
		   	
		   call BucketSortCopyPrtlKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,empty_test_prtl_chunk,chunk_size_test_prtl_gpu,mx,my)
		   
		   
		used_test_prtl_chunk(empty_test_prtl_chunk)=total_pcount_chunk
		
		used_test_prtl_chunk(kc)=0
	    empty_test_prtl_chunk=kc
		!print*,'total prtl count is',pcount_host(kc,xmax1_host+2,ymax1_host+2,k2)-1
	end do 
	

end subroutine ReorderTestPrtlArrGPU  

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
		if(q(n).ne.0) then 
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
		if(q(n).ne.0) then 
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
!   Compaction of the Prtl Arr
!
!--------------------------------------------------------------------------------------------	
subroutine CompactPrtlArrGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: total_pcount_chunk
	integer, device :: pcount_chunk_gpu
	
	integer, dimension(1:chunk_size_prtl_gpu) :: perm_mask
	integer, dimension(1:chunk_size_prtl_gpu), device :: perm_mask_gpu 
	
	call GeneratePermutationFunction(perm_mask,chunk_size_prtl_gpu)
	perm_mask_gpu=perm_mask
	
	
	kc=empty_prtl_chunk
	
	do nn=1,Nchunk_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
		
		   indi=(kc-1)*chunk_size_prtl_gpu+1
		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		   
		   total_pcount_chunk=0
		   pcount_chunk_gpu=total_pcount_chunk
		   
		   call CompactPrtlArrKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,pcount_chunk_gpu,empty_prtl_chunk,chunk_size_prtl_gpu,perm_mask_gpu)

! 		   call CompactPrtlArrKernel<<<ceiling(real(chunk_size_prtl_gpu)/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
! 		   indi,indf,pcount_chunk_gpu,empty_prtl_chunk,chunk_size_prtl_gpu,perm_mask_gpu)
		
		total_pcount_chunk= pcount_chunk_gpu   
		used_prtl_chunk(empty_prtl_chunk)=total_pcount_chunk
		
		used_prtl_chunk(kc)=0
	    empty_prtl_chunk=kc
	end do 
	
end subroutine CompactPrtlArrGPU  

attributes(global) subroutine CompactPrtlArrKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,count,chunk_dest,chunk_size,perm_mask)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
	integer, dimension(:)  :: perm_mask
	integer, value :: chunk_dest,chunk_size
 	integer, value :: i1,i2 
	integer        :: count
	integer :: i,j,k,ind,stat
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
! 	n = blockDim%x * (blockIdx%x - 1) + threadIdx%x
! 	n = perm_mask(n) + (i1-1)
	
	if(n.le.i2) then
		if(q(n).ne.0) then 
			stat=atomicadd(count,1)
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
end subroutine CompactPrtlArrKernel

!--------------------------------------------------------------------------------------------
!
!   Routines for Shuffling Array elements
!
!--------------------------------------------------------------------------------------------

subroutine GeneratePermutationFunction(perm,size)
	integer :: size
	integer, dimension(1:size) :: perm
	integer :: n, temp, ind
	real(kind=8) :: r1
	
	do n=1,size
		perm(n)=n
	end do 
	
	do n=size,2,-1
		call random_number(r1)
		ind=(n-1)*r1+1
		temp=perm(ind)
		perm(ind)=perm(n)
		perm(n)=temp
	end do 
end subroutine GeneratePermutationFunction


!--------------------------------------------------------------------------------------------
!
!   Reorder the Prtl Arr such that particles in different cells are in a group 
!
!--------------------------------------------------------------------------------------------	

subroutine RegroupPrtlArrGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: total_pcount_chunk
	integer, dimension(Ncell_gpu) :: count_temp
	integer, dimension(:), allocatable, device :: prtl_new_cell !group id for the particles 
	integer, dimension(Ncell_gpu), device :: prtl_counted ! group size can be different than ncells
	
	allocate(prtl_new_cell(chunk_size_prtl_gpu))
	
	kc=empty_prtl_chunk
	
	do nn=1,Nchunk_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
		
		   indi=(kc-1)*chunk_size_prtl_gpu+1
		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(cell_count_gpu,Ncell_gpu)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(prtl_counted,Ncell_gpu)
		   
		   call RegroupBucketCountPrtl<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,mx,my,prtl_counted,prtl_new_cell,empty_prtl_chunk,chunk_size_prtl_gpu)
		 	   
		   !call thrustscan(cell_count_gpu,Ncell_gpu)
		   call cpu_scan_temp(Ncell_gpu)
		   total_pcount_chunk=cell_count_gpu(Ncell_gpu)
		   
		   call RegroupCopyPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,empty_prtl_chunk,chunk_size_prtl_gpu,mx,my,prtl_new_cell)
		   
		used_prtl_chunk(empty_prtl_chunk)=total_pcount_chunk
		
		used_prtl_chunk(kc)=0
	    empty_prtl_chunk=kc
	end do 
	if(proc.eq.0) print*,'Total Prtl After Sorting',sum(used_prtl_chunk)
	
	deallocate(prtl_new_cell)
	
end subroutine RegroupPrtlArrGPU 

attributes(global) subroutine RegroupBucketCountPrtl(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,mx_gpu,my_gpu,prtl_counted,prtl_new_cell,chunk_dest,chunk_size)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
    integer, value :: i1,i2 
    integer, value :: arr_size 
    integer, dimension(arr_size) :: Arr 
	integer, dimension(arr_size) :: prtl_counted
	integer, dimension(:)  :: prtl_new_cell 
	integer, value :: mx_gpu,my_gpu
	integer, value :: chunk_dest,chunk_size
	
	integer :: i,j,k		
	integer :: stat
	integer :: ind
	integer :: n, new_cell
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(q(n).ne.0) then 
			i=floor(x(n))
			j=floor(y(n))
			k=floor(z(n))
#ifdef twoD
            k=1  
#endif			
            ind=my_gpu*mx_gpu*(k-1)+mx_gpu*(j-1)+i 
			stat=atomicadd(prtl_counted(ind),1)
			new_cell=mod(stat,arr_size-2)+2
			prtl_new_cell(n-i1+1)=new_cell
			
			
			stat=atomicadd(Arr(new_cell),1)
		end if 
	end if  
end subroutine RegroupBucketCountPrtl

attributes(global) subroutine RegroupCopyPrtlKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,chunk_dest,chunk_size,mx_gpu,my_gpu,prtl_new_cell)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
	integer, value :: chunk_dest,chunk_size
 	integer, value :: i1,i2 
	integer, value :: arr_size
	integer, dimension(arr_size) :: Arr
	integer, dimension(:)  :: prtl_new_cell
	integer, value :: mx_gpu,my_gpu
	integer :: i,j,k,ind,stat,ind1
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(q(n).ne.0) then 
			i=floor(x(n))
			j=floor(y(n))
			k=floor(z(n))
#ifdef twoD 
            k=1
#endif 		
            !ind1=my_gpu*mx_gpu*(k-1)+mx_gpu*(j-1)+i
			ind1=prtl_new_cell(n-i1+1) -1 
			
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
end subroutine RegroupCopyPrtlKernel



!--------------------------------------------------------------------------------------------
!
!   Data Transfer for saving output data 
!
!--------------------------------------------------------------------------------------------	

subroutine CopyGPUDataToHostForOutput
	if((modulo(t,prtl_save_period).eq.0).or.(modulo(t,restart_save_period).eq.0)) then
		call RecvFullDomainPrtlFromGPU
		!call ChunkOffsetGPU
	    !call TempAppendGPUPrtl
		!call RecvFullDomainTestPrtlFromGPU
		!call ChunkOffsetTestPrtlGPU
	    !call TempAppendGPUTestPrtl
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
	used_prtl_arr_size=0
	qp=0
	used_test_prtl_arr_size=0 
	qtp=0
end subroutine ClearPrtlArrCPU


subroutine RecvFullDomainPrtlFromGPU
	integer :: kc, off1,off2
	off1=0
	off2=0
	used_prtl_arr_size=0
	do kc=1,Nchunk_prtl_gpu
	    if(used_prtl_chunk(kc).ne.0) then 
		    qp(off1+1:off1+used_prtl_chunk(kc))=qp_gpu(off2+1:off2+used_prtl_chunk(kc))
			xp(off1+1:off1+used_prtl_chunk(kc))=xp_gpu(off2+1:off2+used_prtl_chunk(kc))
			yp(off1+1:off1+used_prtl_chunk(kc))=yp_gpu(off2+1:off2+used_prtl_chunk(kc))
			zp(off1+1:off1+used_prtl_chunk(kc))=zp_gpu(off2+1:off2+used_prtl_chunk(kc))
			up(off1+1:off1+used_prtl_chunk(kc))=up_gpu(off2+1:off2+used_prtl_chunk(kc))
			vp(off1+1:off1+used_prtl_chunk(kc))=vp_gpu(off2+1:off2+used_prtl_chunk(kc))
			wp(off1+1:off1+used_prtl_chunk(kc))=wp_gpu(off2+1:off2+used_prtl_chunk(kc))
			var1p(off1+1:off1+used_prtl_chunk(kc))=var1p_gpu(off2+1:off2+used_prtl_chunk(kc))
			flvp(off1+1:off1+used_prtl_chunk(kc))=flvp_gpu(off2+1:off2+used_prtl_chunk(kc))
			tagp(off1+1:off1+used_prtl_chunk(kc))=tagp_gpu(off2+1:off2+used_prtl_chunk(kc))
	    end if 
		off1=off1+used_prtl_chunk(kc)
		off2=off2+chunk_size_prtl_gpu
		used_prtl_arr_size=used_prtl_arr_size+used_prtl_chunk(kc)
	end do

end subroutine RecvFullDomainPrtlFromGPU


subroutine RecvFullDomainTestPrtlFromGPU
	integer :: kc, off1,off2
	off1=0
	off2=0
	used_test_prtl_arr_size=0
	do kc=1,Nchunk_test_prtl_gpu
		if(used_test_prtl_chunk(kc).eq.0) cycle
		qtp(off1+1:off1+used_test_prtl_chunk(kc))=qtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		xtp(off1+1:off1+used_test_prtl_chunk(kc))=xtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		ytp(off1+1:off1+used_test_prtl_chunk(kc))=ytp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		ztp(off1+1:off1+used_test_prtl_chunk(kc))=ztp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		utp(off1+1:off1+used_test_prtl_chunk(kc))=utp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		vtp(off1+1:off1+used_test_prtl_chunk(kc))=vtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		wtp(off1+1:off1+used_test_prtl_chunk(kc))=wtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		var1tp(off1+1:off1+used_test_prtl_chunk(kc))=var1tp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		flvtp(off1+1:off1+used_test_prtl_chunk(kc))=flvtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		tagtp(off1+1:off1+used_test_prtl_chunk(kc))=tagtp_gpu(off2+1:off2+used_test_prtl_chunk(kc))
		off1=off1+used_test_prtl_chunk(kc)
		off2=off2+chunk_size_test_prtl_gpu
		used_test_prtl_arr_size=used_test_prtl_arr_size+used_test_prtl_chunk(kc)
	end do
end subroutine RecvFullDomainTestPrtlFromGPU



! subroutine TempAppendGPUPrtl
! 	integer :: n,kc,off1
! 	do kc=1,Nchunk_prtl_gpu
! 		off1=(kc-1)*chunk_size_prtl_gpu
! 	    do n=1,used_prtl_chunk(kc)
! 			qp(n+compact_chunk_offset_cpu(kc))=qp_host_recv(n+off1)
! 			xp(n+compact_chunk_offset_cpu(kc))=xp_host_recv(n+off1)
! 			yp(n+compact_chunk_offset_cpu(kc))=yp_host_recv(n+off1)
! 			zp(n+compact_chunk_offset_cpu(kc))=zp_host_recv(n+off1)
! 			up(n+compact_chunk_offset_cpu(kc))=up_host_recv(n+off1)
! 			vp(n+compact_chunk_offset_cpu(kc))=vp_host_recv(n+off1)
! 			wp(n+compact_chunk_offset_cpu(kc))=wp_host_recv(n+off1)
! 			var1p(n+compact_chunk_offset_cpu(kc))=var1p_host_recv(n+off1)
! 			flvp(n+compact_chunk_offset_cpu(kc))=flvp_host_recv(n+off1)
! 			tagp(n+compact_chunk_offset_cpu(kc))=tagp_host_recv(n+off1)
! 	    end do
!     end do
! 	used_prtl_arr_size=used_prtl_arr_size+compact_chunk_offset_cpu(Nchunk_prtl_gpu)+used_prtl_chunk(Nchunk_prtl_gpu)
! end subroutine TempAppendGPUPrtl
! subroutine DeleteTempGPUPrtl
! 	integer :: n,kc,off1
! 	do kc=1,Nchunk_prtl_gpu
! 		off1=(kc-1)*chunk_size_prtl_gpu
! 	    do n=1,used_prtl_chunk(kc)
! 			qp(n+compact_chunk_offset_cpu(kc))=0
! 			xp(n+compact_chunk_offset_cpu(kc))=xmin+0.5
! 			yp(n+compact_chunk_offset_cpu(kc))=ymin+0.5
! 			zp(n+compact_chunk_offset_cpu(kc))=zmin+0.5
! 			up(n+compact_chunk_offset_cpu(kc))=0
! 			vp(n+compact_chunk_offset_cpu(kc))=0
! 			wp(n+compact_chunk_offset_cpu(kc))=0
! 			var1p(n+compact_chunk_offset_cpu(kc))=0
! 			flvp(n+compact_chunk_offset_cpu(kc))=0
! 			tagp(n+compact_chunk_offset_cpu(kc))=0
! 	    end do
!     end do
! 	used_prtl_arr_size=used_prtl_arr_size-compact_chunk_offset_cpu(Nchunk_prtl_gpu)-used_prtl_chunk(Nchunk_prtl_gpu)
! end subroutine DeleteTempGPUPrtl
! subroutine ChunkOffsetGPU
! 	integer :: kc,off!,count
! 	!count=0
! 	off=used_prtl_arr_size
! 	do kc=1,Nchunk_prtl_gpu
! 		compact_chunk_offset_cpu(kc)=off
! 		off=off+used_prtl_chunk(kc)
! 	    !count=count+used_prtl_chunk(kc)
! 	end do
! end subroutine ChunkOffsetGPU

! subroutine TempAppendGPUTestPrtl
! 	integer :: n,kc,off1
! 	do kc=1,Nchunk_test_prtl_gpu
! 		off1=(kc-1)*chunk_size_test_prtl_gpu
! 	    do n=1,used_test_prtl_chunk(kc)
! 			qtp(n+compact_chunk_offset_test_prtl_cpu(kc))=qtp_host_recv(n+off1)
! 			xtp(n+compact_chunk_offset_test_prtl_cpu(kc))=xtp_host_recv(n+off1)
! 			ytp(n+compact_chunk_offset_test_prtl_cpu(kc))=ytp_host_recv(n+off1)
! 			ztp(n+compact_chunk_offset_test_prtl_cpu(kc))=ztp_host_recv(n+off1)
! 			utp(n+compact_chunk_offset_test_prtl_cpu(kc))=utp_host_recv(n+off1)
! 			vtp(n+compact_chunk_offset_test_prtl_cpu(kc))=vtp_host_recv(n+off1)
! 			wtp(n+compact_chunk_offset_test_prtl_cpu(kc))=wtp_host_recv(n+off1)
! 			var1tp(n+compact_chunk_offset_test_prtl_cpu(kc))=var1tp_host_recv(n+off1)
! 			flvtp(n+compact_chunk_offset_test_prtl_cpu(kc))=flvtp_host_recv(n+off1)
! 			tagtp(n+compact_chunk_offset_test_prtl_cpu(kc))=tagtp_host_recv(n+off1)
! 	    end do
!     end do
! 	used_test_prtl_arr_size=used_test_prtl_arr_size+compact_chunk_offset_test_prtl_cpu(Nchunk_test_prtl_gpu)+used_test_prtl_chunk(Nchunk_test_prtl_gpu)
! end subroutine TempAppendGPUTestPrtl

! subroutine DeleteTempGPUTestPrtl
! 	integer :: n,kc,off1
! 	do kc=1,Nchunk_test_prtl_gpu
! 		off1=(kc-1)*chunk_size_test_prtl_gpu
! 	    do n=1,used_test_prtl_chunk(kc)
! 			qtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			xtp(n+compact_chunk_offset_test_prtl_cpu(kc))=xmin+0.5
! 			ytp(n+compact_chunk_offset_test_prtl_cpu(kc))=ymin+0.5
! 			ztp(n+compact_chunk_offset_test_prtl_cpu(kc))=zmin+0.5
! 			utp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			vtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			wtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			var1tp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			flvtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 			tagtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
! 	    end do
!     end do
! 	used_test_prtl_arr_size=used_test_prtl_arr_size-compact_chunk_offset_test_prtl_cpu(Nchunk_test_prtl_gpu)-used_test_prtl_chunk(Nchunk_test_prtl_gpu)
! end subroutine DeleteTempGPUTestPrtl

! subroutine ChunkOffsetTestPrtlGPU
! 	integer :: kc,off!,count
! 	off=used_test_prtl_arr_size
! 	do kc=1,Nchunk_test_prtl_gpu
! 		compact_chunk_offset_test_prtl_cpu(kc)=off
! 		off=off+used_test_prtl_chunk(kc)
! 	end do
! end subroutine ChunkOffsetTestPrtlGPU


			
end module particles_gpu 

