module loadbalance_gpu
    use parameters
    use vars
	use var_gpu
	use particles_gpu
	use cudafor
contains
	subroutine CalcNumPrtlX_GPU(npx,shift)
		integer :: shift
		real(dbpsn), dimension(-1:nx) :: npx
		integer, dimension(:), allocatable         :: npx_int
		integer, dimension(:), allocatable, device :: npx_gpu
		integer :: kc, indi, indf
		
		allocate(npx_gpu(-1:nx),npx_int(-1:nx))
		npx_int=0
		npx_gpu(-1:nx)=npx_int(-1:nx)
		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call CalcNumPrtlX_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,npx_gpu,nx,shift,1,used_prtl_chunk(kc))
		end do
		npx_int(-1:nx)=npx_gpu(-1:nx) 
		npx=npx_int
		deallocate(npx_gpu,npx_int)
	end subroutine CalcNumPrtlX_GPU
	
	attributes(global) subroutine CalcNumPrtlX_Kernel(q,x,npx,nx,shift,indi,indf)
    	real, dimension(:) :: q,x
		integer, value :: shift,indi,indf, nx
		integer, dimension(-1:nx) :: npx
		integer :: n, ind , stat
		
		n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (indi-1)
		if(n.gt.indf) return
        if(flvp(n).eq.0) return
		ind=x(n)+shift
		stat=atomicAdd(npx(ind),1)
	end subroutine CalcNumPrtlX_Kernel
	
	
	subroutine CountPrtlSegmentsX_GPU(borders_send,np_send,lmost_send,rmost_send)
		integer :: lmost_send,rmost_send
		integer, dimension(lmost_send:rmost_send+1):: borders_send, np_send
		integer, dimension(:), allocatable, device :: borders_send_gpu, np_send_gpu
		integer :: kc, indi, indf
		
		allocate(borders_send_gpu(lmost_send:rmost_send+1),np_send_gpu(lmost_send:rmost_send+1))
		borders_send_gpu=borders_send
		np_send_gpu=np_send
		
		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call CountPrtlSegmentsX_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,lmost_send,rmost_send,borders_send_gpu,np_send_gpu,1,used_prtl_chunk(kc))
		end do 
		np_send=np_send_gpu
		deallocate(borders_send_gpu,np_send_gpu) 
	end subroutine CountPrtlSegmentsX_GPU
	
	attributes(global) subroutine CountPrtlSegmentsX_Kernel(q,x,lmost,rmost,borders_send,np_send,indi,indf)
		real, dimension(:) :: q,x
		integer, value :: lmost,rmost, indi,indf
		integer, dimension(lmost:rmost+1) :: borders_send,np_send
		integer :: n,i,stat
		
		n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		if(n.gt.indf) return
        if(flvp(n).eq.0) return
		
        i=lmost
        do i=lmost, rmost !this can be optimised, if the need be (if large number segments break away from one proc)
             if((x(n).ge.borders_send(i)).and.(x(n).lt.borders_send(i+1))) then
                  stat=atomicAdd(np_send(i),1)                 
             end if 
        end do 
	end subroutine CountPrtlSegmentsX_Kernel  
	
	subroutine MoveSendPrtlToHost(x1,x2,shift)
		integer :: used_arr_size, x1,x2,shift
		integer, device :: np_send_gpu
		integer :: kc, indi, indf
		integer :: np_send_this
		
		do kc=1,Nchunk_prtl_gpu
			np_send_this=0
			np_send_gpu=np_send_this
			
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			
		    call MoveSendPrtlGPU<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,qp_gpu,tagp_gpu,flvp_gpu,var1p_gpu,np_send_gpu,x1,x2,shift,indi,indf,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
			
			np_send_this=np_send_gpu
			
	   		if(np_send_this.gt.0) then 
				call CopyPrtlToHost(np_send_this, used_prtl_arr_size,(empty_prtl_chunk-1)*chunk_size_prtl_gpu)
				used_prtl_arr_size = used_prtl_arr_size + np_send_this
				np_gpu=np_gpu-np_send_this
	   	    end if 

		end do 
		
	end subroutine MoveSendPrtlToHost
	
	attributes(global) subroutine MoveSendPrtlGPU(x,y,z,u,v,w,q,tag,flv,var,np_send,xmin,xmax,shift,indi,indf,off)
		real, dimension(:) :: q,x,y,z,u,v,w,var
		integer, dimension(:) :: tag,flv
		integer, value :: off
		integer, value :: xmin, xmax, shift
		integer :: np_send
		integer, value :: indi, indf
		integer :: n,InsertAt
		 
		n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		if(n.gt.indf) return
		if(flvp(n).eq.0) return
		
		if(x(n).ge.xmin.and.x(n).lt.xmax) then
             x(n)=x(n)+shift
		else 
			 InsertAt=atomicinc(np_send,1000000000)
			 !InsertAt=atomicadd(np_send_gpu,1)
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
			 
			 q(n)=0.0 !delete paricles
			 flv(n)=0.0
		end if 
		
    end subroutine MoveSendPrtlGPU 
	
	
	 
end module loadbalance_gpu