module loadbalance_gpu
    use parameters
    use vars
	use var_gpu
	use cudafor
contains
	subroutine CalcNumPrtlX_GPU(npx,shift)
		integer :: shift
		real(dbpsn), dimension(-1:nx) :: npx
		integer, dimension(:), allocatable         :: npx_int
		integer, dimension(:), allocatable, device :: npx_gpu
		integer :: kc,indi,indf
		
		allocate(npx_gpu(-1:nx),npx_int(-1:nx))
		npx_int=0
		npx_gpu(-1:nx)=npx_int(-1:nx)
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call CalcNumPrtlX_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,npx_gpu,nx,shift,indi,indf)
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
        if(q(n).eq.0) return
		ind=x(n)+shift
		stat=atomicAdd(npx(ind),1)
	end subroutine CalcNumPrtlX_Kernel
	
	
	subroutine CountPrtlSegmentsX_GPU(x,y,z,u,v,w,q,tag,flv,var1,borders_send,np_send,lmost_send,rmost_send)
		real(psn), device, dimension(:) :: x,y,z,u,v,w,q,var1
		integer, device, dimension(:)   :: flv,tag
		integer :: lmost_send,rmost_send
		integer, dimension(lmost_send:rmost_send+1):: borders_send, np_send
		integer, dimension(:), allocatable, device :: borders_send_gpu, np_send_gpu
		integer :: kc,indi,indf
		
		allocate(borders_send_gpu(lmost_send:rmost_send+1),np_send_gpu(lmost_send:rmost_send+1))
		borders_send_gpu=borders_send
		np_send_gpu=np_send
		
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call CountPrtlSegmentsX_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(q,x,lmost_send,rmost_send,borders_send_gpu,np_send_gpu,indi,indf)
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
        if(q(n).eq.0) return
		
        i=lmost
        do i=lmost, rmost !this can be optimised, if the need be (if large number segments break away from one proc)
             if((x(n).ge.borders_send(i)).and.(x(n).lt.borders_send(i+1))) then
                  stat=atomicAdd(np_send(i),1)                 
             end if 
        end do 
	end subroutine CountPrtlSegmentsX_Kernel  
	
	subroutine MoveSendPrtlToHost(x,y,z,u,v,w,q,tag,flv,var1,used_arr_size,x1,x2,shift)
		real(psn), device, dimension(:) :: x,y,z,u,v,w,q,var1
		integer, device, dimension(:)   :: flv,tag
		integer :: used_arr_size, x1,x2,shift
		integer, device :: np_send_gpu
		integer :: kc,indi,indf
		integer :: np_send_this
		
		used_arr_size=0
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			np_send_this=0
			np_send_gpu=np_send_this
			
		    call MoveSendPrtlGPU<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(x,y,z,u,v,w,q,tag,flv,var1,np_send_gpu,x1,x2,shift,indi,indf,xp_send_gpu,yp_send_gpu,&
			zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,qp_send_gpu,tagp_send_gpu,flvp_send_gpu,var1p_send_gpu)
			
			np_send_this=np_send_gpu
			
	   		if(np_send_this.gt.0) then 
	   			qp(used_arr_size+1:used_arr_size+np_send_this)=qp_send_gpu(1:np_send_this)
	   			xp(used_arr_size+1:used_arr_size+np_send_this)=xp_send_gpu(1:np_send_this)
	   			yp(used_arr_size+1:used_arr_size+np_send_this)=yp_send_gpu(1:np_send_this)
	   			zp(used_arr_size+1:used_arr_size+np_send_this)=zp_send_gpu(1:np_send_this)
	   			up(used_arr_size+1:used_arr_size+np_send_this)=up_send_gpu(1:np_send_this)
	   			vp(used_arr_size+1:used_arr_size+np_send_this)=vp_send_gpu(1:np_send_this)
	   			wp(used_arr_size+1:used_arr_size+np_send_this)=wp_send_gpu(1:np_send_this)
	   			var1p(used_arr_size+1:used_arr_size+np_send_this)=var1p_send_gpu(1:np_send_this)
	   			flvp(used_arr_size+1:used_arr_size+np_send_this)=flvp_send_gpu(1:np_send_this)
	   			tagp(used_arr_size+1:used_arr_size+np_send_this)=tagp_send_gpu(1:np_send_this)
	   	    end if 
			used_arr_size=used_arr_size+np_send_this
		    np_gpu=np_gpu-np_send_this
		end do 
		
	end subroutine MoveSendPrtlToHost
	
	attributes(global) subroutine MoveSendPrtlGPU(x,y,z,u,v,w,q,tag,flv,var,np_send,xmin,xmax,shift,indi,indf,x1,y1,z1,u1,v1,w1,q1,tag1,flv1,var1)
		real, dimension(:) :: q,x,y,z,u,v,w,var
		integer, dimension(:) :: tag,flv
        real, dimension(:) ::q1,x1,y1,z1,u1,v1,w1,var1
	    integer, dimension(:) :: flv1,tag1  
		integer, value :: xmin, xmax, shift
		integer :: np_send
		integer, value :: indi, indf
		integer :: n,InsertAt
		 
		n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		if(n.gt.indf) return
		if(q(n).eq.0) return
		
		if(x(n).ge.xmin.and.x(n).lt.xmax) then
             x(n)=x(n)+shift
		else 
			 InsertAt=atomicinc(np_send,1000000000)
			 !InsertAt=atomicadd(np_send_gpu,1)
			 InsertAt=InsertAt+1
								  
			 q1(InsertAt)=q(n)		 
		     x1(InsertAt)=x(n)
			 y1(InsertAt)=y(n)
			 z1(InsertAt)=z(n)
			 u1(InsertAt)=u(n)
			 v1(InsertAt)=v(n)
			 w1(InsertAt)=w(n)
			 var1(InsertAt)=var(n)
			 flv1(InsertAt)=flv(n)
			 tag1(InsertAt)=tag(n)
			 
			 q(n)=0.0 !delete paricles
			 flv(n)=0.0
		end if 
		
    end subroutine MoveSendPrtlGPU 
	
	
	 
end module loadbalance_gpu