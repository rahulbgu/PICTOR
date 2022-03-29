module injector_gpu
	use parameters
	use vars
	use var_gpu
	use particles_gpu
contains 
	
	

	subroutine ClearOldPrtlRightGPU(xinj_local)
		real(dbpsn) :: xinj_local
		integer :: kc,indi,indf
		
		call Reset_pcount_gpu
		
		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call ClearOldPrtlRightGPU_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,up_gpu,vp_gpu,wp_gpu,xinj_local,indi,indf,pcount_gpu)
	    end do 
		np_gpu = np_gpu - gather_pcount_gpu()
	end subroutine ClearOldPrtlRightGPU
	
	attributes(global) subroutine ClearOldPrtlRightGPU_Kernel(qp,xp,up,vp,wp,xinj_local,indi,indf,pcount)
		real, dimension(:) :: qp,xp,up,vp,wp
		integer, dimension(:) :: pcount
		integer, value :: pcount_size
		real(dbpsn), value :: xinj_local
		integer, value :: indi, indf
		integer :: n, ind
		
	    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
	    if(n.gt.indf) return
		
		if((qp(n).ne.0).and.(xp(n).ge.xinj_local)) then 
			qp(n)=0
			up(n)=0
			vp(n)=0
			wp(n)=0
			ind = mod(n-1,pcount_size)+1
	        ind = atomicadd(pcount(ind),1)
		end if 
	
	end subroutine ClearOldPrtlRightGPU_Kernel
	
	
end module injector_gpu