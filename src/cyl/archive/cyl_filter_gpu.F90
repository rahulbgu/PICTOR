module cyl_filter_gpu
    use parameters
    use vars
	use cyl_vars
	use comm_fld_gpu
	implicit none 
contains 
    
	subroutine smoothen_current_gpu_cyl
         integer :: ni
         do ni=1,curr_filter
              call Filter_cyl_gpu(Jx,Jx_gpu,buffJx_gpu,0.5_psn)
              call Filter_cyl_gpu(Jy,Jy_gpu,buffJy_gpu,0.0_psn)
              call Filter_cyl_gpu(Jz,Jz_gpu,buffJz_gpu,0.0_psn)
         end do
    end subroutine smoothen_current_gpu_cyl
	
    subroutine Filter_cyl_gpu(Fld,Fld_gpu,FldTemp_gpu,shift)
   		real(psn), dimension(mx,my,mz)         :: Fld		
   		real(psn), dimension(mx,my,mz), device :: Fld_gpu,FldTemp_gpu   
        real(psn) :: shift
	
  		call ResetMatrixGPU(FldTemp_gpu)
		  
  		call SyncZXedgeGPU(Fld,Fld_gpu)  
  		call AvgY_GPUKernel_cyl<<<grid,tBlock>>>(Fld_gpu,FldTemp_gpu,mx,my,mz,procxind,rborders(procxind),shift)		                     
  		call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)		  
		    
     end subroutine Filter_cyl_gpu
	 
	 attributes(global) subroutine AvgY_GPUKernel_cyl(Fld1,Fld2,mx,my,mz,procx,rmin,shift)
  		integer, value :: mx,my,mz, procx
		real, value  :: rmin, shift
		real         :: wt0,wt1,r
	    real, dimension(mx,my,mz) :: Fld1,Fld2  		
	    integer :: i,j,k,i1
		
		i1=3
		if(procx.eq.0) i1=4
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +2
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +2
#ifndef twoD 		
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +2
       if((i.ge.i1).and.(i.le.mx-3).and.(j.le.my-3).and.(k.le.mz-3)) then    
#else
       k=1
       if((i.ge.i1).and.(i.le.mx-3).and.(j.le.my-3)) then 
#endif
			r=(i-3.0)+rmin+shift
			wt0=1.0 - 0.5/(1.0+r)
			wt1= 0.25/(1.0+r) 
			Fld2(i,j,k)=wt1*Fld1(i,j-1,k) +wt0*Fld1(i,j,k) + wt1*Fld1(i,j+1,k)	 
		end if
	 end subroutine AvgY_GPUKernel_cyl
	 
end module cyl_filter_gpu