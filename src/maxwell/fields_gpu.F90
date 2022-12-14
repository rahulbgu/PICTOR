module fields_gpu
	use parameters
	use vars 
	use var_gpu
	use em_update_gpu
	use comm_fldprtl
	use comm_fld_gpu
	use fields
#ifdef cyl	 
	use cyl_filter_gpu
#endif
	implicit none
contains 
			
	!---------------------------------------------------------------------
	! Transfer Fld data to and from GPU
	!---------------------------------------------------------------------
	
	subroutine SendFullDomainEMFldtoGPU		
		Ex_gpu=Ex; Ey_gpu=Ey; Ez_gpu=Ez
		Bx_gpu=Bx; By_gpu=By; Bz_gpu=Bz
	end subroutine SendFullDomainEMFldtoGPU
	
	subroutine RecvFullDomainEMFldFromGPU

		Ex=Ex_gpu; Ey=Ey_gpu; Ez=Ez_gpu
		Bx=Bx_gpu; By=By_gpu; Bz=Bz_gpu
	end subroutine RecvFullDomainEMFldFromGPU
	
	subroutine SendFullDomainCurrToGPU
		Jx_gpu=Jx; Jy_gpu=Jy; Jz_gpu=Jz
	end subroutine SendFullDomainCurrToGPU
	
	subroutine RecvFullDomainCurrFromGPU
		Jx=Jx_gpu; Jy=Jy_gpu; Jz=Jz_gpu
	end subroutine RecvFullDomainCurrFromGPU
	
	subroutine SendCurrToGPU(i1,i2,j1,j2,k1,k2)
		integer :: i1,i2,j1,j2,k1,k2
#ifdef twoD
        k1=1; k2 = 1;
#endif 
        Jx_gpu(i1:i2,j1:j2,k1:k2) = Jx(i1:i2,j1:j2,k1:k2)
		Jy_gpu(i1:i2,j1:j2,k1:k2) = Jy(i1:i2,j1:j2,k1:k2)
		Jz_gpu(i1:i2,j1:j2,k1:k2) = Jz(i1:i2,j1:j2,k1:k2)
	end subroutine SendCurrToGPU
	
	subroutine RecvCurrFromGPU(i1,i2,j1,j2,k1,k2)
		integer :: i1,i2,j1,j2,k1,k2
#ifdef twoD
        k1=1; k2 = 1;
#endif 
        Jx(i1:i2,j1:j2,k1:k2) = Jx_gpu(i1:i2,j1:j2,k1:k2)
		Jy(i1:i2,j1:j2,k1:k2) = Jy_gpu(i1:i2,j1:j2,k1:k2)
		Jz(i1:i2,j1:j2,k1:k2) = Jz_gpu(i1:i2,j1:j2,k1:k2)
	end subroutine RecvCurrFromGPU
	

	
		
	subroutine ResetCurrentGPU
		 call ResetMatrixGPUKernel<<<grid,tBlock>>>(Jx_gpu,mx,my,mz)
		 call ResetMatrixGPUKernel<<<grid,tBlock>>>(Jy_gpu,mx,my,mz)
		 call ResetMatrixGPUKernel<<<grid,tBlock>>>(Jz_gpu,mx,my,mz)				
	end subroutine ResetCurrentGPU
			
	
	!---------------------------------------------------------------------
	! Subroutines used in filtering
	!---------------------------------------------------------------------
	! Currently gpu+cpu verison in not implemented
	subroutine SetFilteredEfldGPU
		integer :: ni
		call CopyMatrixGPU(TexEx_gpu,FilteredEx_gpu)
		call CopyMatrixGPU(TexEy_gpu,FilteredEy_gpu)
		call CopyMatrixGPU(TexEz_gpu,FilteredEz_gpu)
		
		do ni=1,nMoverEMfilter
           call MovingAverageFilterGPU(FilteredEx,FilteredEx_gpu,buffJx_gpu)
           call MovingAverageFilterGPU(FilteredEy,FilteredEy_gpu,buffJy_gpu)
           call MovingAverageFilterGPU(FilteredEz,FilteredEz_gpu,buffJz_gpu)
	    end do
		
 	   call FldYZEdgeExchangeGPU(FilteredEx,FilteredEy,FilteredEz,FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu)
 	   call FldZXEdgeExchangeGPU(FilteredEx,FilteredEy,FilteredEz,FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu)
 	   
#ifndef twoD  	   
 	   call FldXYEdgeExchangeGPU(FilteredEx,FilteredEy,FilteredEz,FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu)
#endif 
 			
	end subroutine SetFilteredEfldGPU
	
	

	
    subroutine smoothen_current_gpu
        
#ifdef cyl		  
	     call smoothen_current_gpu_cyl
#else 

         if(curr_filter_gpu.eqv..true.) then 
			 call smoothen_current_gpu_OnDevice
		 else 
			 call smoothen_current_gpu_OnHost
		 end if   
		       
#endif
    end subroutine smoothen_current_gpu
	
	subroutine smoothen_current_gpu_OnDevice
		 integer :: ni
         do ni=1,curr_filter
              call MovingAverageFilterGPU(Jx,Jx_gpu,buffJx_gpu)
              call MovingAverageFilterGPU(Jy,Jy_gpu,buffJy_gpu)
              call MovingAverageFilterGPU(Jz,Jz_gpu,buffJz_gpu)
         end do
	end subroutine smoothen_current_gpu_OnDevice
	
	subroutine smoothen_current_gpu_OnHost
		if(curr_filter.eq.0) return
		Jx = Jx_gpu
		Jy = Jy_gpu
		Jz = Jz_gpu
		call smoothen_current
		Jx_gpu = Jx
		Jy_gpu = Jy
		Jz_gpu = Jz
	end subroutine smoothen_current_gpu_OnHost

	 
	 
	subroutine MovingAverageFilterGPU(Fld,Fld_gpu,FldTemp_gpu)
		real(psn), dimension(mx,my,mz)         :: Fld		
		real(psn), dimension(mx,my,mz), device :: Fld_gpu,FldTemp_gpu 				  
		
		call ResetMatrixGPU(FldTemp_gpu)
		
		call SyncYZedgeGPU(Fld,Fld_gpu)
		call AvgX_GPUKernel<<<grid,tBlock>>>(Fld_gpu,FldTemp_gpu,mx,my,mz,wtm1,wt0,wtp1)	
		call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
		  
		call SyncZXedgeGPU(Fld,Fld_gpu)
		call AvgY_GPUKernel<<<grid,tBlock>>>(Fld_gpu,FldTemp_gpu,mx,my,mz,wtm1,wt0,wtp1)
		call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
		  
#ifndef twoD 
		call SyncXYedgeGPU(Fld,Fld_gpu)
		call AvgZ_GPUKernel<<<grid,tBlock>>>(Fld_gpu,FldTemp_gpu,mx,my,mz,wtm1,wt0,wtp1)		                       
		call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
#endif          
	end subroutine MovingAverageFilterGPU
	 
	 
	 attributes(global) subroutine AvgX_GPUKernel(Fld1,Fld2,mx,my,mz,wtm1,wt0,wtp1)
  		integer, value :: mx,my,mz
		real, value  :: wtm1,wt0,wtp1
	    real(psn), dimension(mx,my,mz) :: Fld1,Fld2  		
	    integer :: i,j,k
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +2
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +2
#ifndef twoD 		
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +2
       if((i.le.mx-3).and.(j.le.my-3).and.(k.le.mz-3)) Fld2(i,j,k)=wtm1*Fld1(i-1,j,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i+1,j,k)		     
#else
       k=1
       if((i.le.mx-3).and.(j.le.my-3)) Fld2(i,j,k)=wtm1*Fld1(i-1,j,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i+1,j,k)		     
#endif
	 end subroutine AvgX_GPUKernel
	 
	 attributes(global) subroutine AvgY_GPUKernel(Fld1,Fld2,mx,my,mz,wtm1,wt0,wtp1)
  		integer, value :: mx,my,mz
		real, value  :: wtm1,wt0,wtp1
	    real(psn), dimension(mx,my,mz) :: Fld1,Fld2
	    integer :: i,j,k
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +2
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +2
#ifndef twoD 		
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +2
       if((i.le.mx-3).and.(j.le.my-3).and.(k.le.mz-3)) Fld2(i,j,k)=wtm1*Fld1(i,j-1,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j+1,k)     
#else
       k=1
       if((i.le.mx-3).and.(j.le.my-3)) Fld2(i,j,k)=wtm1*Fld1(i,j-1,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j+1,k) 
#endif
	 end subroutine AvgY_GPUKernel
	 
	 attributes(global) subroutine AvgZ_GPUKernel(Fld1,Fld2,mx,my,mz,wtm1,wt0,wtp1)
  		integer, value :: mx,my,mz
		real, value  :: wtm1,wt0,wtp1
	    real, dimension(mx,my,mz) :: Fld1,Fld2
	    integer :: i,j,k
	
	    i = (blockIdx%x-1)*blockDim%x + threadIdx%x +2
	    j = (blockIdx%y-1)*blockDim%y + threadIdx%y +2
	    k = (blockIdx%z-1)*blockDim%z + threadIdx%z +2
		
        if((i.le.mx-3).and.(j.le.my-3).and.(k.le.mz-3)) Fld2(i,j,k)=wtm1*Fld1(i,j,k-1)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j,k+1)
	 end subroutine AvgZ_GPUKernel
		 	 	 
	
		
	
end module fields_gpu

