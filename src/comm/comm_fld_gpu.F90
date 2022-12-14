module comm_fld_gpu 
	use parameters
	use vars 
	use var_gpu 
	use mem_prtl
	use cudafor
	use communication
	use comm_fldprtl
	use fields
! #ifdef cyl
!     use cyl_comm_fld_gpu
! #endif	

contains
!----------------------------------------------------------------------------------------------------------	
! The following subroutines are used in packing and exchanging data between CPU and GPU 
!
!----------------------------------------------------------------------------------------------------------	
   subroutine EMfldExchangeGPU
	   call FldYZEdgeExchangeGPU(Ex,Ey,Ez,Ex_gpu,Ey_gpu,Ez_gpu)
	   call FldYZEdgeExchangeGPU(Bx,By,Bz,Bx_gpu,By_gpu,Bz_gpu)
	   
	   call FldZXEdgeExchangeGPU(Ex,Ey,Ez,Ex_gpu,Ey_gpu,Ez_gpu)
	   call FldZXEdgeExchangeGPU(Bx,By,Bz,Bx_gpu,By_gpu,Bz_gpu)

#ifndef twoD  	   
	   call FldXYEdgeExchangeGPU(Ex,Ey,Ez,Ex_gpu,Ey_gpu,Ez_gpu)
	   call FldXYEdgeExchangeGPU(Bx,By,Bz,Bx_gpu,By_gpu,Bz_gpu)
#endif 	   
	   
   end subroutine EMfldExchangeGPU 
   
   subroutine UpdateCurrentsAllEdgesGPU
!
! 	   Jx=Jx_gpu
! 	   Jy=Jy_gpu
! 	   Jz=Jz_gpu
!
!        call ExchangeYZEdgeCurrent
!        call AddImportedCurrentYZ
!        call ExchangeZXEdgeCurrent
!        call AddImportedCurrentZX
	      
	   call CurrentYZEdgeUpdateGPU(Jx,Jy,Jz,Jx_gpu,Jy_gpu,Jz_gpu)
	   call CurrentZXEdgeUpdateGPU(Jx,Jy,Jz,Jx_gpu,Jy_gpu,Jz_gpu)
#ifndef twoD
	   call CurrentXYEdgeUpdateGPU(Jx,Jy,Jz,Jx_gpu,Jy_gpu,Jz_gpu)
#endif
! Jx_gpu=Jx
! Jy_gpu=Jy
! Jz_gpu=Jz
	   
   end subroutine UpdateCurrentsAllEdgesGPU

   
   
   
   subroutine FldYZEdgeExchangeGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
   			real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
			real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	
			call CopyYZEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyYZEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mx-4,mx-3)
			call ExchangeYZEdgeField(Fldx,Fldy,Fldz)		
            call CopyYZEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,2)
			call CopyYZEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mx-2,mx)									 		 
   end subroutine FldYZEdgeExchangeGPU
   
   subroutine FldZXEdgeExchangeGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
	        real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
	        real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
			
			call CopyZXEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyZXEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,my-4,my-3)
			call ExchangeZXEdgeField(Fldx,Fldy,Fldz)
            call CopyZXEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,2)
			call CopyZXEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,my-2,my)		 	
   end subroutine FldZXEdgeExchangeGPU
   
   subroutine FldXYEdgeExchangeGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
	        real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
	        real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
			
			call CopyXYEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyXYEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mz-4,mz-3)
			call ExchangeXYEdgeField(Fldx,Fldy,Fldz)
            call CopyXYEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,2)
			call CopyXYEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mz-2,mz)		 	
   end subroutine FldXYEdgeExchangeGPU
   
   
   subroutine CurrentYZEdgeUpdateGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
   			real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
			real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	
			call CopyYZEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,5)
			call CopyYZEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mx-4,mx)
			call ExchangeYZEdgeCurrent
			call AddImportedCurrentYZ
            call CopyYZEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyYZEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mx-4,mx-2)									 		 
   end subroutine CurrentYZEdgeUpdateGPU
   
   
   subroutine CurrentZXEdgeUpdateGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
   			real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
			real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	
			call CopyZXEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,5)
			call CopyZXEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,my-4,my)
			call ExchangeZXEdgeCurrent
			call AddImportedCurrentZX
            call CopyZXEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyZXEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,my-4,my-2)									 		 
   end subroutine CurrentZXEdgeUpdateGPU
   
   subroutine CurrentXYEdgeUpdateGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu)
   			real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
			real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	
			call CopyXYEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,1,5)
			call CopyXYEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mz-4,mz)
			call ExchangeXYEdgeCurrent
			call AddImportedCurrentXY
            call CopyXYEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,3,5)
			call CopyXYEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,mz-4,mz-2)									 		 
   end subroutine CurrentXYEdgeUpdateGPU


   subroutine CopyYZEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx(i1:i2,1:my,1:mz)=Fldx_gpu(i1:i2,1:my,1:mz)
	   Fldy(i1:i2,1:my,1:mz)=Fldy_gpu(i1:i2,1:my,1:mz)
	   Fldz(i1:i2,1:my,1:mz)=Fldz_gpu(i1:i2,1:my,1:mz)
   end subroutine CopyYZEdgeToHost
   
   subroutine CopyYZEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx_gpu(i1:i2,1:my,1:mz)=Fldx(i1:i2,1:my,1:mz)
	   Fldy_gpu(i1:i2,1:my,1:mz)=Fldy(i1:i2,1:my,1:mz)
	   Fldz_gpu(i1:i2,1:my,1:mz)=Fldz(i1:i2,1:my,1:mz)
   end subroutine CopyYZEdgeToGPU
	   
   subroutine CopyZXEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx(1:mx,i1:i2,1:mz)=Fldx_gpu(1:mx,i1:i2,1:mz)
	   Fldy(1:mx,i1:i2,1:mz)=Fldy_gpu(1:mx,i1:i2,1:mz)
	   Fldz(1:mx,i1:i2,1:mz)=Fldz_gpu(1:mx,i1:i2,1:mz)
   end subroutine CopyZXEdgeToHost
   
   subroutine CopyZXEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx_gpu(1:mx,i1:i2,1:mz)=Fldx(1:mx,i1:i2,1:mz)
	   Fldy_gpu(1:mx,i1:i2,1:mz)=Fldy(1:mx,i1:i2,1:mz)
	   Fldz_gpu(1:mx,i1:i2,1:mz)=Fldz(1:mx,i1:i2,1:mz)
   end subroutine CopyZXEdgeToGPU
	   
   subroutine CopyXYEdgeToHost(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx(1:mx,1:my,i1:i2)=Fldx_gpu(1:mx,1:my,i1:i2)
	   Fldy(1:mx,1:my,i1:i2)=Fldy_gpu(1:mx,1:my,i1:i2)
	   Fldz(1:mx,1:my,i1:i2)=Fldz_gpu(1:mx,1:my,i1:i2)
   end subroutine CopyXYEdgeToHost
   
   subroutine CopyXYEdgeToGPU(Fldx,Fldy,Fldz,Fldx_gpu,Fldy_gpu,Fldz_gpu,i1,i2)
	   real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
       real(psn), dimension(mx,my,mz), device :: Fldx_gpu,Fldy_gpu,Fldz_gpu
	   integer :: i1,i2
	    
	   Fldx_gpu(1:mx,1:my,i1:i2)=Fldx(1:mx,1:my,i1:i2)
	   Fldy_gpu(1:mx,1:my,i1:i2)=Fldy(1:mx,1:my,i1:i2)
	   Fldz_gpu(1:mx,1:my,i1:i2)=Fldz(1:mx,1:my,i1:i2)
   end subroutine CopyXYEdgeToGPU

!---------------------------------------------------------------   
! Sync only one layer, used in current filtering on GPU
!---------------------------------------------------------------     
subroutine SyncYZedgeGPU(Fld,Fld_gpu)
	real(psn), dimension(mx,my,mz)         :: Fld		
	real(psn), dimension(mx,my,mz), device :: Fld_gpu
	Fld(3:3, 1:my, 1:mz)=Fld_gpu(3:3, 1:my, 1:mz)
	Fld(mx-3:mx-3 ,1:my, 1:mz)=Fld_gpu(mx-3:mx-3, 1:my, 1:mz)
	call ExchangeYZEdgeCurrent1(Fld)
	Fld_gpu(2:2, 1:my, 1:mz)=Fld(2:2, 1:my, 1:mz)
	Fld_gpu(mx-2:mx-2, 1:my, 1:mz)=Fld(mx-2:mx-2, 1:my, 1:mz)
end subroutine SyncYZedgeGPU 
 
subroutine SyncZXedgeGPU(Fld,Fld_gpu)
	real(psn), dimension(mx,my,mz)         :: Fld		
	real(psn), dimension(mx,my,mz), device :: Fld_gpu
	Fld(1:mx, 3:3, 1:mz)=Fld_gpu(1:mx, 3:3, 1:mz)
	Fld(1:mx, my-3:my-3 , 1:mz)=Fld_gpu(1:mx, my-3:my-3, 1:mz)
	call ExchangeZXEdgeCurrent1(Fld)
	Fld_gpu(1:mx, 2:2, 1:mz)=Fld(1:mx, 2:2 , 1:mz)
	Fld_gpu(1:mx, my-2:my-2, 1:mz)=Fld(1:mx, my-2:my-2, 1:mz)
end subroutine SyncZXedgeGPU  

subroutine SyncXYedgeGPU(Fld,Fld_gpu)
	real(psn), dimension(mx,my,mz)         :: Fld		
	real(psn), dimension(mx,my,mz), device :: Fld_gpu
	Fld(1:mx, 1:my, 3:3)=Fld_gpu(1:mx, 1:my, 3:3)
	Fld(1:mx, 1:my, mz-3:mz-3)=Fld_gpu(1:mx, 1:my, mz-3:mz-3)
	call ExchangeXYEdgeCurrent1(Fld)
	Fld_gpu(1:mx, 1:my, 2:2)=Fld(1:mx, 1:my, 2:2)
	Fld_gpu(1:mx, 1:my, mz-2:mz-2)=Fld(1:mx, 1:my, mz-2:mz-2)
end subroutine SyncXYedgeGPU

!---------------------------------------------------------------   
! Some unilitiy routines: should eventaully be moved/removed if using newwe compilers
!---------------------------------------------------------------     

 	subroutine ResetMatrixGPU(Fld)
		 real(psn), dimension(mx,my,mz), device :: Fld	
         call ResetMatrixGPUKernel<<<grid,tBlock>>>(Fld,mx,my,mz)			
 	end subroutine ResetMatrixGPU	
 	attributes(global) subroutine ResetMatrixGPUKernel(Fld,mx,my,mz)
 		integer, value :: mx,my,mz
 		real(psn), dimension(mx,my,mz) :: Fld 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#else
        k=1
#endif  
        if((i.le.mx).and.(j.le.my).and.(k.le.mz)) Fld(i,j,k)=0.0
 	end subroutine ResetMatrixGPUKernel
	
	
	
 	subroutine CopyMatrixGPU(FromFld,ToFld)	
		 real(psn), dimension(mx,my,mz), device :: FromFld,ToFld
          call CopyMatrixGPUKernel<<<grid,tBlock>>>(FromFld,ToFld,mx,my,mz)			
 	end subroutine CopyMatrixGPU	
 	attributes(global) subroutine CopyMatrixGPUKernel(FromFld,ToFld,mx,my,mz)
 		integer, value :: mx,my,mz
 		real, dimension(mx,my,mz) :: FromFld,ToFld
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#else
        k=1
#endif  
        if((i.le.mx).and.(j.le.my).and.(k.le.mz)) ToFld(i,j,k)=FromFld(i,j,k)
 	end subroutine CopyMatrixGPUKernel
   
	
end module comm_fld_gpu	