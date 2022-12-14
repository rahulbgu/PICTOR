module em_update_gpu
	use parameters
	use vars 
	use var_gpu
	use usc_CFDTD_gpu
contains 
	
     subroutine UpdateEfieldGPU
         call UpdateEfieldGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,mx,my,mz,fldc) 
	end subroutine UpdateEfieldGPU
	
	attributes(global) subroutine UpdateEfieldGPUKernel(Ex,Ey,Ez,Bx,By,Bz,mx,my,mz,fldc)	
	    integer, value :: mx,my,mz	    
		real, dimension(mx,my,mz) :: Ex,Ey,Ez,Bx,By,Bz	
		real, value :: fldc
		integer :: i,j,k		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif  

	    if((i.ge.3).and.(i.le.mx-3)) then 
		    if((j.ge.3).and.(j.le.my-3)) then 					  
#ifndef twoD
                if((k.ge.3).and.(k.le.mz-3)) then 
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))-fldc*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
			    end if 
#else
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))
                    Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
#endif
	              
	           end if 
	       end if
     end subroutine UpdateEfieldGPUKernel
	 
     subroutine UpdateBfieldHalfGPU
		 if(curved_bc) then 
		     call UpdateBfieldHalf_USC_gpu
			 return 
		 end if  
         call UpdateBfieldHalfGPUKernel<<<grid,tBlock>>>(Bx_gpu,By_gpu,Bz_gpu,Ex_gpu,Ey_gpu,Ez_gpu,mx,my,mz,fld_halfc) 
	end subroutine UpdateBfieldHalfGPU
	
	attributes(global) subroutine UpdateBfieldHalfGPUKernel(Bx,By,Bz,Ex,Ey,Ez,mx,my,mz,fld_halfc)	
	    integer, value :: mx,my,mz	
		real, dimension(mx,my,mz) :: Bx,By,Bz,Ex,Ey,Ez 		
		real, value :: fld_halfc
		integer :: i,j,k		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif

			    if(i.le.mx-1) then 
				    if(j.le.my-1) then 
#ifndef twoD						
					   if(k.le.mz-1) then 				 

                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))+fld_halfc*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
					   end if 
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))
                         By(i,j,k)=By(i,j,k)+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
#endif
                      
                   end if 
               end if

     end subroutine UpdateBfieldHalfGPUKernel
	 
	 
 	subroutine AddCurrentGPU
          call AddCurrentGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Jx_gpu,Jy_gpu,Jz_gpu,mx,my,mz) 
 	end subroutine AddCurrentGPU
	

 	attributes(global) subroutine AddCurrentGPUKernel(Ex,Ey,Ez,Jx,Jy,Jz,mx,my,mz)
 	    integer, value :: mx,my,mz	    
 		real, dimension(mx,my,mz) :: Ex,Ey,Ez,Jx,Jy,Jz		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif  
        if((i.le.mx).and.(j.le.my).and.(k.le.mz)) then 
             Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)
             Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)
             Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)
 		end if 	
 	end subroutine AddCurrentGPUKernel 
	
 
	 
end module em_update_gpu 
