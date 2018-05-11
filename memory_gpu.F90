module memory_gpu 
	use parameters
	use vars 
	use var_gpu 
	use memory
	use cudafor
	use communication
	use comm_fldprtl
	use comm_prtl_gpu
	use thrust_module
	use cudadevice
contains 


!------------------------------------------------------------
! Field Subroutines  
!------------------------------------------------------------
    ! The Following subroutines are meant to transfer entire domain between GPU and CPU  
	subroutine SendFullDomainEMFldtoGPU
		integer :: k1,k2 
#ifdef twoD
     k1=1 
     k2=1 
#else
     k1=zmin1_host-2
	 k2=zmax1_host+2
#endif 		
		Ex_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Ex(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		Ey_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Ey(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		Ez_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Ez(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		Bx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Bx(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		By_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=By(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		Bz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Bz(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
	end subroutine SendFullDomainEMFldtoGPU
	
	subroutine RecvFullDomainEMFldFromGPU
		integer :: k1,k2 
#ifdef twoD
     k1=1 
     k2=1 
#else
     k1=zmin1_host
	 k2=zmax1_host-1
#endif 
		Ex(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=Ex_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
		Ey(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=Ey_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
		Ez(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=Ez_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
		Bx(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=Bx_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
		By(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=By_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
		Bz(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)=Bz_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,k1:k2)
	end subroutine RecvFullDomainEMFldFromGPU
		 		
	subroutine RecvFullDomainCurrentFromGPU
		integer :: i,j,k
		integer :: k1,k2 
#ifdef twoD
     k1=1
     k2=1
#else
     k1=zmin1_host-2
	 k2=zmax1_host+2
#endif 		
		!First copy the current from GPU
		buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
		buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)

#ifdef twoD 
        do k=1,1
#else			 		
		do k=zmin1_host-2,zmax1_host+2
#endif			
			do j=ymin1_host-2,ymax1_host+2
				do i=xmin1_host-2,xmax1_host+2
					Jx(i,j,k)=Jx(i,j,k)+buffJx_host(i,j,k)
					Jy(i,j,k)=Jy(i,j,k)+buffJy_host(i,j,k)
					Jz(i,j,k)=Jz(i,j,k)+buffJz_host(i,j,k)
				end do 
			end do 
		end do 		
	end subroutine RecvFullDomainCurrentFromGPU
	
	subroutine ResetCurrentGPU
	     type(dim3) :: grid, tBlock
			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
         call ResetCurrentGPUKernel<<<grid,tBlock>>>(Jx_gpu,Jy_gpu,Jz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)	
		
	end subroutine ResetCurrentGPU
	
	subroutine ExchangeEMfldGPU
		call ExchangeFldEdges_GPU(Ex,Ex_gpu)
		call ExchangeFldEdges_GPU(Ey,Ey_gpu)
		call ExchangeFldEdges_GPU(Ez,Ez_gpu)

		call ExchangeFldEdges_GPU(Bx,Bx_gpu)
		call ExchangeFldEdges_GPU(By,By_gpu)
		call ExchangeFldEdges_GPU(Bz,Bz_gpu)
	end subroutine ExchangeEMfldGPU
	
	subroutine ExchangeFldEdges_GPU(Fld,Fld_gpu)
		real(psn), dimension(mx,my,mz) :: Fld
#ifdef twoD		
		real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 
        real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif		
        integer :: i1,i2,j1,j2,k1,k2,z1,z2
		i1=xmin1_host-2
		i2=xmax1_host+2
		j1=ymin1_host-2
		j2=ymax1_host+2
#ifndef twoD		
		k1=zmin1_host-2
		k2=zmax1_host+2
		z1=zmin1_host
		z2=zmax1_host-1
#else
        k1=1
        k2=1
		z1=1
		z2=1
#endif
		
			Fld(xmin1_host:xmin1_host+2,ymin1_host:ymax1_host-1,z1:z2)=Fld_gpu(xmin1_host:xmin1_host+2,ymin1_host:ymax1_host-1,z1:z2) 
			Fld_gpu(xmin1_host-2:xmin1_host-1,j1:j2,k1:k2)=Fld(xmin1_host-2:xmin1_host-1,j1:j2,k1:k2)
			Fld(xmax1_host-2:xmax1_host-1,ymin1_host:ymax1_host-1,z1:z2)=Fld_gpu(xmax1_host-2:xmax1_host-1,ymin1_host:ymax1_host-1,z1:z2)
			Fld_gpu(xmax1_host:xmax1_host+2,j1:j2,k1:k2)=Fld(xmax1_host:xmax1_host+2,j1:j2,k1:k2)
			
			Fld(xmin1_host:xmax1_host-1,ymin1_host:ymin1_host+2,z1:z2)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymin1_host+2,z1:z2) 
			Fld_gpu(i1:i2,ymin1_host-2:ymin1_host-1,k1:k2)=Fld(i1:i2,ymin1_host-2:ymin1_host-1,k1:k2)
			Fld(xmin1_host:xmax1_host-1,ymax1_host-2:ymax1_host-1,z1:z2)=Fld_gpu(xmin1_host:xmax1_host-1,ymax1_host-2:ymax1_host-1,z1:z2)
			Fld_gpu(i1:i2,ymax1_host:ymax1_host+2,k1:k2)=Fld(i1:i2,ymax1_host:ymax1_host+2,k1:k2)
#ifndef twoD 			
			Fld(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmin1_host:zmin1_host+2)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmin1_host:zmin1_host+2) 
			Fld_gpu(i1:i2,j1:j2,zmin1_host-2:zmin1_host-1)=Fld(i1:i2,j1:j2,zmin1_host-2:zmin1_host-1)
			Fld(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmax1_host-2:zmax1_host-1)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmax1_host-2:zmax1_host-1)
			Fld_gpu(i1:i2,j1:j2,zmax1_host:zmax1_host+2)=Fld(i1:i2,j1:j2,zmax1_host:zmax1_host+2)
#endif						
	end subroutine ExchangeFldEdges_GPU
	
	subroutine UpdateCurrentsAllEdgesGPU
		integer :: i,j,k
		
		call ExchangeCurrentGPU
		call AddFullDomainBuffCurrentGPU
		!Add current on host
#ifdef twoD 
        do k=1,1 
#else 			 		
		do k=zmin1_host-2,zmax1_host+2
#endif 			
			do j=ymin1_host-2,ymax1_host+2
				do i=xmin1_host-2,xmax1_host+2
					Jx(i,j,k)=Jx(i,j,k)+buffJx_host(i,j,k)
					Jy(i,j,k)=Jy(i,j,k)+buffJy_host(i,j,k)
					Jz(i,j,k)=Jz(i,j,k)+buffJz_host(i,j,k)
				end do 
			end do 
		end do 
	end subroutine UpdateCurrentsAllEdgesGPU
	
	subroutine ExchangeCurrentEdges_GPU(Fld,Fld_gpu,Fld_buff,Fld_buff_gpu)
			real(psn), dimension(mx,my,mz) :: Fld
#ifdef twoD			
			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu,Fld_buff_gpu
			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1)         :: Fld_buff
#else 			
			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu,Fld_buff_gpu
			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2)        :: Fld_buff
#endif 		
            integer :: i1,i2,j1,j2,k1,k2 
		i1=xmin1_host-2
		i2=xmax1_host+2
		j1=ymin1_host-2
		j2=ymax1_host+2
#ifndef twoD		
		k1=zmin1_host-2
		k2=zmax1_host+2
#else
        k1=1
        k2=1
#endif	

				Fld_buff(xmin1_host-1:xmin1_host-1,j1:j2,k1:k2)=Fld_gpu(xmin1_host-1:xmin1_host-1,j1:j2,k1:k2) 
				Fld_buff_gpu(xmin1_host:xmin1_host,j1:j2,k1:k2)=Fld(xmin1_host:xmin1_host,j1:j2,k1:k2)
				Fld_buff(xmax1_host:xmax1_host,j1:j2,k1:k2)=Fld_gpu(xmax1_host:xmax1_host,j1:j2,k1:k2)
				Fld_buff_gpu(xmax1_host-1:xmax1_host-1,j1:j2,k1:k2)=Fld(xmax1_host-1:xmax1_host-1,j1:j2,k1:k2)
			
				Fld_buff(i1:i2,ymin1_host-1:ymin1_host-1,k1:k2)=Fld_gpu(i1:i2,ymin1_host-1:ymin1_host-1,k1:k2) 
				Fld_buff_gpu(i1:i2,ymin1_host:ymin1_host,k1:k2)=Fld(i1:i2,ymin1_host:ymin1_host,k1:k2)
				Fld_buff(i1:i2,ymax1_host:ymax1_host,k1:k2)=Fld_gpu(i1:i2,ymax1_host:ymax1_host,k1:k2)
				Fld_buff_gpu(i1:i2,ymax1_host-1:ymax1_host-1,k1:k2)=Fld(i1:i2,ymax1_host-1:ymax1_host-1,k1:k2)
#ifndef twoD 			
				Fld_buff(i1:i2,j1:j2,zmin1_host-1:zmin1_host-1)=Fld_gpu(i1:i2,j1:j2,zmin1_host-1:zmin1_host-1) 
				Fld_buff_gpu(i1:i2,j1:j2,zmin1_host:zmin1_host)=Fld(i1:i2,j1:j2,zmin1_host:zmin1_host)
				Fld_buff(i1:i2,j1:j2,zmax1_host:zmax1_host)=Fld_gpu(i1:i2,j1:j2,zmax1_host:zmax1_host)
				Fld_buff_gpu(i1:i2,j1:j2,zmax1_host-1:zmax1_host-1)=Fld(i1:i2,j1:j2,zmax1_host-1:zmax1_host-1)		
#endif 					

!
! 				Fld_buff(xmin1_host:xmin1_host+2,j1:j2,k1:k2)=Fld_gpu(xmin1_host:xmin1_host+2,j1:j2,k1:k2)
! 				Fld_buff_gpu(xmin1_host-2:xmin1_host-1,j1:j2,k1:k2)=Fld(xmin1_host-2:xmin1_host-1,j1:j2,k1:k2)
! 				Fld_buff(xmax1_host-2:xmax1_host-1,j1:j2,k1:k2)=Fld_gpu(xmax1_host-2:xmax1_host-1,j1:j2,k1:k2)
! 				Fld_buff_gpu(xmax1_host:xmax1_host+2,j1:j2,k1:k2)=Fld(xmax1_host:xmax1_host+2,j1:j2,k1:k2)
!
! 				Fld_buff(i1:i2,ymin1_host:ymin1_host+2,k1:k2)=Fld_gpu(i1:i2,ymin1_host:ymin1_host+2,k1:k2)
! 				Fld_buff_gpu(i1:i2,ymin1_host-2:ymin1_host-1,k1:k2)=Fld(i1:i2,ymin1_host-2:ymin1_host-1,k1:k2)
! 				Fld_buff(i1:i2,ymax1_host-2:ymax1_host-1,k1:k2)=Fld_gpu(i1:i2,ymax1_host-2:ymax1_host-1,k1:k2)
! 				Fld_buff_gpu(i1:i2,ymax1_host:ymax1_host+2,k1:k2)=Fld(i1:i2,ymax1_host:ymax1_host+2,k1:k2)
! #ifndef twoD
! 				Fld_buff(i1:i2,j1:j2,zmin1_host:zmin1_host+2)=Fld_gpu(i1:i2,j1:j2,zmin1_host:zmin1_host+2)
! 				Fld_buff_gpu(i1:i2,j1:j2,zmin1_host-2:zmin1_host-1)=Fld(i1:i2,j1:j2,zmin1_host-2:zmin1_host-1)
! 				Fld_buff(i1:i2,j1:j2,zmax1_host-2:zmax1_host-1)=Fld_gpu(i1:i2,j1:j2,zmax1_host-2:zmax1_host-1)
! 				Fld_buff_gpu(i1:i2,j1:j2,zmax1_host:zmax1_host+2)=Fld(i1:i2,j1:j2,zmax1_host:zmax1_host+2)
! #endif
	end subroutine ExchangeCurrentEdges_GPU
	
	subroutine ExchangeCurrentGPU
		integer :: k1,k2 
#ifdef twoD
     k1=1 
     k2=1 
#else
     k1=zmin1_host-2
	 k2=zmax1_host+2
#endif 		
		if(FullCurrentExchange) then 
			buffJx_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
			buffJy_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
			buffJz_host(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
			buffJx_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jx(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
			buffJy_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jy(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
			buffJz_gpu(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=Jz(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)			
		else 
			call ExchangeCurrentEdges_GPU(Jx,Jx_gpu,buffJx_host,buffJx_gpu)
			call ExchangeCurrentEdges_GPU(Jy,Jy_gpu,buffJy_host,buffJy_gpu)
			call ExchangeCurrentEdges_GPU(Jz,Jz_gpu,buffJz_host,buffJz_gpu)
		end if 
	end subroutine ExchangeCurrentGPU 
	subroutine AddFullDomainBuffCurrentGPU
	     type(dim3) :: grid, tBlock
			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
         call AddFullDomainBuffCurrentGPUKernel<<<grid,tBlock>>>(Jx_gpu,Jy_gpu,Jz_gpu,buffJx_gpu,buffJy_gpu,buffJz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)	
	end subroutine AddFullDomainBuffCurrentGPU

	
 
	
	!------GPU Kernels -------! 
	
	
	attributes(global) subroutine AddFullDomainBuffCurrentGPUKernel(Jx,Jy,Jz,buffJx,buffJy,buffJz,x1,x2,y1,y2,z1,z2)	
	    integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Jx,Jy,Jz
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: buffJx,buffJy,buffJz 
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: Jx,Jy,Jz
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: buffJx,buffJy,buffJz 
#endif  		
		integer :: i,j,k		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif
			    if((i.ge.x1-2).and.(i.le.x2+2)) then 
				    if((j.ge.y1-2).and.(j.le.y2+2)) then 
					    if((k.ge.z1-2).and.(k.le.z2+2)) then 	
							Jx(i,j,k)=Jx(i,j,k)+buffJx(i,j,k)
							Jy(i,j,k)=Jy(i,j,k)+buffJy(i,j,k)
							Jz(i,j,k)=Jz(i,j,k)+buffJz(i,j,k)			  
                        end if 
                   end if 
               end if

     end subroutine AddFullDomainBuffCurrentGPUKernel
	
	
	
	attributes(global) subroutine ResetCurrentGPUKernel(Jx,Jy,Jz,x1,x2,y1,y2,z1,z2)
		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
		real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Jx,Jy,Jz
#else 
        real(psn), dimension(x1-2:x2+2,y1-2:y2+2,1) :: Jx,Jy,Jz
#endif  		
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  
        if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) then 
			Jx(i,j,k)=0.0
			Jy(i,j,k)=0.0
			Jz(i,j,k)=0.0
		end if 	
	end subroutine ResetCurrentGPUKernel 
	

!------------------------------------------------------------
! Particle Subroutines  
!------------------------------------------------------------
    
	subroutine RecvFullDomainPrtlFromGPU	
		qp_host_recv=qp_gpu
		xp_host_recv=xp_gpu
		yp_host_recv=yp_gpu
		zp_host_recv=zp_gpu
		up_host_recv=up_gpu
		vp_host_recv=vp_gpu
		wp_host_recv=wp_gpu
		var1p_host_recv=var1p_gpu
		flvp_host_recv=flvp_gpu
		tagp_host_recv=tagp_gpu
		!particles are still active particles on GPU 
		!np_gpu is not updated, caller of this routine should take care of it 	
	end subroutine RecvFullDomainPrtlFromGPU
	
	subroutine SendFullDomainPrtltoGPU
		integer :: n,kc,off
		!This is simply a direct copy of all the particles from buffer to GPU prtl Arr
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
		!now update the used_index (on Host) for further prtl sends
		do kc=1,Nchunk_prtl_gpu !OpenMP candidate
			off=(kc-1)*chunk_size_prtl_gpu
			do n=1,chunk_size_prtl_gpu
				if(qp_host_send(n+off).ne.0) used_prtl_chunk(kc)=used_prtl_chunk(kc)+1
			end do 
		end do 	
		np_gpu=np_gpu+np_send_host		
	end subroutine SendFullDomainPrtltoGPU
	
	subroutine ExchangePrtlGPU
		
		call StartTimer(32)
		call LoadGPUOutlier
		call LoadHostOutLier
		call StopTimer(32)
		!print*,'total time in loading the outliers', real(exec_time(32))
		
		call StartTimer(33)
		!GPU To CPU   
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
		call StopTimer(33)
		!print*,'total time in GPU to CPU data transfer:', real(exec_time(33))
		
!implement this in the test particle segment as well 		
!  	   if((np_send_host.lt.AvailGPUSlots).and.(np_send_host.le.buff_size_prtl_gpu)) then
!  		   FullCurrentExchange=.false.
!  		   SpillOverSendHost=.true.
!  	   else
!  		   SpillOverSendHost=.false. !used in load balancing
!  		   FullCurrentExchange=.true.
!  	   end if
		!CPU to GPU 
	    call StartTimer(34)
		if(np_send_host.gt.0) then 
			qp_recv_gpu(1:np_send_host)=qp_host_send(1:np_send_host)
			xp_recv_gpu(1:np_send_host)=xp_host_send(1:np_send_host)
			yp_recv_gpu(1:np_send_host)=yp_host_send(1:np_send_host)
			zp_recv_gpu(1:np_send_host)=zp_host_send(1:np_send_host)
			up_recv_gpu(1:np_send_host)=up_host_send(1:np_send_host)
			vp_recv_gpu(1:np_send_host)=vp_host_send(1:np_send_host)
			wp_recv_gpu(1:np_send_host)=wp_host_send(1:np_send_host)
			var1p_recv_gpu(1:np_send_host)=var1p_host_send(1:np_send_host)
			flvp_recv_gpu(1:np_send_host)=flvp_host_send(1:np_send_host)
			tagp_recv_gpu(1:np_send_host)=tagp_host_send(1:np_send_host)
	    end if 
		np_gpu=np_gpu+np_send_host
		call StartTimer(34)
		!print*,'total time in CPU to GPU data transfer:', real(exec_time(34))
		
				
		
	end subroutine ExchangePrtlGPU
	
	subroutine ExchangePrtlGPU_Exclusive
		integer :: n !used in temporary fix 
		call StartTimer(32)
		call LoadGPUOutlier
	    call LoadGPUOutlierTestPrtl
		call StopTimer(32)
		!print*,'total time in loading the outliers', real(exec_time(32))
		
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
           !call UpdateTransferInSize
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
	   		 
		  call StopTimer(36)
		  ! print*,'total time in unloading data', real(exec_time(36))	
		
		
		
		  
  		if(np_send_host.gt.0) then 
  			qp_recv_gpu(1:np_send_host)=qp_host_send(1:np_send_host)
  			xp_recv_gpu(1:np_send_host)=xp_host_send(1:np_send_host)
  			yp_recv_gpu(1:np_send_host)=yp_host_send(1:np_send_host)
  			zp_recv_gpu(1:np_send_host)=zp_host_send(1:np_send_host)
  			up_recv_gpu(1:np_send_host)=up_host_send(1:np_send_host)
  			vp_recv_gpu(1:np_send_host)=vp_host_send(1:np_send_host)
  			wp_recv_gpu(1:np_send_host)=wp_host_send(1:np_send_host)
  			var1p_recv_gpu(1:np_send_host)=var1p_host_send(1:np_send_host)
  			flvp_recv_gpu(1:np_send_host)=flvp_host_send(1:np_send_host)
  			tagp_recv_gpu(1:np_send_host)=tagp_host_send(1:np_send_host)
  	    end if 
  		np_gpu=np_gpu+np_send_host
		np=np+np_send_host
		  
  		if(ntp_send_host.gt.0) then
  			qtp_recv_gpu(1:ntp_send_host)=qtp_host_send(1:ntp_send_host)
  			xtp_recv_gpu(1:ntp_send_host)=xtp_host_send(1:ntp_send_host)
  			ytp_recv_gpu(1:ntp_send_host)=ytp_host_send(1:ntp_send_host)
  			ztp_recv_gpu(1:ntp_send_host)=ztp_host_send(1:ntp_send_host)
  			utp_recv_gpu(1:ntp_send_host)=utp_host_send(1:ntp_send_host)
  			vtp_recv_gpu(1:ntp_send_host)=vtp_host_send(1:ntp_send_host)
  			wtp_recv_gpu(1:ntp_send_host)=wtp_host_send(1:ntp_send_host)
  			var1tp_recv_gpu(1:ntp_send_host)=var1tp_host_send(1:ntp_send_host)
  			flvtp_recv_gpu(1:ntp_send_host)=flvtp_host_send(1:ntp_send_host)
  			tagtp_recv_gpu(1:ntp_send_host)=tagtp_host_send(1:ntp_send_host)
  	    end if
  		ntp_gpu=ntp_gpu+np_send_host
		ntp=ntp+ntp_send_host
		   
	end subroutine ExchangePrtlGPU_Exclusive
	
	subroutine AppendHostPrtlArr
		integer :: n,off
        off=used_prtl_arr_size
		do n=1,np_recv_host
			   qp(n+off)=qp_host_recv(n)
 			   xp(n+off)=xp_host_recv(n)
			   yp(n+off)=yp_host_recv(n)
			   zp(n+off)=zp_host_recv(n)
			   up(n+off)=up_host_recv(n)
			   vp(n+off)=vp_host_recv(n)
			   wp(n+off)=wp_host_recv(n)
			   qp(n+off)=qp_host_recv(n)
			   tagp(n+off)=tagp_host_recv(n)
			   flvp(n+off)=flvp_host_recv(n)
			   var1p(n+off)=var1p_host_recv(n)
		end do
		used_prtl_arr_size=used_prtl_arr_size+np_recv_host
		prtl_random_insert_index=prtl_random_insert_index+np_recv_host		 					
			   !call InsertNewPrtl(xp_host_recv(n),yp_host_recv(n),zp_host_recv(n),up_host_recv(n),vp_host_recv(n),wp_host_recv(n),qp_host_recv(n),tagp_host_recv(n),flvp_host_recv(n),var1p_host_recv(n))	
	end subroutine AppendHostPrtlArr
![HOST] The following subroutine load the particles to be sent to GPU in to the buffer array	
	subroutine LoadHostOutlier
		integer :: i,j,count 
		count=0
		j=1  
        do i=1,used_prtl_arr_size
			
             if(qp(i).eq.0) cycle
             if((xp(i).ge.xmin_host).and.(xp(i).le.xmax_host)) then
	             if((yp(i).ge.ymin_host).and.(yp(i).le.ymax_host)) then
#ifndef twoD					 
		             if((zp(i).ge.zmin_host).and.(zp(i).le.zmax_host)) then
#endif 						 
						 qp_host_send(j)=qp(i)
						 xp_host_send(j)=xp(i)
						 yp_host_send(j)=yp(i)
						 zp_host_send(j)=zp(i)
						 up_host_send(j)=up(i)
						 vp_host_send(j)=vp(i)
						 wp_host_send(j)=wp(i)
						 var1p_host_send(j)=var1p(i)
						 flvp_host_send(j)=flvp(i)
						 tagp_host_send(j)=tagp(i)
						 
						 !Delete the particle on host  					 
					     qp(i)=0
						 xp(i)=xmin+0.5
						 yp(i)=ymin+0.5
						 zp(i)=zmin+0.5
						 up(i)=0
						 vp(i)=0
						 wp(i)=0
						 var1p(i)=0
						 flvp(i)=0
						 tagp(i)=0
						 
						 count=count+1
						 j=j+1
#ifndef twoD 						 
					 end if
#endif 					 
				 end if 
			  end if 				                               
       end do
	   np_send_host=count
	end subroutine LoadHostOutlier
	
	subroutine LoadGPUOutlier
		integer :: kc,indi,indf
		integer :: np_send_gpu_temp
		np_send_gpu_temp=0 
		np_send_gpu=np_send_gpu_temp

		do kc=1,Nchunk_prtl_gpu
		    indi=(kc-1)*chunk_size_prtl_gpu+1
		    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		    call LoadGPUOutlierKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		    qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu,flvp_send_gpu,tagp_send_gpu,&
		    xmin,xmax,ymin,ymax,zmin,zmax,np_send_gpu,buff_size_prtl_gpu,indi,indf)
	    end do 
		
		np_recv_host=np_send_gpu !update the incoming prtl count on the host 	
	end subroutine LoadGPUOutlier

!----Old subroutine -----particles are always inserted starting from same bucket 	
	subroutine AppendGPUPrtlArr
		integer :: kc !Kernel Call
		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted

		Inserted=0
		do kc=1,Nchunk_prtl_gpu
			if(kc.eq.empty_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
			indi_prtl=used_prtl_chunk(kc)+1+(kc-1)*chunk_size_prtl_gpu
			indf_prtl=chunk_size_prtl_gpu+(kc-1)*chunk_size_prtl_gpu
			indi_buff=Inserted+1
			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),np_send_host)
				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
				qp_recv_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,var1p_recv_gpu,flvp_recv_gpu,tagp_recv_gpu,&
				indi_prtl,indf_prtl,indi_buff,indf_buff)

		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles
		    used_prtl_chunk(kc)=used_prtl_chunk(kc)+(indf_buff-indi_buff)+1
			if(Inserted.ge.np_send_host) exit ! no need to coninue
 		end do

		!call UpdateAvailGPUSlots
	end subroutine AppendGPUPrtlArr
!----- New subroutines --------	
! 	subroutine AppendGPUPrtlArr
! 		integer :: nn,kc !Kernel Call
! 		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted
!
! 		Inserted=0
! 		kc=empty_prtl_chunk
! 		do nn=1,Nchunk_prtl_gpu-1
! 			kc=kc+1
! 		    if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
! 			if(kc.eq.empty_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
!
! 			indi_prtl=used_prtl_chunk(kc)+1+(kc-1)*chunk_size_prtl_gpu
! 			indf_prtl=chunk_size_prtl_gpu+(kc-1)*chunk_size_prtl_gpu
! 			indi_buff=Inserted+1
! 			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),np_send_host)
! 				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
! 				qp_recv_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,var1p_recv_gpu,flvp_recv_gpu,tagp_recv_gpu,&
! 				indi_prtl,indf_prtl,indi_buff,indf_buff)
!
! 		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles
! 		    used_prtl_chunk(kc)=used_prtl_chunk(kc)+(indf_buff-indi_buff)+1
! 			if(Inserted.ge.np_send_host) exit ! no need to coninue
!  		end do
!
! 		!call UpdateAvailGPUSlots
! 	end subroutine AppendGPUPrtlArr
	
	
	attributes(global) subroutine  AppendGPUPrtlArrKernel(q1,x1,y1,z1,u1,v1,w1,var1,flv1,tag1,q2,x2,y2,z2,u2,v2,w2,var2,flv2,tag2,i1,j1,i2,j2)
	     real, dimension(:) ::q1,x1,y1,z1,u1,v1,w1,var1
		 integer, dimension(:) :: flv1,tag1  
	     real, dimension(:) ::q2,x2,y2,z2,u2,v2,w2,var2
		 integer, dimension(:) :: flv2,tag2  
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
		 
		 
		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
		 
		 if(n.gt.i2) return !to avoid out of bound memory access 
		 if((q1(n).ne.0)) then	

#ifndef twoD			 
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1).or.(z1(n).lt.zmin1).or.(z1(n).gt.zmax1)) then
#else  					 
         if((x1(n).lt.xmin1).or.(x1(n).gt.xmax1).or.(y1(n).lt.ymin1).or.(y1(n).gt.ymax1)) then
#endif 						 
						 InsertAt=atomicinc(np_send_gpu,400000000)
						 !InsertAt=atomicadd(np_send_gpu,1)
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
						 x1(n)=xmin1+0.5
						 y1(n)=ymin1+0.5
						 z1(n)=zmin1+0.5
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
						 x1(n)=xmin1+0.5
						 y1(n)=ymin1+0.5
						 z1(n)=zmin1+0.5
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
	
	
	subroutine UpdateAvailGPUSlots
		integer :: kc 
		AvailGPUSlots=0
		do kc=1,Nchunk_prtl_gpu
			if(kc.eq.empty_prtl_chunk) cycle
			AvailGPUSlots=AvailGPUSlots+(chunk_size_prtl_gpu-used_prtl_chunk(kc))
		end do 
	end subroutine UpdateAvailGPUSlots
	!------------------------------------------------------------
	! Test Particle Subroutines  
	!------------------------------------------------------------
	subroutine RecvFullDomainTestPrtlFromGPU	
		qtp_host_recv=qtp_gpu
		xtp_host_recv=xtp_gpu
		ytp_host_recv=ytp_gpu
		ztp_host_recv=ztp_gpu
		utp_host_recv=utp_gpu
		vtp_host_recv=vtp_gpu
		wtp_host_recv=wtp_gpu
		var1tp_host_recv=var1tp_gpu
		flvtp_host_recv=flvtp_gpu
		tagtp_host_recv=tagtp_gpu
		!particles are still active particles on GPU 
		!np_gpu is not updated, caller of this routine should take care of it 	
	end subroutine RecvFullDomainTestPrtlFromGPU
	subroutine SendFullDomainTestPrtltoGPU
		integer :: n,kc,off
		!This is simply a direct copy of all the particles from buffer to GPU prtl Arr
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
		!now update the used_index (on Host) for further prtl sends
		do kc=1,Nchunk_test_prtl_gpu !OpenMP candidate
			off=(kc-1)*chunk_size_test_prtl_gpu
			do n=1,chunk_size_test_prtl_gpu
				if(qtp_host_send(n+off).ne.0) used_test_prtl_chunk(kc)=used_test_prtl_chunk(kc)+1
			end do 
		end do 	
		ntp_gpu=ntp_gpu+ntp_send_host		
	end subroutine SendFullDomainTestPrtltoGPU
	
	subroutine ExchangeTestPrtlGPU
		
		
		call LoadGPUOutlierTestPrtl
		call LoadHostOutLierTestPrtl
		
		
		!GPU To CPU   
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
		ntp_gpu=ntp_gpu-ntp_recv_host
		
		
!implement this in the test particle segment as well 		
!  	   if((np_send_host.lt.AvailGPUSlots).and.(np_send_host.le.buff_size_prtl_gpu)) then
!  		   FullCurrentExchange=.false.
!  		   SpillOverSendHost=.true.
!  	   else
!  		   SpillOverSendHost=.false. !used in load balancing
!  		   FullCurrentExchange=.true.
!  	   end if
		!CPU to GPU 
		if(ntp_send_host.gt.0) then 
			qtp_recv_gpu(1:ntp_send_host)=qtp_host_send(1:ntp_send_host)
			xtp_recv_gpu(1:ntp_send_host)=xtp_host_send(1:ntp_send_host)
			ytp_recv_gpu(1:ntp_send_host)=ytp_host_send(1:ntp_send_host)
			ztp_recv_gpu(1:ntp_send_host)=ztp_host_send(1:ntp_send_host)
			utp_recv_gpu(1:ntp_send_host)=utp_host_send(1:ntp_send_host)
			vtp_recv_gpu(1:ntp_send_host)=vtp_host_send(1:ntp_send_host)
			wtp_recv_gpu(1:ntp_send_host)=wtp_host_send(1:ntp_send_host)
			var1tp_recv_gpu(1:ntp_send_host)=var1tp_host_send(1:ntp_send_host)
			flvtp_recv_gpu(1:ntp_send_host)=flvtp_host_send(1:ntp_send_host)
			tagtp_recv_gpu(1:ntp_send_host)=tagtp_host_send(1:ntp_send_host)
	    end if 
		ntp_gpu=ntp_gpu+np_send_host
				
		
	end subroutine ExchangeTestPrtlGPU
	
	
	subroutine AppendHostTestPrtlArr
		integer :: n,off
        off=used_test_prtl_arr_size
		do n=1,ntp_recv_host
			   qtp(n+off)=qtp_host_recv(n)
 			   xtp(n+off)=xtp_host_recv(n)
			   ytp(n+off)=ytp_host_recv(n)
			   ztp(n+off)=ztp_host_recv(n)
			   utp(n+off)=utp_host_recv(n)
			   vtp(n+off)=vtp_host_recv(n)
			   wtp(n+off)=wtp_host_recv(n)
			   qtp(n+off)=qtp_host_recv(n)
			   tagtp(n+off)=tagtp_host_recv(n)
			   flvtp(n+off)=flvtp_host_recv(n)
			   var1tp(n+off)=var1tp_host_recv(n)
		end do
		used_test_prtl_arr_size=used_test_prtl_arr_size+ntp_recv_host
		test_prtl_random_insert_index=test_prtl_random_insert_index+ntp_recv_host		 						
	end subroutine AppendHostTestPrtlArr
	
	subroutine AppendGPUTestPrtlArr
		integer :: kc !Kernel Call
		integer :: indi_prtl,indi_buff,indf_prtl,indf_buff,Inserted  
		
		Inserted=0
		do kc=1,Nchunk_test_prtl_gpu
			if(kc.eq.empty_test_prtl_chunk) cycle ! This chunk is left empty for sorting, put particles elsewhere
			indi_prtl=used_test_prtl_chunk(kc)+1+(kc-1)*chunk_size_test_prtl_gpu
			indf_prtl=chunk_size_test_prtl_gpu+(kc-1)*chunk_size_test_prtl_gpu
			indi_buff=Inserted+1
			indf_buff=min(indi_buff+(indf_prtl-indi_prtl),ntp_send_host)	
				call AppendGPUPrtlArrKernel<<<ceiling(real((indf_prtl-indi_prtl+1))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
				qtp_recv_gpu,xtp_recv_gpu,ytp_recv_gpu,ztp_recv_gpu,utp_recv_gpu,vtp_recv_gpu,wtp_recv_gpu,var1tp_recv_gpu,flvtp_recv_gpu,tagtp_recv_gpu,&
				indi_prtl,indf_prtl,indi_buff,indf_buff)
				
		    Inserted=Inserted+(indf_buff-indi_buff)+1 !Now update the value of inserted particles 
		    used_test_prtl_chunk(kc)=used_test_prtl_chunk(kc)+(indf_buff-indi_buff)+1
			if(Inserted.ge.ntp_send_host) exit ! no need to coninue 
 		end do 
			
		!call UpdateAvailGPUSlotsTestPrtl
	end subroutine AppendGPUTestPrtlArr
	
	
	subroutine LoadHostOutlierTestPrtl
		integer :: i,j,count 
		count=0
		j=1  
        do i=1,used_test_prtl_arr_size
			
             if(qtp(i).eq.0) cycle
             if((xtp(i).gt.xmin_host).and.(xtp(i).lt.xmax_host)) then
	             if((ytp(i).gt.ymin_host).and.(ytp(i).lt.ymax_host)) then
#ifndef twoD					 
		             if((ztp(i).gt.zmin_host).and.(ztp(i).lt.zmax_host)) then
#endif 						 
						 qtp_host_send(j)=qtp(i)
						 xtp_host_send(j)=xtp(i)
						 ytp_host_send(j)=ytp(i)
						 ztp_host_send(j)=ztp(i)
						 utp_host_send(j)=utp(i)
						 vtp_host_send(j)=vtp(i)
						 wtp_host_send(j)=wtp(i)
						 var1tp_host_send(j)=var1tp(i)
						 flvtp_host_send(j)=flvtp(i)
						 tagtp_host_send(j)=tagtp(i)
						 
						 !Delete the particle on host  					 
					     qtp(i)=0
						 xtp(i)=xmin+0.5
						 ytp(i)=ymin+0.5
						 ztp(i)=zmin+0.5
						 utp(i)=0
						 vtp(i)=0
						 wtp(i)=0
						 var1tp(i)=0
						 flvtp(i)=0
						 tagtp(i)=0
						 
						 count=count+1
						 j=j+1
#ifndef twoD 						 
					 end if
#endif 					 
				 end if 
			  end if 				                               
       end do
	   ntp_send_host=count
	end subroutine LoadHostOutlierTestPrtl
	
	subroutine LoadGPUOutlierTestPrtl
		integer :: kc,indi,indf
		integer :: ntp_send_gpu_temp
		ntp_send_gpu_temp=0 
		ntp_send_gpu=ntp_send_gpu_temp

		do kc=1,Nchunk_test_prtl_gpu
		    indi=(kc-1)*chunk_size_test_prtl_gpu+1
		    indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
			!prtl kernel is used for test prtl as well
		    call LoadGPUOutlierKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		    qtp_send_gpu,xtp_send_gpu,ytp_send_gpu,ztp_send_gpu,utp_send_gpu,vtp_send_gpu,wtp_send_gpu,var1tp_send_gpu,flvtp_send_gpu,tagtp_send_gpu,&
		    xmin,xmax,ymin,ymax,zmin,zmax,ntp_send_gpu,buff_size_test_prtl_gpu,indi,indf)		
		end do 
		ntp_recv_host=ntp_send_gpu !update the incoming prtl count on the host 	
	end subroutine LoadGPUOutlierTestPrtl
	
	subroutine UpdateAvailGPUSlotsTestPrtl
		integer :: kc 
		AvailGPUSlotsTestPrtl=0
		do kc=1,Nchunk_test_prtl_gpu
			if(kc.eq.empty_test_prtl_chunk) cycle
			AvailGPUSlotsTestPrtl=AvailGPUSlotsTestPrtl+(chunk_size_test_prtl_gpu-used_test_prtl_chunk(kc))
		end do 
	end subroutine UpdateAvailGPUSlotsTestPrtl
	
!--------------------------------------------------------------------------------------------
!
!   Particle Sorting 
!
!--------------------------------------------------------------------------------------------	
	
subroutine ReorderPrtlGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: k1,k2
	integer :: total_pcount_chunk
	integer, dimension(Ncell_gpu) :: count_temp
	
#ifndef twoD		
		k1=zmin1_host-2
		k2=zmax1_host+2
#else
        k1=1
        k2=1
#endif

	if(modulo(t,prtl_reorder_period).ne.0) return 
	kc=empty_prtl_chunk
	!print*,'before sorting used prtl chunk',used_prtl_chunk(:) 
	
	do nn=1,Nchunk_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
		
		   indi=(kc-1)*chunk_size_prtl_gpu+1
		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(cell_count_gpu,Ncell_gpu)
		   
		   call BucketCountPrtl<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,xmin1_gpu,ymin1_gpu,zmin1_gpu,mx_gpu_host,my_gpu_host)
		 	   
		   call thrustscan(cell_count_gpu,Ncell_gpu)
		   total_pcount_chunk=cell_count_gpu(Ncell_gpu)
		   
		   call BucketSortCopyPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,empty_prtl_chunk,chunk_size_prtl_gpu,xmin1_gpu,ymin1_gpu,zmin1_gpu,mx_gpu_host,my_gpu_host)
		   
		used_prtl_chunk(empty_prtl_chunk)=total_pcount_chunk
		
		used_prtl_chunk(kc)=0
	    empty_prtl_chunk=kc
	end do 
	!print*,'After Sorting used prtl chunk',used_prtl_chunk(:)
	!print*,'Total Prtl After Sorting',sum(used_prtl_chunk)
	
	!call UpdateAvailGPUSlots !not needed at this point
end subroutine ReorderPrtlGPU  

	
subroutine ReorderTestPrtlGPU
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: k1,k2
	integer :: total_pcount_chunk
	integer, dimension(Ncell_gpu) :: count_temp
	
#ifndef twoD		
		k1=zmin1_host-2
		k2=zmax1_host+2
#else
        k1=1
        k2=1
#endif

	if(modulo(t,prtl_reorder_period).ne.0) return 
	kc=empty_test_prtl_chunk
	
	do nn=1,Nchunk_test_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_test_prtl_gpu) kc=kc-Nchunk_test_prtl_gpu
		
		   indi=(kc-1)*chunk_size_test_prtl_gpu+1
		   indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
		   !call Reset_pcount_gpu(kc)
		   call ResetCellCount<<<ceiling(real(Ncell_gpu)/NthreadsGPU), NthreadsGPU >>>(cell_count_gpu,Ncell_gpu)

		   call BucketCountPrtl<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,xmin1_gpu,ymin1_gpu,zmin1_gpu,mx_gpu_host,my_gpu_host)
		   
		   call thrustscan(cell_count_gpu,Ncell_gpu)
		   total_pcount_chunk=cell_count_gpu(Ncell_gpu)
		   	
		   call BucketSortCopyPrtlKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,flvtp_gpu,tagtp_gpu,&
		   indi,indf,cell_count_gpu,Ncell_gpu,empty_test_prtl_chunk,chunk_size_test_prtl_gpu,xmin1_gpu,ymin1_gpu,zmin1_gpu,mx_gpu_host,my_gpu_host)
		   
		   
		used_test_prtl_chunk(empty_test_prtl_chunk)=total_pcount_chunk
		
		used_test_prtl_chunk(kc)=0
	    empty_test_prtl_chunk=kc
		!print*,'total prtl count is',pcount_host(kc,xmax1_host+2,ymax1_host+2,k2)-1
	end do 
	!print*,'Total Test Prtl After Sorting',sum(used_test_prtl_chunk)
	
	
	!call UpdateAvailGPUSlots Not needed at this point
end subroutine ReorderTestPrtlGPU  


subroutine ReorderPrtlGPU_OLD
	integer :: nn,kc !chunks
	integer :: indi,indf
	integer :: k1,k2
	integer :: total_pcount_chunk
	integer, dimension(Ncell_gpu) :: count_temp
	
#ifndef twoD		
		k1=zmin1_host-2
		k2=zmax1_host+2
#else
        k1=1
        k2=1
#endif

	if(modulo(t,prtl_reorder_period).ne.0) return 
	kc=empty_prtl_chunk
	!print*,'before sorting used prtl chunk',used_prtl_chunk(:) 
	do nn=1,Nchunk_prtl_gpu-1
		   kc=kc+1
		   if(kc.gt.Nchunk_prtl_gpu) kc=kc-Nchunk_prtl_gpu
		
		   indi=(kc-1)*chunk_size_prtl_gpu+1
		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		   call Reset_pcount_gpu(kc)
		   
		   
		   call SortPrtlOffsetKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,pcount_gpu,kc,Nchunk_prtl_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)
		   
		   pcount_host(kc,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=pcount_gpu(kc,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)
 		   call SetAllPrtlBinOffset_Host(kc)
		   pcount_gpu(kc,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)=pcount_host(kc,xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,k1:k2)


		   call SortCopyPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,flvp_gpu,tagp_gpu,&
		   indi,indf,pcount_gpu,kc,empty_prtl_chunk,chunk_size_prtl_gpu,Nchunk_prtl_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)


		used_prtl_chunk(empty_prtl_chunk)=pcount_host(kc,xmax1_host+2,ymax1_host+2,k2)-1

		!used_prtl_chunk(empty_prtl_chunk)=total_pcount_chunk
		
		used_prtl_chunk(kc)=0
	    empty_prtl_chunk=kc
		!print*,'total prtl count is',pcount_host(kc,xmax1_host+2,ymax1_host+2,k2)-1
		!print*,'empty chunk is',kc,nn
	end do 
	!print*,'After Sorting used prtl chunk',used_prtl_chunk(:)
	
	call UpdateAvailGPUSlots
end subroutine ReorderPrtlGPU_OLD


subroutine SetAllPrtlBinOffset_Host(chunkID)
	integer :: chunkID
	integer :: i,j,k,offset,pcount_this
	
	offset=1
#ifndef twoD 	
	do k=zmin1_host-2,zmax1_host+2
#else 
    do k=1,1
#endif 		
		do j=ymin1_host-2,ymax1_host+2
			do i=xmin1_host-2,xmax1_host+2
                pcount_this=pcount_host(chunkID,i,j,k)
				pcount_host(chunkID,i,j,k)=offset
                offset=offset+pcount_this
			end do 
		end do 
	end do 	
end subroutine SetAllPrtlBinOffset_Host

attributes(global) subroutine SortPrtlOffsetKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,pcount,chunkID,Nc,x1,x2,y1,y2,z1,z2)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
        integer, value :: Nc 
		integer :: x1,x2,y1,y2,z1,z2 
#ifndef twoD 		
		integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,z1-2:z1+2) :: pcount
#else
        integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,1) :: pcount
#endif 		
	integer, value :: chunkID
	integer, value :: i1,i2 
	integer :: i,j,k,stat
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
			stat=atomicadd(pcount(chunkID,i,j,k),1)
		end if 
	end if  
end subroutine SortPrtlOffsetKernel

attributes(global) subroutine SortCopyPrtlKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,pcount,chunkID,chunk_dest,chunk_size,Nc,x1,x2,y1,y2,z1,z2)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
        integer, value :: Nc 
		integer :: x1,x2,y1,y2,z1,z2 
		integer, value :: chunk_dest,chunk_size
#ifndef twoD 		
		integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: pcount
#else 
        integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,1) :: pcount
#endif 	
		integer,value :: chunkID
		integer, value :: i1,i2 
		integer :: i,j,k,ind,stat
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
			stat=atomicadd(pcount(chunkID,i,j,k),1)
			ind=stat+(chunk_dest-1)*chunk_size
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
			
			!delete old particle
			q(n)=0
			flv(n)=0
			u(n)=0
			v(n)=0
			w(n)=0
			tag(n)=0
			var1(n)=0
		end if 
	end if  
end subroutine SortCopyPrtlKernel




	subroutine Reset_pcount_gpu(chunkID)
		 integer :: chunkID 
	     type(dim3) :: grid, tBlock
			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
         call Reset_pcount_gpu_Kernel<<<grid,tBlock>>>(pcount_gpu,chunkID,Nchunk_prtl_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)		
	end subroutine Reset_pcount_gpu

attributes(global) subroutine Reset_pcount_gpu_Kernel(Arr,chunkID,Nc,x1,x2,y1,y2,z1,z2) 
        integer, value :: Nc 
		integer :: x1,x2,y1,y2,z1,z2 
#ifndef twoD 		
		integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Arr
#else 
        integer, dimension(Nc,x1-2:x2+2,y1-2:y2+2,1) :: Arr
#endif 		
		integer, value :: chunkID
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) Arr(chunkID,i,j,k)=0
end subroutine Reset_pcount_gpu_Kernel

attributes(global) subroutine ResetCellCount(Arr,N)
         integer, value :: N 
		 integer, dimension(N) :: Arr
		 integer :: i
		 i = blockDim%x * (blockIdx%x - 1) + threadIdx%x
		 if(i.le.N) Arr(i)=0
end subroutine ResetCellCount  

attributes(global) subroutine BucketCountPrtl(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,x1,y1,z1,mx_gpu,my_gpu)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
    integer, value :: i1,i2 
    integer, value :: arr_size 
    integer, dimension(arr_size) :: Arr 
	integer :: x1,y1,z1
	integer, value :: mx_gpu,my_gpu
	
	integer :: i,j,k		
	integer :: stat
	integer :: ind
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(q(n).ne.0) then 
			i=floor(x(n)) -(x1-3)
			j=floor(y(n)) -(y1-3)
			k=floor(z(n)) -(z1-3)
#ifdef twoD
            k=1  
#endif			
            ind=my_gpu*mx_gpu*(k-1)+mx_gpu*(j-1)+i 
			stat=atomicadd(Arr(ind),1)
		end if 
	end if  
end subroutine BucketCountPrtl

attributes(global) subroutine BucketSortCopyPrtlKernel(q,x,y,z,u,v,w,var1,flv,tag,i1,i2,Arr,arr_size,chunk_dest,chunk_size,x1,y1,z1,mx_gpu,my_gpu)
	real, dimension(:) :: q,x,y,z,u,v,w,var1
	integer, dimension(:)  :: flv,tag
	integer, value :: chunk_dest,chunk_size
 	integer, value :: i1,i2 
	integer, value :: arr_size
	integer, dimension(arr_size) :: Arr
	integer :: x1,y1,z1
	integer, value :: mx_gpu,my_gpu
	integer :: i,j,k,ind,stat,ind1
	integer :: n
	 
    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1)
	if((n.ge.i1).and.(n.le.i2)) then
		if(q(n).ne.0) then 
			i=floor(x(n)) -(x1-3)
			j=floor(y(n)) -(y1-3)
			k=floor(z(n)) -(z1-3)
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

			!delete old particle
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
!   CPU-GPU load balancing
!
!--------------------------------------------------------------------------------------------	

subroutine AdjustDomainGPU
    if(SpillOverSendHost) then
		call CopyBackGPUSendPrtl
		xmin_host=xmin_host+0.5
		FullCurrentExchange=.true.
	end if  
	
	if(xmin_host.gt.xmin0_host) xmin_host=max(xmin_host-0.1,xmin0_host) 	
    if(xmin_host.le.xmin0_host) FullCurrentExchange=.false. 
end subroutine AdjustDomainGPU 

subroutine CopyBackGPUSendPrtl
	integer :: n,off
	off=used_prtl_arr_size
	do n=1,np_send_host
	   qp(n+off)=qp_host_send(n)
	   xp(n+off)=xp_host_send(n)
	   yp(n+off)=yp_host_send(n)
	   zp(n+off)=zp_host_send(n)
	   up(n+off)=up_host_send(n)
	   vp(n+off)=vp_host_send(n)
	   wp(n+off)=wp_host_send(n)
	   qp(n+off)=qp_host_send(n)
	   tagp(n+off)=tagp_host_send(n)
	   flvp(n+off)=flvp_host_send(n)
	   var1p(n+off)=var1p_host_send(n)
	end do
	used_prtl_arr_size=used_prtl_arr_size+np_send_host
	prtl_random_insert_index=prtl_random_insert_index+np_send_host	
	! 		if(qp_host_recv(n).ne.0) then
	! 		   call InsertNewPrtl(xp_host_send(n),yp_host_send(n),zp_host_send(n),up_host_send(n),vp_host_send(n),wp_host_send(n),qp_host_send(n),tagp_host_send(n),flvp_host_send(n),var1p_host_send(n))
	! 		end if
end subroutine CopyBackGPUSendPrtl 

!--------------------------------------------------------------------------------------------
!
!   Data Transfer for saving output data 
!
!--------------------------------------------------------------------------------------------	

subroutine CopyGPUDataToHostForOutput
	if((modulo(t,prtl_save_period).eq.0)) then
    !if((modulo(t,prtl_save_period).eq.0).or.(modulo(t,fld_save_period).eq.0)) then	 
		call RecvFullDomainPrtlFromGPU
		call ChunkOffsetGPU
	    call TempAppendGPUPrtl
		call RecvFullDomainTestPrtlFromGPU
		call ChunkOffsetTestPrtlGPU
	    call TempAppendGPUTestPrtl
	end if 	
	
	if((fld_save_period.gt.0).and.(modulo(t,fld_save_period).eq.0)) then 
#ifdef gpu	
		call RecvFullDomainEMFldFromGPU
#endif
	
#ifdef GPU_EXCLUSIVE
        Ex=Ex_gpu
		Ey=Ey_gpu
		Ez=Ez_gpu
		Bx=Bx_gpu
		By=By_gpu
		Bz=Bz_gpu
		
		Jx=Jx_gpu
		Jy=Jy_gpu
		Jz=Jz_gpu
#endif
#ifdef gpu
!Copy part of the current ! not implementd yet  
#endif
	end if 
end subroutine CopyGPUDataToHostForOutput

subroutine FinishSaveDataGPU
	if((modulo(t,prtl_save_period).eq.0)) then 
    !if((modulo(t,prtl_save_period).eq.0).or.(modulo(t,fld_save_period).eq.0)) then 	
		call DeleteTempGPUPrtl
		call DeleteTempGPUTestPrtl
	end if
end subroutine FinishSaveDataGPU

subroutine TempAppendGPUPrtl 
	integer :: n,kc,off1  
	do kc=1,Nchunk_prtl_gpu
		off1=(kc-1)*chunk_size_prtl_gpu
	    do n=1,used_prtl_chunk(kc)
			qp(n+compact_chunk_offset_cpu(kc))=qp_host_recv(n+off1)
			xp(n+compact_chunk_offset_cpu(kc))=xp_host_recv(n+off1)
			yp(n+compact_chunk_offset_cpu(kc))=yp_host_recv(n+off1)
			zp(n+compact_chunk_offset_cpu(kc))=zp_host_recv(n+off1)
			up(n+compact_chunk_offset_cpu(kc))=up_host_recv(n+off1)
			vp(n+compact_chunk_offset_cpu(kc))=vp_host_recv(n+off1)
			wp(n+compact_chunk_offset_cpu(kc))=wp_host_recv(n+off1)
			var1p(n+compact_chunk_offset_cpu(kc))=var1p_host_recv(n+off1)
			flvp(n+compact_chunk_offset_cpu(kc))=flvp_host_recv(n+off1)
			tagp(n+compact_chunk_offset_cpu(kc))=tagp_host_recv(n+off1)
	    end do 
    end do
	used_prtl_arr_size=used_prtl_arr_size+compact_chunk_offset_cpu(Nchunk_prtl_gpu)+used_prtl_chunk(Nchunk_prtl_gpu)
end subroutine TempAppendGPUPrtl 
subroutine DeleteTempGPUPrtl 
	integer :: n,kc,off1  
	do kc=1,Nchunk_prtl_gpu
		off1=(kc-1)*chunk_size_prtl_gpu
	    do n=1,used_prtl_chunk(kc)
			qp(n+compact_chunk_offset_cpu(kc))=0
			xp(n+compact_chunk_offset_cpu(kc))=xmin+0.5
			yp(n+compact_chunk_offset_cpu(kc))=ymin+0.5
			zp(n+compact_chunk_offset_cpu(kc))=zmin+0.5
			up(n+compact_chunk_offset_cpu(kc))=0
			vp(n+compact_chunk_offset_cpu(kc))=0
			wp(n+compact_chunk_offset_cpu(kc))=0
			var1p(n+compact_chunk_offset_cpu(kc))=0
			flvp(n+compact_chunk_offset_cpu(kc))=0
			tagp(n+compact_chunk_offset_cpu(kc))=0
	    end do 
    end do 
	used_prtl_arr_size=used_prtl_arr_size-compact_chunk_offset_cpu(Nchunk_prtl_gpu)-used_prtl_chunk(Nchunk_prtl_gpu)
end subroutine DeleteTempGPUPrtl
subroutine ChunkOffsetGPU
	integer :: kc,off!,count
	!count=0
	off=used_prtl_arr_size
	do kc=1,Nchunk_prtl_gpu
		compact_chunk_offset_cpu(kc)=off
		off=off+used_prtl_chunk(kc)
	    !count=count+used_prtl_chunk(kc)
	end do 
end subroutine ChunkOffsetGPU

subroutine TempAppendGPUTestPrtl 
	integer :: n,kc,off1  
	do kc=1,Nchunk_test_prtl_gpu
		off1=(kc-1)*chunk_size_test_prtl_gpu
	    do n=1,used_test_prtl_chunk(kc)
			qtp(n+compact_chunk_offset_test_prtl_cpu(kc))=qtp_host_recv(n+off1)
			xtp(n+compact_chunk_offset_test_prtl_cpu(kc))=xtp_host_recv(n+off1)
			ytp(n+compact_chunk_offset_test_prtl_cpu(kc))=ytp_host_recv(n+off1)
			ztp(n+compact_chunk_offset_test_prtl_cpu(kc))=ztp_host_recv(n+off1)
			utp(n+compact_chunk_offset_test_prtl_cpu(kc))=utp_host_recv(n+off1)
			vtp(n+compact_chunk_offset_test_prtl_cpu(kc))=vtp_host_recv(n+off1)
			wtp(n+compact_chunk_offset_test_prtl_cpu(kc))=wtp_host_recv(n+off1)
			var1tp(n+compact_chunk_offset_test_prtl_cpu(kc))=var1tp_host_recv(n+off1)
			flvtp(n+compact_chunk_offset_test_prtl_cpu(kc))=flvtp_host_recv(n+off1)
			tagtp(n+compact_chunk_offset_test_prtl_cpu(kc))=tagtp_host_recv(n+off1)
	    end do 
    end do
	used_test_prtl_arr_size=used_test_prtl_arr_size+compact_chunk_offset_test_prtl_cpu(Nchunk_test_prtl_gpu)+used_test_prtl_chunk(Nchunk_test_prtl_gpu)
end subroutine TempAppendGPUTestPrtl 

subroutine DeleteTempGPUTestPrtl 
	integer :: n,kc,off1  
	do kc=1,Nchunk_test_prtl_gpu
		off1=(kc-1)*chunk_size_test_prtl_gpu
	    do n=1,used_test_prtl_chunk(kc)
			qtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			xtp(n+compact_chunk_offset_test_prtl_cpu(kc))=xmin+0.5
			ytp(n+compact_chunk_offset_test_prtl_cpu(kc))=ymin+0.5
			ztp(n+compact_chunk_offset_test_prtl_cpu(kc))=zmin+0.5
			utp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			vtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			wtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			var1tp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			flvtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
			tagtp(n+compact_chunk_offset_test_prtl_cpu(kc))=0
	    end do 
    end do 
	used_test_prtl_arr_size=used_test_prtl_arr_size-compact_chunk_offset_test_prtl_cpu(Nchunk_test_prtl_gpu)-used_test_prtl_chunk(Nchunk_test_prtl_gpu)
end subroutine DeleteTempGPUTestPrtl

subroutine ChunkOffsetTestPrtlGPU
	integer :: kc,off!,count
	off=used_test_prtl_arr_size
	do kc=1,Nchunk_test_prtl_gpu
		compact_chunk_offset_test_prtl_cpu(kc)=off
		off=off+used_test_prtl_chunk(kc)
	end do 
end subroutine ChunkOffsetTestPrtlGPU




			
end module memory_gpu 