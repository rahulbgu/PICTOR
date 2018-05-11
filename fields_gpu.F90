module fields_gpu
	use parameters
	use vars 
	use var_gpu
	use comm_fldprtl
	use comm_fld_gpu
contains 
	
     subroutine UpdateEfieldGPU

	     type(dim3) :: grid, tBlock			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
         call UpdateEfieldGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,fldc) 
	end subroutine UpdateEfieldGPU
	
	attributes(global) subroutine UpdateEfieldGPUKernel(Ex,Ey,Ez,Bx,By,Bz,x1,x2,y1,y2,z1,z2,fldc)	
	    integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 	    
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Ex,Ey,Ez,Bx,By,Bz
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1:1) :: Ex,Ey,Ez,Bx,By,Bz
#endif 		
		real, value :: fldc
		integer :: i,j,k		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  

			    if((i.ge.x1).and.(i.lt.x2)) then 
				    if((j.ge.y1).and.(j.lt.y2)) then 
					    if((k.ge.z1).and.(k.lt.z2)) then 				  
					  
#ifndef twoD
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))-fldc*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
#else
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))
                    Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
#endif
                       end if 
                   end if 
               end if
     end subroutine UpdateEfieldGPUKernel
	 
     subroutine UpdateBfieldHalfGPU
	     type(dim3) :: grid, tBlock			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
         call UpdateBfieldHalfGPUKernel<<<grid,tBlock>>>(Bx_gpu,By_gpu,Bz_gpu,Ex_gpu,Ey_gpu,Ez_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,fld_halfc) 
	end subroutine UpdateBfieldHalfGPU
	
	attributes(global) subroutine UpdateBfieldHalfGPUKernel(Bx,By,Bz,Ex,Ey,Ez,x1,x2,y1,y2,z1,z2,fld_halfc)	
	    integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Bx,By,Bz,Ex,Ey,Ez 
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: Bx,By,Bz,Ex,Ey,Ez 
#endif 		
		real, value :: fld_halfc
		integer :: i,j,k		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif
			    if((i.ge.x1-2).and.(i.le.x2+1)) then 
				    if((j.ge.y1-2).and.(j.le.y2+1)) then 
					    if((k.ge.z1-2).and.(k.le.z2+1)) then 				  
					  
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))+fld_halfc*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))
                         By(i,j,k)=By(i,j,k)+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
#endif
                       end if 
                   end if 
               end if

     end subroutine UpdateBfieldHalfGPUKernel
	 
	 
	 
	subroutine AddCurrentGPU
	     type(dim3) :: grid, tBlock			 
#ifdef twoD			 
		 tBlock = dim3 (16 ,16 ,1)	
		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+1)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+1)/tBlock%y), 1) 
#else
         tBlock = dim3 (8 ,8 ,4)
         grid = dim3(ceiling(real(xmax1_host-xmin1_host+1)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+1)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+1)/tBlock%z)) 
#endif  
         call AddCurrentGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Jx_gpu,Jy_gpu,Jz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu) 
	end subroutine AddCurrentGPU
	

	attributes(global) subroutine AddCurrentGPUKernel(Ex,Ey,Ez,Jx,Jy,Jz,x1,x2,y1,y2,z1,z2)
	    integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 	    
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Ex,Ey,Ez,Jx,Jy,Jz
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: Ex,Ey,Ez,Jx,Jy,Jz
#endif 		
		integer :: i,j,k
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-1
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-1
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-1
#else
        k=1
#endif  
        if((i.le.x2).and.(j.le.y2).and.(k.le.z2)) then 
            Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)
            Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)
            Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)
		end if 	
	end subroutine AddCurrentGPUKernel 
	

	!---------------------------------------------------------------------
	! Subroutines used in filtering
	!---------------------------------------------------------------------
	! Currently gpu+cpu verison in not implemented
	subroutine SetFilteredEfldGPU_Exclusive
		integer :: ni
		call CopyMatrixGPU(TexEx_gpu,FilteredEx_gpu)
		call CopyMatrixGPU(TexEy_gpu,FilteredEy_gpu)
		call CopyMatrixGPU(TexEz_gpu,FilteredEz_gpu)
		
		do ni=1,nMoverEMfilter
           call MovingAverageFilterGPU_Exclusive(FilteredEx_gpu,buffJx_gpu)
           call MovingAverageFilterGPU_Exclusive(FilteredEy_gpu,buffJy_gpu)
           call MovingAverageFilterGPU_Exclusive(FilteredEz_gpu,buffJz_gpu)
	    end do
		
 	   call FldYZEdgeExchangeGPU_Exclusive(Ex,lEx_send_host,rEx_send_host,lEx_recv_host,rEx_recv_host,FilteredEx_gpu,lEx_send_gpu,rEx_send_gpu,lEx_recv_gpu,rEx_recv_gpu)
 	   call FldYZEdgeExchangeGPU_Exclusive(Ey,lEy_send_host,rEy_send_host,lEy_recv_host,rEy_recv_host,FilteredEy_gpu,lEy_send_gpu,rEy_send_gpu,lEy_recv_gpu,rEy_recv_gpu)
 	   call FldYZEdgeExchangeGPU_Exclusive(Ez,lEz_send_host,rEz_send_host,lEz_recv_host,rEz_recv_host,FilteredEz_gpu,lEz_send_gpu,rEz_send_gpu,lEz_recv_gpu,rEz_recv_gpu)
	   
 	   call FldZXEdgeExchangeGPU_Exclusive(Ex,bEx_send_host,tEx_send_host,bEx_recv_host,tEx_recv_host,FilteredEx_gpu,bEx_send_gpu,tEx_send_gpu,bEx_recv_gpu,tEx_recv_gpu)
 	   call FldZXEdgeExchangeGPU_Exclusive(Ey,bEy_send_host,tEy_send_host,bEy_recv_host,tEy_recv_host,FilteredEy_gpu,bEy_send_gpu,tEy_send_gpu,bEy_recv_gpu,tEy_recv_gpu)
 	   call FldZXEdgeExchangeGPU_Exclusive(Ez,bEz_send_host,tEz_send_host,bEz_recv_host,tEz_recv_host,FilteredEz_gpu,bEz_send_gpu,tEz_send_gpu,bEz_recv_gpu,tEz_recv_gpu)
#ifndef twoD
 	   call FldXYEdgeExchangeGPU_Exclusive(Ex,dEx_send_host,uEx_send_host,dEx_recv_host,uEx_recv_host,FilteredEx_gpu,dEx_send_gpu,uEx_send_gpu,dEx_recv_gpu,uEx_recv_gpu)
 	   call FldXYEdgeExchangeGPU_Exclusive(Ey,dEy_send_host,uEy_send_host,dEy_recv_host,uEy_recv_host,FilteredEy_gpu,dEy_send_gpu,uEy_send_gpu,dEy_recv_gpu,uEy_recv_gpu)
 	   call FldXYEdgeExchangeGPU_Exclusive(Ez,dEz_send_host,uEz_send_host,dEz_recv_host,uEz_recv_host,FilteredEz_gpu,dEz_send_gpu,uEz_send_gpu,dEz_recv_gpu,uEz_recv_gpu)
#endif 			
	end subroutine SetFilteredEfldGPU_Exclusive 	
	
	
    subroutine smoothen_current_gpu
         integer :: ni
         do ni=1,curr_filter
              call MovingAverageFilterGPU(Jx,Jx_gpu,buffJx_gpu)
              call MovingAverageFilterGPU(Jy,Jy_gpu,buffJy_gpu)
              call MovingAverageFilterGPU(Jz,Jz_gpu,buffJz_gpu)
         end do
    end subroutine smoothen_current_gpu	
	
    subroutine smoothen_current_gpu_Exclusive
         integer :: ni
         do ni=1,curr_filter
              call MovingAverageFilterGPU_Exclusive(Jx_gpu,buffJx_gpu)
              call MovingAverageFilterGPU_Exclusive(Jy_gpu,buffJy_gpu)
              call MovingAverageFilterGPU_Exclusive(Jz_gpu,buffJz_gpu)
         end do
    end subroutine smoothen_current_gpu_Exclusive	
	
     subroutine MovingAverageFilterGPU(Fld,Fld_gpu,FldTemp_gpu)
          real(psn),dimension(mx,my,mz) :: Fld
          real(psn),dimension(mx,my,mz) :: FldTemp
#ifdef twoD			
		  real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu,FldTemp_gpu
#else 			
		  real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu,FldTemp_gpu
#endif 				  
          integer :: i,j,k  
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
		
        call ResetMatrixGPU(FldTemp_gpu)
		FldTemp=0.0_psn
		
		call ExchangeYZEdgeCurrent1(Fld) 
		Fld(xmin1_host:xmin1_host,ymin1_host:ymax1_host-1,z1:z2)=Fld_gpu(xmin1_host:xmin1_host+2,ymin1_host:ymax1_host-1,z1:z2) 
		Fld_gpu(xmin1_host-1:xmin1_host-1,j1:j2,k1:k2)=Fld(xmin1_host-1:xmin1_host-1,j1:j2,k1:k2)
		Fld(xmax1_host-1:xmax1_host-1,ymin1_host:ymax1_host-1,z1:z2)=Fld_gpu(xmax1_host-1:xmax1_host-1,ymin1_host:ymax1_host-1,z1:z2)
		Fld_gpu(xmax1_host:xmax1_host,j1:j2,k1:k2)=Fld(xmax1_host:xmax1_host,j1:j2,k1:k2)
		
		call AvgX_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)	
#ifndef twoD
         do k=3,mz-3
#else
        do k=1,1
#endif
              do j=3,my-3
                  do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*Fld(i-1,j,k)+wt0*Fld(i,j,k)+wtp1*Fld(i+1,j,k)
                    end do
               end do
          end do
		  call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
          Fld=FldTemp
		  
          call ExchangeZXEdgeCurrent1(Fld) 
		  Fld(xmin1_host:xmax1_host-1,ymin1_host:ymin1_host,z1:z2)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymin1_host,z1:z2) 
		  Fld_gpu(i1:i2,ymin1_host-1:ymin1_host-1,k1:k2)=Fld(i1:i2,ymin1_host-1:ymin1_host-1,k1:k2)
		  Fld(xmin1_host:xmax1_host-1,ymax1_host-1:ymax1_host-1,z1:z2)=Fld_gpu(xmin1_host:xmax1_host-1,ymax1_host-1:ymax1_host-1,z1:z2)
		  Fld_gpu(i1:i2,ymax1_host:ymax1_host,k1:k2)=Fld(i1:i2,ymax1_host:ymax1_host,k1:k2)
  		  
		  call AvgY_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)		                     
#ifndef twoD
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                   do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*Fld(i,j-1,k)+wt0*Fld(i,j,k)+wtp1*Fld(i,j+1,k)
                    end do
               end do
          end do
		  call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
          Fld=FldTemp
		  
		  
		  
#ifndef twoD          
          call ExchangeXYEdgeCurrent1(Fld)
		  Fld(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmin1_host:zmin1_host)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmin1_host:zmin1_host) 
		  Fld_gpu(i1:i2,j1:j2,zmin1_host-1:zmin1_host-1)=Fld(i1:i2,j1:j2,zmin1_host-1:zmin1_host-1)
		  Fld(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmax1_host-1:zmax1_host-1)=Fld_gpu(xmin1_host:xmax1_host-1,ymin1_host:ymax1_host-1,zmax1_host-1:zmax1_host-1)
		  Fld_gpu(i1:i2,j1:j2,zmax1_host:zmax1_host)=Fld(i1:i2,j1:j2,zmax1_host:zmax1_host)
		  
		  call AvgZ_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)		                     
          do k=3,mz-3
               do j=3,my-3
                    do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*Fld(i,j,k-1)+wt0*Fld(i,j,k)+wtp1*Fld(i,j,k+1)
                    end do
               end do
          end do		  
		  call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
          Fld=FldTemp
#endif          
     end subroutine MovingAverageFilterGPU
	 
	 
     subroutine MovingAverageFilterGPU_Exclusive(Fld_gpu,FldTemp_gpu)
#ifdef twoD			
		  real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu,FldTemp_gpu
#else 			
		  real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu,FldTemp_gpu
#endif 				  
          integer :: i,j,k  
          integer :: i1,i2,j1,j2,k1,k2,z1,z2
          integer :: dcount1,dcount2,mpi_err
          integer :: stat(MPI_STATUS_SIZE)
		  
		  
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
		
        call ResetMatrixGPU(FldTemp_gpu)
        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     lCurr1_send_gpu,1,my_gpu_host,mz_gpu_host,xmin1_host,ymin1_host-2,zmin1_host-2)
        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     rCurr1_send_gpu,1,my_gpu_host,mz_gpu_host,xmax1_host-1,ymin1_host-2,zmin1_host-2)
		
		 lCurr1_send_host=lCurr1_send_gpu 
		 rCurr1_send_host=rCurr1_send_gpu 
		 dcount1=my_gpu_host*mz_gpu_host
         dcount2=my_gpu_host*mz_gpu_host
         call MPI_SENDRECV(lCurr1_send_host,dcount1,mpi_psn,lproc,1,rCurr1_recv_host,dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
         call MPI_SENDRECV(rCurr1_send_host,dcount2,mpi_psn,rproc,2,lCurr1_recv_host,dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
		 lCurr1_recv_gpu=lCurr1_recv_host
		 rCurr1_recv_gpu=rCurr1_recv_host
		 
         call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     lCurr1_recv_gpu,1,my_gpu_host,mz_gpu_host,xmin1_host-1,ymin1_host-2,zmin1_host-2)
         call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
 			 rCurr1_recv_gpu,1,my_gpu_host,mz_gpu_host,xmax1_host,ymin1_host-2,zmin1_host-2)	
			 
		 call AvgX_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)	
		 call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
		  
         call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     bCurr1_send_gpu,mx_gpu_host,1,mz_gpu_host,xmin1_host-2,ymin1_host,zmin1_host-2)
         call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     tCurr1_send_gpu,mx_gpu_host,1,mz_gpu_host,xmin1_host-2,ymax1_host-1,zmin1_host-2)
		 
		 bCurr1_send_host=bCurr1_send_gpu 
		 tCurr1_send_host=tCurr1_send_gpu 
		 dcount1=mx_gpu_host*mz_gpu_host
         dcount2=mx_gpu_host*mz_gpu_host
         call MPI_SENDRECV(bCurr1_send_host,dcount1,mpi_psn,bproc,1,tCurr1_recv_host,dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
         call MPI_SENDRECV(tCurr1_send_host,dcount2,mpi_psn,tproc,2,bCurr1_recv_host,dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)	
		 bCurr1_recv_gpu=bCurr1_recv_host
		 tCurr1_recv_gpu=tCurr1_recv_host
         
		 call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     bCurr1_recv_gpu,mx_gpu_host,1,mz_gpu_host,xmin1_host-2,ymin1_host-1,zmin1_host-2)
         call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
 			 tCurr1_recv_gpu,mx_gpu_host,1,mz_gpu_host,xmin1_host-2,ymax1_host,zmin1_host-2)	
  		  
		  call AvgY_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)		                     
		  call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
		  
		  
#ifndef twoD 
		  call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		       dCurr1_send_gpu,mx_gpu_host,my_gpu_host,1,xmin1_host-2,ymin1_host-2,zmin1_host)
		  call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		       uCurr1_send_gpu,mx_gpu_host,my_gpu_host,1,xmin1_host-2,ymin1_host-2,zmax1_host-1)    
  		  dCurr1_send_host=dCurr1_send_gpu
  		  uCurr1_send_host=uCurr1_send_gpu 
  		  dcount1=mx_gpu_host*my_gpu_host
          dcount2=mx_gpu_host*my_gpu_host      
          call MPI_SENDRECV(dCurr1_send_host,dcount1,mpi_psn,dproc,1,uCurr1_recv_host,dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(uCurr1_send_host,dcount2,mpi_psn,uproc,2,dCurr1_recv_host,dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
 		  dCurr1_recv_gpu=dCurr1_recv_host
 		  uCurr1_recv_gpu=uCurr1_recv_host
          call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
		     dCurr1_recv_gpu,mx_gpu_host,my_gpu_host,1,xmin1_host-2,ymin1_host-2,zmin1_host-1)
          call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
 			 uCurr1_recv_gpu,mx_gpu_host,my_gpu_host,1,xmin1_host-2,ymin1_host-2,zmax1_host)	
		  
		  call AvgZ_GPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(Fld_gpu,FldTemp_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,wtm1_gpu,wt0_gpu,wtp1_gpu)		                       
		  call CopyMatrixGPU(FldTemp_gpu,Fld_gpu)
#endif          
     end subroutine MovingAverageFilterGPU_Exclusive
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 attributes(global) subroutine AvgX_GPUKernel(Fld1,Fld2,x1,x2,y1,y2,z1,z2,wtm1,wt0,wtp1)
  		integer :: x1,x2,y1,y2,z1,z2
		real  :: wtm1,wt0,wtp1
#ifndef twoD 		
	   real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld1,Fld2
#else 
       real(psn), dimension(x1-2:x2+2,y1-2:y2+2,1) :: Fld1,Fld2
#endif  		
	   integer :: i,j,k
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-1
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-1
#ifndef twoD 		
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-1
       if((i.lt.x2).and.(j.lt.y2).and.(k.lt.z2)) Fld2(i,j,k)=wtm1*Fld1(i-1,j,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i+1,j,k)		     
#else
       k=1
       if((i.lt.x2).and.(j.lt.y2)) Fld2(i,j,k)=wtm1*Fld1(i-1,j,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i+1,j,k)		     
#endif
	 end subroutine AvgX_GPUKernel
	 
	 attributes(global) subroutine AvgY_GPUKernel(Fld1,Fld2,x1,x2,y1,y2,z1,z2,wtm1,wt0,wtp1)
  		integer :: x1,x2,y1,y2,z1,z2
		real  :: wtm1,wt0,wtp1
#ifndef twoD 		
	   real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld1,Fld2
#else 
       real(psn), dimension(x1-2:x2+2,y1-2:y2+2,1) :: Fld1,Fld2
#endif  		
	   integer :: i,j,k
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-1
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-1
#ifndef twoD 		
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-1
       if((i.lt.x2).and.(j.lt.y2).and.(k.lt.z2)) Fld2(i,j,k)=wtm1*Fld1(i,j-1,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j+1,k)     
#else
       k=1
       if((i.lt.x2).and.(j.lt.y2)) Fld2(i,j,k)=wtm1*Fld1(i,j-1,k)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j+1,k) 
#endif
	 end subroutine AvgY_GPUKernel
	 
	 attributes(global) subroutine AvgZ_GPUKernel(Fld1,Fld2,x1,x2,y1,y2,z1,z2,wtm1,wt0,wtp1)
  		integer :: x1,x2,y1,y2,z1,z2
		real  :: wtm1,wt0,wtp1
	    real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld1,Fld2
	    integer :: i,j,k
	
	   i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-1
	   j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-1
	   k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-1
       if((i.lt.x2).and.(j.lt.y2).and.(k.lt.z2)) Fld2(i,j,k)=wtm1*Fld1(i,j,k-1)+wt0*Fld1(i,j,k)+wtp1*Fld1(i,j,k+1)
	 end subroutine AvgZ_GPUKernel
		 	 	 
	 
	 
 	subroutine ResetMatrixGPU(Fld)
#ifdef twoD			
		 real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld
#else 			
		 real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld
#endif 		
 	     type(dim3) :: grid, tBlock

			 
#ifdef twoD			 
 		 tBlock = dim3 (16 ,16 ,1)	
 		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
          tBlock = dim3 (8 ,8 ,4)
          grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
          call ResetMatrixGPUKernel<<<grid,tBlock>>>(Fld,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)			
 	end subroutine ResetMatrixGPU	
 	attributes(global) subroutine ResetMatrixGPUKernel(Fld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld
#else 
        real(psn), dimension(x1-2:x2+2,y1-2:y2+2,1) :: Fld
#endif  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  
         if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) Fld(i,j,k)=0.0
 	end subroutine ResetMatrixGPUKernel
	
	
	
 	subroutine CopyMatrixGPU(FromFld,ToFld)
#ifdef twoD			
		 real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: FromFld,ToFld
#else 			
		 real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: FromFld,ToFld
#endif 		
 	     type(dim3) :: grid, tBlock

			 
#ifdef twoD			 
 		 tBlock = dim3 (16 ,16 ,1)	
 		 grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), 1) 
#else
          tBlock = dim3 (8 ,8 ,4)
          grid = dim3(ceiling(real(xmax1_host-xmin1_host+5)/tBlock%x), ceiling(real(ymax1_host-ymin1_host+5)/tBlock%y), ceiling(real(zmax1_host-zmin1_host+5)/tBlock%z)) 
#endif  
          call CopyMatrixGPUKernel<<<grid,tBlock>>>(FromFld,ToFld,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)			
 	end subroutine CopyMatrixGPU	
 	attributes(global) subroutine CopyMatrixGPUKernel(FromFld,ToFld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real(psn), dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: FromFld,ToFld
#else 
         real(psn), dimension(x1-2:x2+2,y1-2:y2+2,1) :: FromFld,ToFld
#endif  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  
         if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) ToFld(i,j,k)=FromFld(i,j,k)
 	end subroutine CopyMatrixGPUKernel
		
	
end module fields_gpu
	