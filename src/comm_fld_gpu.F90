module comm_fld_gpu 
	use parameters
	use vars 
	use var_gpu 
	use memory
	use cudafor
	use communication

contains
!----------------------------------------------------------------------------------------------------------	
! The following subroutines are used in packing and exchanging data between CPU and GPU 
!
!----------------------------------------------------------------------------------------------------------	
   subroutine EMfldExchangeGPU_Exclusive
	   call FldYZEdgeExchangeGPU_Exclusive(Ex,lEx_send_host,rEx_send_host,lEx_recv_host,rEx_recv_host,Ex_gpu,lEx_send_gpu,rEx_send_gpu,lEx_recv_gpu,rEx_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Ey,lEy_send_host,rEy_send_host,lEy_recv_host,rEy_recv_host,Ey_gpu,lEy_send_gpu,rEy_send_gpu,lEy_recv_gpu,rEy_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Ez,lEz_send_host,rEz_send_host,lEz_recv_host,rEz_recv_host,Ez_gpu,lEz_send_gpu,rEz_send_gpu,lEz_recv_gpu,rEz_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Bx,lBx_send_host,rBx_send_host,lBx_recv_host,rBx_recv_host,Bx_gpu,lBx_send_gpu,rBx_send_gpu,lBx_recv_gpu,rBx_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(By,lBy_send_host,rBy_send_host,lBy_recv_host,rBy_recv_host,By_gpu,lBy_send_gpu,rBy_send_gpu,lBy_recv_gpu,rBy_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Bz,lBz_send_host,rBz_send_host,lBz_recv_host,rBz_recv_host,Bz_gpu,lBz_send_gpu,rBz_send_gpu,lBz_recv_gpu,rBz_recv_gpu)
	   
	   call FldZXEdgeExchangeGPU_Exclusive(Ex,bEx_send_host,tEx_send_host,bEx_recv_host,tEx_recv_host,Ex_gpu,bEx_send_gpu,tEx_send_gpu,bEx_recv_gpu,tEx_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Ey,bEy_send_host,tEy_send_host,bEy_recv_host,tEy_recv_host,Ey_gpu,bEy_send_gpu,tEy_send_gpu,bEy_recv_gpu,tEy_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Ez,bEz_send_host,tEz_send_host,bEz_recv_host,tEz_recv_host,Ez_gpu,bEz_send_gpu,tEz_send_gpu,bEz_recv_gpu,tEz_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Bx,bBx_send_host,tBx_send_host,bBx_recv_host,tBx_recv_host,Bx_gpu,bBx_send_gpu,tBx_send_gpu,bBx_recv_gpu,tBx_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(By,bBy_send_host,tBy_send_host,bBy_recv_host,tBy_recv_host,By_gpu,bBy_send_gpu,tBy_send_gpu,bBy_recv_gpu,tBy_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Bz,bBz_send_host,tBz_send_host,bBz_recv_host,tBz_recv_host,Bz_gpu,bBz_send_gpu,tBz_send_gpu,bBz_recv_gpu,tBz_recv_gpu)

#ifndef twoD  	   
	   call FldXYEdgeExchangeGPU_Exclusive(Ex,dEx_send_host,uEx_send_host,dEx_recv_host,uEx_recv_host,Ex_gpu,dEx_send_gpu,uEx_send_gpu,dEx_recv_gpu,uEx_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Ey,dEy_send_host,uEy_send_host,dEy_recv_host,uEy_recv_host,Ey_gpu,dEy_send_gpu,uEy_send_gpu,dEy_recv_gpu,uEy_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Ez,dEz_send_host,uEz_send_host,dEz_recv_host,uEz_recv_host,Ez_gpu,dEz_send_gpu,uEz_send_gpu,dEz_recv_gpu,uEz_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Bx,dBx_send_host,uBx_send_host,dBx_recv_host,uBx_recv_host,Bx_gpu,dBx_send_gpu,uBx_send_gpu,dBx_recv_gpu,uBx_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(By,dBy_send_host,uBy_send_host,dBy_recv_host,uBy_recv_host,By_gpu,dBy_send_gpu,uBy_send_gpu,dBy_recv_gpu,uBy_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Bz,dBz_send_host,uBz_send_host,dBz_recv_host,uBz_recv_host,Bz_gpu,dBz_send_gpu,uBz_send_gpu,dBz_recv_gpu,uBz_recv_gpu)
#endif 	   
	   
   end subroutine EMfldExchangeGPU_Exclusive 
   
   subroutine UpdateCurrentsAllEdgesGPU_Exclusive
	   call CurrentYZEdgeUpdateGPU_Exclusive(Jx,lJx_send_host,rJx_send_host,lJx_recv_host,rJx_recv_host,Jx_gpu,lJx_send_gpu,rJx_send_gpu,lJx_recv_gpu,rJx_recv_gpu)
	   call CurrentYZEdgeUpdateGPU_Exclusive(Jy,lJy_send_host,rJy_send_host,lJy_recv_host,rJy_recv_host,Jy_gpu,lJy_send_gpu,rJy_send_gpu,lJy_recv_gpu,rJy_recv_gpu)
	   call CurrentYZEdgeUpdateGPU_Exclusive(Jz,lJz_send_host,rJz_send_host,lJz_recv_host,rJz_recv_host,Jz_gpu,lJz_send_gpu,rJz_send_gpu,lJz_recv_gpu,rJz_recv_gpu)
	   
	   call CurrentZXEdgeUpdateGPU_Exclusive(Jx,bJx_send_host,tJx_send_host,bJx_recv_host,tJx_recv_host,Jx_gpu,bJx_send_gpu,tJx_send_gpu,bJx_recv_gpu,tJx_recv_gpu)
	   call CurrentZXEdgeUpdateGPU_Exclusive(Jy,bJy_send_host,tJy_send_host,bJy_recv_host,tJy_recv_host,Jy_gpu,bJy_send_gpu,tJy_send_gpu,bJy_recv_gpu,tJy_recv_gpu)
	   call CurrentZXEdgeUpdateGPU_Exclusive(Jz,bJz_send_host,tJz_send_host,bJz_recv_host,tJz_recv_host,Jz_gpu,bJz_send_gpu,tJz_send_gpu,bJz_recv_gpu,tJz_recv_gpu)

#ifndef twoD  	   
	   call CurrentXYEdgeUpdateGPU_Exclusive(Jx,dJx_send_host,uJx_send_host,dJx_recv_host,uJx_recv_host,Jx_gpu,dJx_send_gpu,uJx_send_gpu,dJx_recv_gpu,uJx_recv_gpu)
	   call CurrentXYEdgeUpdateGPU_Exclusive(Jy,dJy_send_host,uJy_send_host,dJy_recv_host,uJy_recv_host,Jy_gpu,dJy_send_gpu,uJy_send_gpu,dJy_recv_gpu,uJy_recv_gpu)
	   call CurrentXYEdgeUpdateGPU_Exclusive(Jz,dJz_send_host,uJz_send_host,dJz_recv_host,uJz_recv_host,Jz_gpu,dJz_send_gpu,uJz_send_gpu,dJz_recv_gpu,uJz_recv_gpu)
#endif 	   
   end subroutine UpdateCurrentsAllEdgesGPU_Exclusive
   
   subroutine SyncCurrentEdgesGPU_Exclusive
	   call FldYZEdgeExchangeGPU_Exclusive(Jx,lEx_send_host,rEx_send_host,lEx_recv_host,rEx_recv_host,Jx_gpu,lEx_send_gpu,rEx_send_gpu,lEx_recv_gpu,rEx_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Jy,lEy_send_host,rEy_send_host,lEy_recv_host,rEy_recv_host,Jy_gpu,lEy_send_gpu,rEy_send_gpu,lEy_recv_gpu,rEy_recv_gpu)
	   call FldYZEdgeExchangeGPU_Exclusive(Jz,lEz_send_host,rEz_send_host,lEz_recv_host,rEz_recv_host,Jz_gpu,lEz_send_gpu,rEz_send_gpu,lEz_recv_gpu,rEz_recv_gpu)
	   
	   call FldZXEdgeExchangeGPU_Exclusive(Jx,bEx_send_host,tEx_send_host,bEx_recv_host,tEx_recv_host,Jx_gpu,bEx_send_gpu,tEx_send_gpu,bEx_recv_gpu,tEx_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Jy,bEy_send_host,tEy_send_host,bEy_recv_host,tEy_recv_host,Jy_gpu,bEy_send_gpu,tEy_send_gpu,bEy_recv_gpu,tEy_recv_gpu)
	   call FldZXEdgeExchangeGPU_Exclusive(Jz,bEz_send_host,tEz_send_host,bEz_recv_host,tEz_recv_host,Jz_gpu,bEz_send_gpu,tEz_send_gpu,bEz_recv_gpu,tEz_recv_gpu)

#ifndef twoD  	   
	   call FldXYEdgeExchangeGPU_Exclusive(Jx,dEx_send_host,uEx_send_host,dEx_recv_host,uEx_recv_host,Jx_gpu,dEx_send_gpu,uEx_send_gpu,dEx_recv_gpu,uEx_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Jy,dEy_send_host,uEy_send_host,dEy_recv_host,uEy_recv_host,Jy_gpu,dEy_send_gpu,uEy_send_gpu,dEy_recv_gpu,uEy_recv_gpu)
	   call FldXYEdgeExchangeGPU_Exclusive(Jz,dEz_send_host,uEz_send_host,dEz_recv_host,uEz_recv_host,Jz_gpu,dEz_send_gpu,uEz_send_gpu,dEz_recv_gpu,uEz_recv_gpu)
#endif 
   end subroutine SyncCurrentEdgesGPU_Exclusive
   
   
   
   subroutine FldYZEdgeExchangeGPU_Exclusive(Fld,lFld_send,rFld_send,lFld_recv,rFld_recv,Fld_gpu,lFld_send_gpu,rFld_send_gpu,lFld_recv_gpu,rFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host)         :: lFld_send,rFld_recv
   			real(psn), dimension(2,my_gpu_host,mz_gpu_host)         :: rFld_send,lFld_recv
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host), device :: lFld_send_gpu,rFld_recv_gpu
   			real(psn), dimension(2,my_gpu_host,mz_gpu_host), device :: rFld_send_gpu,lFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     lFld_send_gpu,3,my_gpu_host,mz_gpu_host,xmin1_host,ymin1_host-2,zmin1_host-2)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     rFld_send_gpu,2,my_gpu_host,mz_gpu_host,xmax1_host-2,ymin1_host-2,zmin1_host-2)
			lFld_send=lFld_send_gpu
			rFld_send=rFld_send_gpu
			
            dcount1=3*my_gpu_host*mz_gpu_host
            dcount2=2*my_gpu_host*mz_gpu_host
            call MPI_SENDRECV(lFld_send,dcount1,mpi_psn,lproc,1,rFld_recv,dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(rFld_send,dcount2,mpi_psn,rproc,2,lFld_recv,dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			lFld_recv_gpu=lFld_recv
			rFld_recv_gpu=rFld_recv	
            call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     lFld_recv_gpu,2,my_gpu_host,mz_gpu_host,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 rFld_recv_gpu,3,my_gpu_host,mz_gpu_host,xmax1_host,ymin1_host-2,zmin1_host-2)						 		 	
   end subroutine FldYZEdgeExchangeGPU_Exclusive
   
   subroutine FldZXEdgeExchangeGPU_Exclusive(Fld,bFld_send,tFld_send,bFld_recv,tFld_recv,Fld_gpu,bFld_send_gpu,tFld_send_gpu,bFld_recv_gpu,tFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host)         :: bFld_send,tFld_recv
   			real(psn), dimension(mx_gpu_host,2,mz_gpu_host)         :: tFld_send,bFld_recv
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host), device :: bFld_send_gpu,tFld_recv_gpu
   			real(psn), dimension(mx_gpu_host,2,mz_gpu_host), device :: tFld_send_gpu,bFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     bFld_send_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymin1_host,zmin1_host-2)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     tFld_send_gpu,mx_gpu_host,2,mz_gpu_host,xmin1_host-2,ymax1_host-2,zmin1_host-2)
			bFld_send=bFld_send_gpu
			tFld_send=tFld_send_gpu
			
            dcount1=3*mx_gpu_host*mz_gpu_host
            dcount2=2*mx_gpu_host*mz_gpu_host
            call MPI_SENDRECV(bFld_send,dcount1,mpi_psn,bproc,1,tFld_recv,dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(tFld_send,dcount2,mpi_psn,tproc,2,bFld_recv,dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			bFld_recv_gpu=bFld_recv
			tFld_recv_gpu=tFld_recv	
            call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     bFld_recv_gpu,mx_gpu_host,2,mz_gpu_host,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 tFld_recv_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymax1_host,zmin1_host-2)					 		 	
   end subroutine FldZXEdgeExchangeGPU_Exclusive
   
   subroutine FldXYEdgeExchangeGPU_Exclusive(Fld,dFld_send,uFld_send,dFld_recv,uFld_recv,Fld_gpu,dFld_send_gpu,uFld_send_gpu,dFld_recv_gpu,uFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3)         :: dFld_send,uFld_recv
   			real(psn), dimension(mx_gpu_host,my_gpu_host,2)         :: uFld_send,dFld_recv
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3), device :: dFld_send_gpu,uFld_recv_gpu
   			real(psn), dimension(mx_gpu_host,my_gpu_host,2), device :: uFld_send_gpu,dFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     dFld_send_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmin1_host)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     uFld_send_gpu,mx_gpu_host,my_gpu_host,2,xmin1_host-2,ymin1_host-2,zmax1_host-2)
			dFld_send=dFld_send_gpu
			uFld_send=uFld_send_gpu
			
            dcount1=3*mx_gpu_host*my_gpu_host
            dcount2=2*mx_gpu_host*my_gpu_host
            call MPI_SENDRECV(dFld_send,dcount1,mpi_psn,dproc,1,uFld_recv,dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(uFld_send,dcount2,mpi_psn,uproc,2,dFld_recv,dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			dFld_recv_gpu=dFld_recv
			uFld_recv_gpu=uFld_recv	
            call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     dFld_recv_gpu,mx_gpu_host,my_gpu_host,2,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyFromTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 uFld_recv_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmax1_host)					 		 	
   end subroutine FldXYEdgeExchangeGPU_Exclusive
   
   
   
   subroutine CurrentYZEdgeUpdateGPU_Exclusive(Fld,lFld_send,rFld_send,lFld_recv,rFld_recv,Fld_gpu,lFld_send_gpu,rFld_send_gpu,lFld_recv_gpu,rFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host)         :: lFld_send,rFld_recv
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host)         :: rFld_send,lFld_recv
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host), device :: lFld_send_gpu,rFld_recv_gpu
   			real(psn), dimension(3,my_gpu_host,mz_gpu_host), device :: rFld_send_gpu,lFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     lFld_send_gpu,3,my_gpu_host,mz_gpu_host,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     rFld_send_gpu,3,my_gpu_host,mz_gpu_host,xmax1_host,ymin1_host-2,zmin1_host-2)
			lFld_send=lFld_send_gpu
			rFld_send=rFld_send_gpu
			
            dcount1=3*my_gpu_host*mz_gpu_host
            dcount2=3*my_gpu_host*mz_gpu_host
            call MPI_SENDRECV(lFld_send,dcount1,mpi_psn,lproc,1,rFld_recv,dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(rFld_send,dcount2,mpi_psn,rproc,2,lFld_recv,dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			lFld_recv_gpu=lFld_recv
			rFld_recv_gpu=rFld_recv	
            call AddTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     lFld_recv_gpu,3,my_gpu_host,mz_gpu_host,xmin1_host,ymin1_host-2,zmin1_host-2)
	        call AddTransferMatrixGPUKernel<<<tGrid_gpu_YZedge,tBlock_gpu_YZedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 rFld_recv_gpu,3,my_gpu_host,mz_gpu_host,xmax1_host-2,ymin1_host-2,zmin1_host-2)						 		 	
   end subroutine CurrentYZEdgeUpdateGPU_Exclusive
   
   subroutine CurrentZXEdgeUpdateGPU_Exclusive(Fld,bFld_send,tFld_send,bFld_recv,tFld_recv,Fld_gpu,bFld_send_gpu,tFld_send_gpu,bFld_recv_gpu,tFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host)         :: bFld_send,tFld_recv
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host)         :: tFld_send,bFld_recv
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host), device :: bFld_send_gpu,tFld_recv_gpu
   			real(psn), dimension(mx_gpu_host,3,mz_gpu_host), device :: tFld_send_gpu,bFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     bFld_send_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     tFld_send_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymax1_host,zmin1_host-2)
			bFld_send=bFld_send_gpu
			tFld_send=tFld_send_gpu
			
            dcount1=3*mx_gpu_host*mz_gpu_host
            dcount2=3*mx_gpu_host*mz_gpu_host
            call MPI_SENDRECV(bFld_send,dcount1,mpi_psn,bproc,1,tFld_recv,dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(tFld_send,dcount2,mpi_psn,tproc,2,bFld_recv,dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			bFld_recv_gpu=bFld_recv
			tFld_recv_gpu=tFld_recv	
            call AddTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     bFld_recv_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymin1_host,zmin1_host-2)
	        call AddTransferMatrixGPUKernel<<<tGrid_gpu_ZXedge,tBlock_gpu_ZXedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 tFld_recv_gpu,mx_gpu_host,3,mz_gpu_host,xmin1_host-2,ymax1_host-2,zmin1_host-2)					 		 	
   end subroutine CurrentZXEdgeUpdateGPU_Exclusive
   
   subroutine CurrentXYEdgeUpdateGPU_Exclusive(Fld,dFld_send,uFld_send,dFld_recv,uFld_recv,Fld_gpu,dFld_send_gpu,uFld_send_gpu,dFld_recv_gpu,uFld_recv_gpu)
   			real(psn), dimension(mx,my,mz) :: Fld
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3)         :: dFld_send,uFld_recv
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3)         :: uFld_send,dFld_recv
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3), device :: dFld_send_gpu,uFld_recv_gpu
   			real(psn), dimension(mx_gpu_host,my_gpu_host,3), device :: uFld_send_gpu,dFld_recv_gpu
            integer :: dcount1,dcount2,mpi_err
            integer :: stat(MPI_STATUS_SIZE)
			
#ifdef twoD			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,1:1), device :: Fld_gpu
#else 			
   			real(psn), dimension(xmin1_host-2:xmax1_host+2,ymin1_host-2:ymax1_host+2,zmin1_host-2:zmax1_host+2), device :: Fld_gpu
#endif
            call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     dFld_send_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmin1_host-2)
	        call CopyToTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     uFld_send_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmax1_host)
			dFld_send=dFld_send_gpu
			uFld_send=uFld_send_gpu
			
            dcount1=3*mx_gpu_host*my_gpu_host
            dcount2=3*mx_gpu_host*my_gpu_host
            call MPI_SENDRECV(dFld_send,dcount1,mpi_psn,dproc,1,uFld_recv,dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
            call MPI_SENDRECV(uFld_send,dcount2,mpi_psn,uproc,2,dFld_recv,dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
            
			dFld_recv_gpu=dFld_recv
			uFld_recv_gpu=uFld_recv	
            call AddTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
			     dFld_recv_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmin1_host)
	        call AddTransferMatrixGPUKernel<<<tGrid_gpu_XYedge,tBlock_gpu_XYedge>>>(Fld_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,&
	 			 uFld_recv_gpu,mx_gpu_host,my_gpu_host,3,xmin1_host-2,ymin1_host-2,zmax1_host-2)					 		 	
   end subroutine CurrentXYEdgeUpdateGPU_Exclusive
   
   
   
   
   

	
	subroutine CopyToTransferMatrixCPU(Fld,Fld_buff,sx,sy,sz,off1,off2,off3)
		integer :: sx,sy,sz
		integer :: off1,off2,off3
		real(psn), dimension(mx,my,mz) :: Fld
		real(psn), dimension(sx,sy,sz) :: Fld_buff		
        integer :: i,j,k 

        do k=1,sz
			do j=1,sy
				do i=1,sx
		            Fld_buff(i,j,k)=Fld(i+off1-1,j+off2-1,k+off3-1) 
				end do  
			end do 
		end do 
	end subroutine CopyToTransferMatrixCPU
	
	subroutine CopyFromTransferMatrixCPU(Fld,Fld_buff,sx,sy,sz,off1,off2,off3)
		integer :: sx,sy,sz
		integer :: off1,off2,off3
		real(psn), dimension(mx,my,mz) :: Fld
		real(psn), dimension(sx,sy,sz) :: Fld_buff		
        integer :: i,j,k 

        do k=1,sz
			do j=1,sy 
				do i=1,sz
		            Fld(i+off1-1,j+off2-1,k+off3-1)=Fld_buff(i,j,k) 
				end do  
			end do 
		end do 
	end subroutine CopyFromTransferMatrixCPU
	attributes(global) subroutine CopyToTransferMatrixGPUKernel(Fld,x1,x2,y1,y2,z1,z2,Fld_buff,sx,sy,sz,off1,off2,off3)
		integer :: x1,x2,y1,y2,z1,z2
		integer,value :: sx,sy,sz
		integer, value :: off1,off2,off3
		real, dimension(sx,sy,sz) :: Fld_buff
#ifndef twoD 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1:1) :: Fld
#endif 
        integer :: i,j,k 

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        if((i.le.sx).and.(j.le.sy).and.(k.le.sz)) Fld_buff(i,j,k)=Fld(i+off1-1,j+off2-1,k+off3-1) 
#else
        k=1
        if((i.le.sx).and.(j.le.sy)) Fld_buff(i,j,1)=Fld(i+off1-1,j+off2-1,1) 
#endif 
	end subroutine CopyToTransferMatrixGPUKernel 
	
	attributes(global) subroutine CopyFromTransferMatrixGPUKernel(Fld,x1,x2,y1,y2,z1,z2,Fld_buff,sx,sy,sz,off1,off2,off3)
		integer :: x1,x2,y1,y2,z1,z2
		integer,value :: sx,sy,sz
		integer, value :: off1,off2,off3
		real, dimension(sx,sy,sz) :: Fld_buff
#ifndef twoD 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1:1) :: Fld
#endif 
        integer :: i,j,k 

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        if((i.le.sx).and.(j.le.sy).and.(k.le.sz)) Fld(i+off1-1,j+off2-1,k+off3-1)=Fld_buff(i,j,k)
#else
        k=1
        if((i.le.sx).and.(j.le.sy)) Fld(i+off1-1,j+off2-1,1)=Fld_buff(i,j,1)
#endif 
	end subroutine CopyFromTransferMatrixGPUKernel 
	
	attributes(global) subroutine AddTransferMatrixGPUKernel(Fld,x1,x2,y1,y2,z1,z2,Fld_buff,sx,sy,sz,off1,off2,off3)
		integer :: x1,x2,y1,y2,z1,z2
		integer,value :: sx,sy,sz
		integer, value :: off1,off2,off3
		real, dimension(sx,sy,sz) :: Fld_buff
#ifndef twoD 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1:1) :: Fld
#endif 
        integer :: i,j,k 

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
        if((i.le.sx).and.(j.le.sy).and.(k.le.sz)) Fld(i+off1-1,j+off2-1,k+off3-1)=Fld(i+off1-1,j+off2-1,k+off3-1)+Fld_buff(i,j,k)
#else
        k=1
        if((i.le.sx).and.(j.le.sy)) Fld(i+off1-1,j+off2-1,1)=Fld(i+off1-1,j+off2-1,1)+Fld_buff(i,j,1)
#endif 
	end subroutine AddTransferMatrixGPUKernel 
	
end module comm_fld_gpu	