module bc_gpu
    use parameters
	use vars
	use var_gpu
	implicit none

contains 
	
	!---------------------------------------------------------------------------------
	!  BC  for EM Fld
	!---------------------------------------------------------------------------------
	
	subroutine CondBC_Fld_Top_GPU(yind)
		integer :: yind
		Ex_gpu(1:mx,yind:my,1:mz)=Ex(1:mx,yind:my,1:mz)
		Ez_gpu(1:mx,yind:my,1:mz)=Ez(1:mx,yind:my,1:mz)
	end subroutine CondBC_Fld_Top_GPU
	
	subroutine CondBC_Fld_Bottom_GPU(yind)
		integer :: yind
		Ex_gpu(1:mx,1:yind,1:mz)=Ex(1:mx,1:yind,1:mz)
		Ez_gpu(1:mx,1:yind,1:mz)=Ez(1:mx,1:yind,1:mz)
	end subroutine CondBC_Fld_Bottom_GPU
	
	subroutine CondBC_Fld_Left_GPU(xind)
		integer :: xind
		Ey_gpu(1:xind,1:my,1:mz)=Ey(1:xind,1:my,1:mz)
		Ez_gpu(1:xind,1:my,1:mz)=Ez(1:xind,1:my,1:mz)
	end subroutine CondBC_Fld_Left_GPU

	
	!---------------------------------------------------------------------------------
	!  BC  for Prtl
	!---------------------------------------------------------------------------------
	
	!---- Reflecting boundary condition
	subroutine RefBC_Prtl_Top_GPU(ymax_local,PrtlWall_local)
		real(psn) :: ymax_local
		integer :: PrtlWall_local
		integer :: kc,indi,indf
	
	
		if(ymax_local.lt.my-1) then 
			
			do kc=1,Nchunk_prtl_gpu
			    indi=(kc-1)*chunk_size_prtl_gpu+1
			    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			    call RefBC_Prtl_Top_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,yp_gpu,vp_gpu,ymax_local,indi,indf)
		    end do 
	
		end if
	
		if(PrtlWall_local.ge.2.and.PrtlWall_local.le.my-2) then
			
			Jx(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jx_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
			Jy(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jy_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
			Jz(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jz_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
			
			Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local+1,:)
		 	Jy(:,PrtlWall_local-1,:)=Jy(:,PrtlWall_local-1,:)-Jy(:,PrtlWall_local,:)
		 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local+1,:)
			Jx(:,PrtlWall_local+1,:)=0.0_psn
			Jy(:,PrtlWall_local,:)=0.0_psn
			Jz(:,PrtlWall_local+1,:)=0.0_psn
			
			Jx_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jx(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
			Jy_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jy(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
			Jz_gpu(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz) = Jz(1:mx,PrtlWall_local-1:PrtlWall_local+1,1:mz)
	    
		end if
	
	end subroutine RefBC_Prtl_Top_GPU
	
	attributes(global) subroutine RefBC_Prtl_Top_Kernel(qp,yp,vp,ymax_local,indi,indf)
		real, dimension(:) :: qp,yp,vp
		real, value :: ymax_local
		integer, value :: indi, indf
		integer :: n
		
	    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
	    if(n.gt.indf) return
		
		if((qp(n).ne.0).and.(yp(n).gt.ymax_local)) then 
			 vp(n)=-vp(n)
             yp(n)=ymax_local-(yp(n)-ymax_local)
		end if 
	
	end subroutine RefBC_Prtl_Top_Kernel
	
	subroutine RefBC_Prtl_Bottom_GPU(ymin_local,PrtlWall_local)
		real(psn) :: ymin_local
		integer :: PrtlWall_local
		integer :: kc, indi, indf
	
		if(ymin_local.gt.1) then
			
			do kc=1,Nchunk_prtl_gpu
			    indi=(kc-1)*chunk_size_prtl_gpu+1
			    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			    call RefBC_Prtl_Left_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,yp_gpu,vp_gpu,ymin_local,indi,indf)
		    end do 
			
		end if
	
		if(PrtlWall_local.ge.3.and.PrtlWall_local.le.my-2) then
			
			Jx(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jx_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
			Jy(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jy_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
			Jz(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jz_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
			
			Jx(:,PrtlWall_local,:)=Jx(:,PrtlWall_local,:)+Jx(:,PrtlWall_local-1,:)
		 	Jy(:,PrtlWall_local,:)=Jy(:,PrtlWall_local,:)-Jy(:,PrtlWall_local-1,:)
		 	Jz(:,PrtlWall_local,:)=Jz(:,PrtlWall_local,:)+Jz(:,PrtlWall_local-1,:)
			Jx(:,PrtlWall_local-1,:)=0.0_psn
			Jy(:,PrtlWall_local-1,:)=0.0_psn
			Jz(:,PrtlWall_local-1,:)=0.0_psn
			
			Jx_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jx(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
			Jy_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jy(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
			Jz_gpu(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz) = Jz(1:mx,PrtlWall_local-1:PrtlWall_local,1:mz)
		
		end if
	end subroutine RefBC_Prtl_Bottom_GPU
	
	
	subroutine RefBC_Prtl_Left_GPU(xmin_local,PrtlWall_local)
		real(psn) :: xmin_local
		integer :: PrtlWall_local
		integer :: kc, indi, indf
	
		if(xmin_local.gt.1) then
			
			do kc=1,Nchunk_prtl_gpu
			    indi=(kc-1)*chunk_size_prtl_gpu+1
			    indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			    call RefBC_Prtl_Left_Kernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,up_gpu,xmin_local,indi,indf)
		    end do 
			
		end if
	
		if(PrtlWall_local.ge.3.and.PrtlWall_local.le.mx-2) then
			
			Jx(PrtlWall_local:PrtlWall_local,1:my,1:mz) = Jx_gpu(PrtlWall_local:PrtlWall_local,1:my,1:mz)
			Jy(PrtlWall_local:PrtlWall_local+1,1:my,1:mz) = Jy_gpu(PrtlWall_local:PrtlWall_local+1,1:my,1:mz)
			Jz(PrtlWall_local:PrtlWall_local+1,1:my,1:mz) = Jz_gpu(PrtlWall_local:PrtlWall_local+1,1:my,1:mz)
			
			Jy(PrtlWall_local+1,:,:)=Jy(PrtlWall_local+1,:,:)+Jy(PrtlWall_local,:,:)
		 	Jx(PrtlWall_local,:,:)=0
		 	Jz(PrtlWall_local+1,:,:)=Jz(PrtlWall_local+1,:,:)+Jz(PrtlWall_local,:,:)
			Jy(PrtlWall_local,:,:)=0.0_psn
			Jz(PrtlWall_local,:,:)=0.0_psn
			
			Jx_gpu(PrtlWall_local:PrtlWall_local,1:my,1:mz) = Jx(PrtlWall_local:PrtlWall_local,1:my,1:mz)
			Jy_gpu(PrtlWall_local:PrtlWall_local+1,1:my,1:mz) = Jy(PrtlWall_local:PrtlWall_local+1,1:my,1:mz)
			Jz_gpu(PrtlWall_local:PrtlWall_local+1,1:my,1:mz) = Jz(PrtlWall_local:PrtlWall_local+1,1:my,1:mz)
			
			
! 			Jx(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jx_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
! 			Jy(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jy_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
! 			Jz(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jz_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
!
! 			Jx(PrtlWall_local,:,:)=Jx(PrtlWall_local,:,:)-Jx(PrtlWall_local-1,:,:)
! 		 	Jy(PrtlWall_local,:,:)=Jy(PrtlWall_local,:,:)+Jy(PrtlWall_local-1,:,:)
! 		 	Jz(PrtlWall_local,:,:)=Jz(PrtlWall_local,:,:)+Jz(PrtlWall_local-1,:,:)
! 			Jx(PrtlWall_local-1,:,:)=0.0_psn
! 			Jy(PrtlWall_local-1,:,:)=0.0_psn
! 			Jz(PrtlWall_local-1,:,:)=0.0_psn
!
! 			Jx_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jx(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
! 			Jy_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jy(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
! 			Jz_gpu(PrtlWall_local-1:PrtlWall_local,1:my,1:mz) = Jz(PrtlWall_local-1:PrtlWall_local,1:my,1:mz)
		
		end if
	end subroutine RefBC_Prtl_Left_GPU
	
	attributes(global) subroutine RefBC_Prtl_Left_Kernel(qp,xp,up,xmin_local,indi,indf)
		real, dimension(:) :: qp,xp,up
		real, value :: xmin_local
		integer, value :: indi, indf
		integer :: n
		
	    n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
	    if(n.gt.indf) return
		
		if((qp(n).ne.0).and.(xp(n).lt.xmin_local)) then 
			up(n)=-up(n)
			xp(n)=xmin_local+(xmin_local-xp(n))
		end if 
	
	end subroutine RefBC_Prtl_Left_Kernel
	
	
end module bc_gpu