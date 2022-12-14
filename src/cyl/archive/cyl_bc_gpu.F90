module cyl_bc_gpu
    use parameters
    use vars
	use var_gpu
	use cyl_vars
	implicit none 
contains 
!---------------------------------------------------------------------------------
!  BC  for EM Fld
!---------------------------------------------------------------------------------

	subroutine ConductingOuterBC_Fld(rmax)
		real(psn) :: rmax
		integer   :: rmax_local
		integer  :: i
		rmax_local=rmax-rborders(procxind)+3
		rmax_local=max(rmax_local,1)
		if(rmax_local.le.mx) then 
			call SetFldToZero_Xrange_GPU<<<grid,tBlock>>>(Ey_gpu,mx_gpu,my_gpu,mz_gpu,rmax_local,mx)
			call SetFldToZero_Xrange_GPU<<<grid,tBlock>>>(Ez_gpu,mx_gpu,my_gpu,mz_gpu,rmax_local,mx)
		end if 
	end subroutine ConductingOuterBC_Fld

	subroutine ConductingInnerBC_Fld(rmin) !Note: place the conducting wall a few cells away from the particle wall
		real(psn) :: rmin
		integer   :: rmin_local
		rmin_local=rmin-rborders(procxind)+3
		rmin_local=min(rmin_local,mx)
		if(rmin_local.ge.1) then 
			call SetFldToZero_Xrange_GPU<<<grid,tBlock>>>(Ey_gpu,mx_gpu,my_gpu,mz_gpu,1,rmin_local)
			call SetFldToZero_Xrange_GPU<<<grid,tBlock>>>(Ez_gpu,mx_gpu,my_gpu,mz_gpu,1,rmin_local)
		end if 
	end subroutine ConductingInnerBC_Fld
	
	attributes(global) subroutine SetFldToZero_Xrange_GPU(Fld,mx,my,mz,i1,i2)
	    integer :: mx,my,mz
	    real, dimension(mx,my,mz) :: Fld
		integer, value :: i1,i2 
		integer        :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 
        if(i.ge.i1.and.i.le.i2) then 
			if(j.le.my.and.k.le.mz) then
			    Fld(i,j,k)=0
		    end if 
		 end if 
    end subroutine SetFldToZero_Xrange_GPU
	
	attributes(global) subroutine SetFldToZero_Xindex_GPU(Fld,mx,my,mz,ind)
	    integer :: mx,my,mz
		real, dimension(mx,my,mz) :: Fld
		integer, value :: ind 
		integer        :: j,k
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 
		if(j.le.my.and.k.le.mz) Fld(ind,j,k)=0
    end subroutine SetFldToZero_Xindex_GPU

	subroutine AxisCurrentBC
		 if(procxind.ne.0) return 
         call AxisCurrentBC_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jx_gpu,Jy_gpu,Jz_gpu,mx,my,mz,ny)
	end subroutine AxisCurrentBC

    attributes(global) subroutine AxisCurrentBC_GPU(Jx,Jy,Jz,mx,my,mz,ny)
	    integer, value :: mx,my,mz
		integer, value :: ny
		integer :: j,k,jp,jm,jpp 
	    real, dimension(mx,my,mz) :: Jx,Jy,Jz
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 	
        if(j.le.my.and.k.le.mz) then 
			
! 			jp=j + ny/4
! 			jm=j - ny/4
! 			jpp = j + ny/2
! 			if(jp.gt.my-3) jp= jp - ny
! 			if(jpp.gt.my-3) jpp= jpp - ny
! 			if(jm.lt.3) jm = jm + ny
!
! 			Jy(3,jp,k) = Jy(3,jp,k) - 0.5_psn*Jx(3,j,k)
! 			Jy(3,jm,k) = Jy(3,jm,k) + 0.5_psn*Jx(3,j,k)
!
! 			Jy(3,jpp,k) = Jy(3,jpp,k) - Jy(2,jpp,k)
! 			Jz(3,jpp,k) = Jz(3,jpp,k) + Jz(2,jpp,k)

			Jy(4,j,k)=Jy(4,j,k)-Jy(3,j,k)
			Jz(4,j,k)=Jz(4,j,k)-Jz(3,j,k)
			
			Jx(3,j,k) = 0.0_psn
			Jy(2,j,k) = 0.0_psn
			Jz(2,j,k) = 0.0_psn
			
	    end if 
    end subroutine AxisCurrentBC_GPU 
	
!--------------------------------------------------------------
! Update all EM Flds at the Axis
!--------------------------------------------------------------	
! 	subroutine UpdateFldAxis
! 		if((inc_axis.eq..true.).and.(procxind.eq.0)) then
! 		     call UpdateFldAxisGPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,mx,my,mz)
! 		end if
! 	end subroutine UpdateFldAxis
!
! 	attributes(global) subroutine UpdateFldAxisGPU(Ex,Ey,Ez,Bx,By,Bz,mx,my,mz)
! 		integer, value :: mx,my,mz
! 		real, dimension(mx,my,mz) :: Ex,Ey,Ez,Bx,By,Bz
! 		integer :: j,k
!
! 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
! #ifndef twoD
! 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
! #else
!         k=1
! #endif
!
! 		if((j.le.my).and.(k.le.mz)) then
! 	   		 Ex(3,j,k)=0.0
! 	   		 Ex(2,j,k)=-Ex(4,j,k)
!
! 	   		 Ey(3,j,k)=-Ey(4,j,k)
! 	   		 Ez(3,j,k)=-Ez(4,j,k)
!
! 	   		 Bx(3,j,k)=-Bx(4,j,k)
! 	   		 By(3,j,k)=0.0
! 	   		 By(2,j,k)=-By(4,j,k)
! 	   		 !Bz(2,j,k)=2*Bz(3,j,k)-Bz(4,j,k)
! 			 Bz(2,j,k)=Bz(3,j,k)
! 		end if
!     end subroutine UpdateFldAxisGPU

!---------------------------------------------------------------------------------
!  BC  for Prtl
!---------------------------------------------------------------------------------
	subroutine RefInnerBC_Prtl(rmin) 
		real(psn) :: rmin,rmin_local,f
		integer   :: ind,n,kc
		integer   :: indi,indf
		
	    rmin_local=rmin-rborders(procxind)+xmin
		ind=rmin_local
		if(rmin_local.lt.0) return
			
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			call RefPrtl_InnerBC_GPU<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,up_gpu,rmin_local,indi,indf)
		end do 
		do kc=1,Nchunk_test_prtl_gpu
			indi=(kc-1)*chunk_size_test_prtl_gpu+1
			indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
			call RefPrtl_InnerBC_GPU<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xtp_gpu,utp_gpu,rmin_local,indi,indf)
		end do 
		
		if(ind.ge.3.and.ind.le.mx-2) then
			!To ensure that the current is deposited on right place for the reflected particles
			f=1.0 ! needs to be changed
            !call SubstractFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jx_gpu,mx_gpu,my_gpu,mz_gpu,ind-1,ind)
			call AddFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jy_gpu,mx_gpu,my_gpu,mz_gpu,ind-1,ind,f)
			call AddFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jz_gpu,mx_gpu,my_gpu,mz_gpu,ind-1,ind,f)
			
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jx_gpu,mx_gpu,my_gpu,mz_gpu,ind-1)
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jy_gpu,mx_gpu,my_gpu,mz_gpu,ind-1)
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jz_gpu,mx_gpu,my_gpu,mz_gpu,ind-1)
		end if
	
	end subroutine RefInnerBC_Prtl

	subroutine RefOuterBC_Prtl(rmax) 
		real(psn) :: rmax,rmax_local,f
		integer   :: ind,n,kc
		integer   :: indi,indf
		
	    rmax_local=rmax-rborders(procxind)+xmin
		ind=rmax_local
		if(rmax_local.gt.mx) return
		
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			call RefPrtl_OuterBC_GPU<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,up_gpu,rmax_local,indi,indf)
		end do 
		do kc=1,Nchunk_test_prtl_gpu
			indi=(kc-1)*chunk_size_test_prtl_gpu+1
			indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
			call RefPrtl_OuterBC_GPU<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xtp_gpu,utp_gpu,rmax_local,indi,indf)
		end do 
		
		call CurrentBC_Outer
		
	end subroutine RefOuterBC_Prtl
	
	subroutine CurrentBC_Outer
		integer :: ind
		real(psn) :: f
		ind= BC_Rmax_Prtl-rborders(procxind)+xmin
		if(ind.ge.2.and.ind.le.mx-2) then
			!To ensure that the current is deposited on right place for the reflected particles	
			!f=(BC_Rmax_Prtl+0.5_psn)/(BC_Rmax_Prtl-0.5_psn)
			f=1.0 	
            !call SubstractFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jx_gpu,mx_gpu,my_gpu,mz_gpu,ind,ind-1)
			call AddFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jy_gpu,mx_gpu,my_gpu,mz_gpu,ind+1,ind,f)
			call AddFld_XLayer_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jz_gpu,mx_gpu,my_gpu,mz_gpu,ind+1,ind,f)
			
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jx_gpu,mx_gpu,my_gpu,mz_gpu,ind)
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jy_gpu,mx_gpu,my_gpu,mz_gpu,ind+1)
			call SetFldToZero_Xindex_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Jz_gpu,mx_gpu,my_gpu,mz_gpu,ind+1)
		end if
	end subroutine CurrentBC_Outer
	
	
	attributes(global) subroutine RefPrtl_InnerBC_GPU(x,u,wall,indi,indf)
	    real, dimension(:) :: x,u
		real, value        :: wall
		integer, value :: indi,indf
		integer        :: n
        n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
        if(n.gt.indf) return
		if(x(n).lt.wall) then 
			x(n)=wall+(wall-x(n))
			u(n)=-u(n)
		end if 
    end subroutine RefPrtl_InnerBC_GPU
	attributes(global) subroutine RefPrtl_OuterBC_GPU(x,u,wall,indi,indf)
	    real, dimension(:) :: x,u
		real, value        :: wall
		integer, value :: indi,indf
		integer        :: n
        n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
        if(n.gt.indf) return
		if(x(n).gt.wall) then 
			 x(n)=wall-(x(n)-wall)
			 u(n)=-u(n)
		end if 
    end subroutine RefPrtl_OuterBC_GPU	

	
	attributes(global) subroutine AddFld_XLayer_GPU(Fld,mx,my,mz,ind_src,ind_dest,f)
	    integer :: mx,my,mz
		real, dimension(mx,my,mz) :: Fld
	    real, value :: f
		integer, value :: ind_src, ind_dest 
		integer        :: j,k
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 
		if(j.le.my.and.k.le.mz) Fld(ind_dest,j,k)=Fld(ind_dest,j,k)+f*Fld(ind_src,j,k)
    end subroutine AddFld_Xlayer_GPU
	attributes(global) subroutine SubstractFld_XLayer_GPU(Fld,mx,my,mz,ind_src,ind_dest)
	    integer :: mx,my,mz
	    real, dimension(mx,my,mz) :: Fld
		integer, value :: ind_src, ind_dest 
		integer        :: j,k
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 
		if(j.le.my.and.k.le.mz) Fld(ind_dest,j,k)=Fld(ind_dest,j,k)-Fld(ind_src,j,k)
    end subroutine SubstractFld_Xlayer_GPU
end module cyl_bc_gpu