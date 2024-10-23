module em_update_gpu
	use parameters
	use vars 
	use cyl_vars
	use var_gpu
	implicit none 
contains 
     subroutine UpdateEfieldGPU
		 integer  :: ioff
		 real(psn):: rmin
		 ioff=3
		 if(procxind.eq.0) ioff=4
		 rmin=rborders(procxind)
		 
         call UpdateEfieldGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,mx,my,mz,fldc,dtheta,ioff,rmin) 
		 
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) then 
			 call UpdateEFldAxis_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Ex_gpu,Ey_gpu,Ez_gpu,mx,my,mz,ny)
		 end if
		 
	end subroutine UpdateEfieldGPU
	
	attributes(global) subroutine UpdateEfieldGPUKernel(Ex,Ey,Ez,Bx,By,Bz,mx,my,mz,fldc, dtheta, ioff, rmin)	
	    integer, value :: mx,my,mz   
		real, dimension(mx,my,mz) :: Ex,Ey,Ez,Bx,By,Bz		
		real, value    :: fldc, dtheta, rmin
		integer, value :: ioff 
		integer :: i,j,k, kmax
		real    :: r,rp_half,rm_half		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x+2 
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y+2 
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z+2
		kmax=mz-3
#else
        k=1
		kmax=1
#endif  

			    if((i.ge.ioff).and.(i.le.mx-3)) then 
				    if(j.le.my-3) then					 
					    if(k.le.kmax) then 				  
				  
							r=(i-3.0)+rmin
							rp_half=r+0.5
							rm_half=r-0.5
#ifndef twoD
		                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))/(rp_half*dtheta) - fldc*(By(i,j,k)-By(i,j,k-1))
		                    Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
		                    Ez(i,j,k)=Ez(i,j,k)+fldc*(rp_half*By(i,j,k)-rm_half*By(i-1,j,k))/r - fldc*(Bx(i,j,k)-Bx(i,j-1,k))/(r*dtheta)
#else
		                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))/(rp_half*dtheta)
		                    Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
		                    Ez(i,j,k)=Ez(i,j,k)+fldc*( (rp_half*By(i,j,k)-rm_half*By(i-1,j,k)) - (Bx(i,j,k)-Bx(i,j-1,k))/dtheta  )/r
#endif
                       end if 
                   end if 
               end if
     end subroutine UpdateEfieldGPUKernel
	 
     subroutine UpdateBfieldHalfGPU
		 integer  :: ioff
		 real(psn):: rmin
		 integer  :: j,k 
		 
		 ioff=1
		 if(procxind.eq.0) ioff=4
		 rmin=rborders(procxind)
		
         call UpdateBfieldHalfGPUKernel<<<grid,tBlock>>>(Bx_gpu,By_gpu,Bz_gpu,Ex_gpu,Ey_gpu,Ez_gpu,mx,my,mz,fld_halfc,dtheta,ioff,rmin) 
	     
! 		 if(inc_axis) then
! 		    if(procxind.eq.0) then
! 				call CopyEyAxisGPU<<<grid_ax,tBlock_ax>>>(Ey_gpu,Ey_ax_gpu,mx,my,mz)
! 				Ey_ax_host=Ey_ax_gpu
! 			    call UpdateBzAxis_Host
! 				Bz_ax_gpu=Bz_ax ! had only db/dt
! 				call UpdateBzAxisGPU<<<grid_ax,tBlock_ax>>>(Bz_gpu,Bz_ax_gpu,mx,my,mz)  ! B_updated = B+dB/dt
! 				call CopyBzAxisGPU<<<grid_ax,tBlock_ax>>>(Bz_gpu,Bz_ax_gpu,mx,my,mz)
! 		     end if
! 	      end if 

			if(inc_axis) then
			  if(procxind.eq.0) then
				  
				  Bz=Bz_gpu
				  Ey=Ey_gpu
				  call UpdateBzAxis

#ifdef twoD
                  do k=1,1
#else					  				  
				  do k=1,mz-1
#endif					  
					  do j=1,my-1
						   Bz(3,j,k)=Bz_ax(k)
					  end do
				  end do

			      Bz_gpu=Bz
				  
				  call UpdateByAxis_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(By_gpu,Ex_gpu,Ez_gpu,mx,my,mz,fld_halfc)
				  
				  call UpdateBFldAxis_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Bx_gpu,By_gpu,Bz_gpu,mx,my,mz,ny)
				  
			  end if
			end if

	end subroutine UpdateBfieldHalfGPU
	
	 subroutine UpdateBzAxis
		 integer :: j,k
		 
		 Bz_ax=0
#ifdef twoD
         do k=1,1 
#else 		 
		 do k=1,mz-1
#endif			 		 
			 do j=3,my-3 !First del Bz
				 Bz_ax(k)=Bz_ax(k)-fld_halfc*ax_perm_area*Ey(4,j,k) 
			 end do 
		     Bz_ax(k)=Bz(3,3,k)+Bz_ax(k)
		 end do 	
	 end subroutine UpdateBzAxis

	
	attributes(global) subroutine UpdateBfieldHalfGPUKernel(Bx,By,Bz,Ex,Ey,Ez,mx,my,mz,fld_halfc, dtheta, ioff, rmin)	
	    integer, value :: mx,my,mz		
		real, dimension(mx,my,mz) :: Bx,By,Bz,Ex,Ey,Ez 
		real, value    :: fld_halfc, dtheta, rmin
		integer, value :: ioff
		integer :: i,j,k, kmax
		real    :: r,rp_half,rp		 		    

		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
		kmax=mz-1
#else
        k=1
		kmax=1
#endif

			    if((i.ge.ioff).and.(i.le.mx-1)) then 
				    if(j.le.my-1) then					 
					    if(k.le.kmax) then 				  				  

						 r=(i-3.0)+rmin
						 rp_half=r+0.5
						 rp=r+1.0
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))/(r*dtheta) + fld_halfc*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*( (rp*Ey(i+1,j,k)-r*Ey(i,j,k)) - (Ex(i,j+1,k)-Ex(i,j,k))/dtheta )/rp_half
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))/(r*dtheta)
                         By(i,j,k)=By(i,j,k)+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
						 Bz(i,j,k)=Bz(i,j,k)-fld_halfc*( (rp*Ey(i+1,j,k)-r*Ey(i,j,k)) - (Ex(i,j+1,k)-Ex(i,j,k))/dtheta )/rp_half
#endif
                       end if 
                   end if 
               end if

     end subroutine UpdateBfieldHalfGPUKernel
	 
	 
	subroutine AddCurrentGPU
		 integer  :: ioff
		 real(psn):: rmin
		 
		 ioff=3
		 if(procxind.eq.0) ioff=4
		 rmin=rborders(procxind)
		 call AddCurrentGPUKernel<<<grid,tBlock>>>(Ex_gpu,Ey_gpu,Ez_gpu,Jx_gpu,Jy_gpu,Jz_gpu, mx,my,mz, dtheta, ioff, rmin) 
		 
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) then 
			 call UpdateEFldAxis_GPU<<<tGrid_gpu_YZedge1,tBlock_gpu_YZedge1>>>(Ex_gpu,Ey_gpu,Ez_gpu,mx,my,mz,ny)
		 end if
	end subroutine AddCurrentGPU
	

	attributes(global) subroutine AddCurrentGPUKernel(Ex,Ey,Ez,Jx,Jy,Jz, mx,my,mz, dtheta, ioff, rmin)
	    integer, value :: mx,my,mz    
		real, dimension(mx,my,mz) :: Ex,Ey,Ez,Jx,Jy,Jz		
        real, value    :: dtheta, rmin
		real    :: r,rp_half	
		integer, value :: ioff
		integer :: i,j,k,kmax
		
		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
        kmax=mz-3
#else
        k=1
		kmax=1
#endif  
        if((i.ge.ioff).and.(i.le.mx-3).and.(j.le.my-3).and.(k.le.kmax)) then 
	        r=(i-3.0)+rmin
		    rp_half=r+0.5
            Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)/(rp_half*dtheta)
            Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)
            Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)/(r*dtheta)
		end if 	
	end subroutine AddCurrentGPUKernel 
	
	
    attributes(global) subroutine UpdateEFldAxis_GPU(Ex,Ey,Ez,mx,my,mz,ny)
	    integer, value :: mx,my,mz
		integer, value :: ny
		integer :: j,k,jp,jm,jpp 
	    real, dimension(mx,my,mz) :: Ex,Ey,Ez
		
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
! 			Ex(3,j,k) = 0.5_psn*(Ey(4,jm,k) - Ey(4,jp,k))
! 			Ex(2,j,k) = -Ex(4,jpp,k)
!
! 			Ey(3,j,k) = -Ey(4,jpp,k)
! 			Ez(3,j,k) = Ez(4,jpp,k)

			Ex(3,j,k)=0.0
			Ex(2,j,k)=-Ex(4,j,k)

			Ey(3,j,k)=-Ey(4,j,k)
			Ez(3,j,k)=-Ez(4,j,k)
			
	    end if 
    end subroutine UpdateEFldAxis_GPU 
	
    attributes(global) subroutine UpdateBFldAxis_GPU(Bx,By,Bz,mx,my,mz,ny)
	    integer, value :: mx,my,mz
		integer, value :: ny
	    real, dimension(mx,my,mz) :: Bx,By,Bz
		integer        :: j,k, jpp
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 	
        if(j.le.my.and.k.le.mz) then 
			


! 			jpp = j + ny/2
! 			if(jpp.gt.my-3) jpp= jpp - ny
!
! 			Bx(3,j,k) = -Bx(4,jpp,k)
! 			By(2,j,k) = -By(4,jpp,k)
! 			Bz(2,j,k) = Bz(4,jpp,k)
			
			! setting the fld at the axis (as seen by the particles) to zero (a special case of extrapolation)
		    Bx(3,j,k)=-Bx(4,j,k)
		    By(3,j,k)=0.0_psn
		    By(2,j,k)=-By(4,j,k)
		    !Bz(2,j,k)=2*Bz(3,j,k)-Bz(4,j,k)
		    Bz(2,j,k)=Bz(3,j,k)
								
								
	    end if 
    end subroutine UpdateBFldAxis_GPU 
	
    attributes(global) subroutine UpdateByAxis_GPU(By,Ex,Ez,mx,my,mz,fld_halfc)
	    integer, value :: mx,my,mz
	    real, dimension(mx,my,mz) :: Ex,By,Ez
		real, value    :: fld_halfc
		integer        :: j,k
		
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 	
        if(j.le.my.and.k.le.mz) then 
			
                By(3,j,k)=By(3,j,k)-fld_halfc*(Ex(3,j,k+1)-Ex(3,j,k))+fld_halfc*(Ez(4,j,k)-Ez(3,j,k))
								
	    end if 
    end subroutine UpdateByAxis_GPU 
	
!--------------------------------------------------------------
! Update B_z at the axis
!--------------------------------------------------------------	
	
	
	 subroutine UpdateBzAxis_Host
		 integer :: j,k
		 Bz_ax=0
#ifdef twoD
         do k=1,1
#else
		 do k=1,mz-1
#endif
			 do j=3,my-3 !First del Bz
				 Bz_ax(k)=Bz_ax(k)-fld_halfc*ax_perm_area*Ey_ax_host(j,k)
			 end do
		 end do
	 end subroutine UpdateBzAxis_host

	attributes(global) subroutine CopyEyAxisGPU(Ey,Ey_ax,mx,my,mz)
	    integer, value :: mx,my,mz
		real, dimension(mx,my,mz) :: Ey
		real, dimension(my,mz)    :: Ey_ax
		integer :: j,k

		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#else
        k=1
#endif
        if((j.le.my).and.(k.le.mz)) then
			Ey_ax(j,k)=Ey(4,j,k)
		end if

	end subroutine CopyEyAxisGPU
	attributes(global) subroutine CopyBzAxisGPU(Bz,Bz_ax,mx,my,mz)
	    integer, value :: mx,my,mz
		real, dimension(mx,my,mz) :: Bz
		real, dimension(mz)    :: Bz_ax
		integer :: j,k

		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#else
        k=1
#endif
        if((j.le.my).and.(k.le.mz)) then
			Bz(3,j,k)=Bz_ax(k)
		end if

	end subroutine CopyBzAxisGPU
	attributes(global) subroutine UpdateBzAxisGPU(Bz,Bz_ax,mx,my,mz)
	    integer, value :: mx,my,mz
		real, dimension(mx,my,mz) :: Bz
		real, dimension(mz)    :: Bz_ax
		integer :: j,k

		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD
		k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#else
        k=1
#endif
        if((j.eq.1).and.(k.le.mz)) then
			Bz_ax(k)= Bz_ax(k) + Bz(3,3,k)
		end if

	end subroutine UpdateBzAxisGPU
	
	
end module em_update_gpu