module deposit_gpu 
     use parameters
     use vars
	 use var_gpu
     implicit none 
contains 
	
	subroutine ResetJwide
		call ResetJwideKernel<<<grid,tBlock>>>(Jx_wide_gpu,mx,my,mz,Jwidth_gpu)
		call ResetJwideKernel<<<grid,tBlock>>>(Jy_wide_gpu,mx,my,mz,Jwidth_gpu)
		call ResetJwideKernel<<<grid,tBlock>>>(Jz_wide_gpu,mx,my,mz,Jwidth_gpu)
	end subroutine ResetJwide
	
	subroutine ReduceJwide
		call MergeJwideKernel<<<grid,tBlock>>>(Jx_gpu,Jx_wide_gpu,mx,my,mz,Jwidth_gpu)
		call MergeJwideKernel<<<grid,tBlock>>>(Jy_gpu,Jy_wide_gpu,mx,my,mz,Jwidth_gpu)
		call MergeJwideKernel<<<grid,tBlock>>>(Jz_gpu,Jz_wide_gpu,mx,my,mz,Jwidth_gpu)
	end subroutine ReduceJwide
	
	attributes(global) subroutine ResetJwideKernel(J_wide,mx,my,mz,n)
		integer, value :: mx,my,mz,n
		real, dimension(mx,my,mz,n) :: J_wide
		integer :: i,j,k,ni
		
	    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
	    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
	    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#ifdef twoD
        k=1
#endif		
		
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) then 
			do ni=1,n
				J_wide(i,j,k,ni)=0.0
			end do 
		end if 
	
    end subroutine ResetJwideKernel
	
	attributes(global) subroutine MergeJwideKernel(J_main,J_wide,mx,my,mz,n)
		integer, value :: mx,my,mz,n
		real, dimension(mx,my,mz) :: J_main
		real, dimension(mx,my,mz,n) :: J_wide
		integer :: i,j,k,ni
		
	    i = (blockIdx%x-1)*blockDim%x + threadIdx%x
	    j = (blockIdx%y-1)*blockDim%y + threadIdx%y
	    k = (blockIdx%z-1)*blockDim%z + threadIdx%z
#ifdef twoD
        k=1
#endif	

		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) then
			
			do ni=1,n
				J_main(i,j,k)=J_main(i,j,k)+J_wide(i,j,k,ni)
			end do
			
		end if 
	
    end subroutine MergeJwideKernel
	
	 attributes(device) subroutine DepositCurrentPIC_GPU(x0,y0,z0,x,y,z,q,n,qi, Jx,Jy,Jz,mx,my,mz,jwidth)
		 real :: x0,y0,z0,x,y,z,q,qi
		 integer :: n,mx,my,mz,jwidth
		 real, dimension(mx,my,mz,jwidth) :: Jx,Jy,Jz
 		 integer :: i1,j1,k1,i2,j2,k2
         real ::Fx1, Fx2, Fy1, Fy2, Fz1, Fz2
         real ::Wx1, Wx2, Wy1, Wy2, Wz1, Wz2
         real ::xr,yr,zr
         real ::qthis 
#ifdef twoD
         real,dimension(2) :: Jx1,Jy1,Jx2,Jy2
#else 
         real,dimension(4) :: Jx1,Jy1,Jx2,Jy2
#endif 	
         real,dimension(4) :: Jz1,Jz2
 		 integer :: j_ind
		 integer :: stat
		 

          qthis=q*qi   ! q = particle's weight X sign of the charge

          i1=floor(x0)
          i2=floor(x)
          j1=floor(y0)
          j2=floor(y)
          k1=floor(z0)
          k2=floor(z)

          xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),0.5*(x0+x)))
          yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),0.5*(y0+y)))
          zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),0.5*(z0+z)))
			   

          Fx1=qthis*(xr-x0)
          Fy1=qthis*(yr-y0)
          Fz1=qthis*(zr-z0)

          Wx1=0.5*(x0+xr)-i1
          Wy1=0.5*(y0+yr)-j1

          Wx2=0.5*(x+xr)-i2
          Wy2=0.5*(y+yr)-j2

#ifdef twoD
          Wz1=0.0
          Wz2=0.0
#else
          Wz1=0.5*(z0+zr)-k1
          Wz2=0.5*(z+zr)-k2
#endif


          Fx2=qthis*(x-xr)
          Fy2=qthis*(y-yr)
          Fz2=qthis*(z-zr)

#ifdef twoD
          k1=1
          k2=1
#endif


		  j_ind=mod(n-1,jwidth)+1


          Jx1(1)=Fx1 * (1.0-Wy1)*(1.0-Wz1)
          Jx1(2)=Fx1 *  Wy1    *(1.0-Wz1)
#ifndef twoD
          Jx1(3)= Fx1 * (1.0-Wy1) * Wz1
          Jx1(4)= Fx1 *  Wy1    * Wz1
#endif

          stat=atomicAdd(Jx(i1,j1,  k1 ,j_ind  ), Jx1(1) )
          stat=atomicAdd(Jx(i1,j1+1,k1 ,j_ind ), Jx1(2) )
#ifndef twoD
          stat=atomicAdd(Jx(i1,j1,  k1+1,j_ind), Jx1(3))
          stat=atomicAdd(Jx(i1,j1+1,k1+1,j_ind), Jx1(4))
#endif

          Jx2(1)=Fx2 * (1.0-Wy2)*(1.0-Wz2)
          Jx2(2)=Fx2 *  Wy2    *(1.0-Wz2)
#ifndef twoD
          Jx2(3)=Fx2 * (1.0-Wy2)* Wz2
          Jx2(4)=Fx2 *  Wy2    * Wz2
#endif

          stat=atomicAdd(Jx(i2,j2,  k2 ,j_ind ), Jx2(1) )
          stat=atomicAdd(Jx(i2,j2+1,k2 ,j_ind ), Jx2(2) )
#ifndef twoD
          stat= atomicAdd(Jx(i2,j2,  k2+1,j_ind), Jx2(3))
          stat= atomicAdd(Jx(i2,j2+1,k2+1,j_ind), Jx2(4))
#endif


          Jy1(1)=Fy1 * (1.0-Wx1)*(1.0-Wz1)
          Jy1(2)=Fy1 *  Wx1    *(1.0-Wz1)
#ifndef twoD
          Jy1(3)=Fy1 * (1.0-Wx1)* Wz1
          Jy1(4)=Fy1 *  Wx1    * Wz1
#endif

          stat=atomicAdd(Jy(i1  ,j1,k1  ,j_ind), Jy1(1))
          stat=atomicAdd(Jy(i1+1,j1,k1  ,j_ind), Jy1(2))
#ifndef twoD
          stat=atomicAdd(Jy(i1  ,j1,k1+1,j_ind), Jy1(3))
          stat=atomicAdd(Jy(i1+1,j1,k1+1,j_ind), Jy1(4))
#endif


          Jy2(1)=Fy2 * (1.0-Wx2)*(1.0-Wz2)
          Jy2(2)=Fy2 *  Wx2    *(1.0-Wz2)
#ifndef twoD
          Jy2(3)= Fy2 * (1.0-Wx2)* Wz2
          Jy2(4)= Fy2 *  Wx2    * Wz2
#endif

          stat=atomicAdd(Jy(i2  ,j2,k2  ,j_ind), Jy2(1))
          stat=atomicAdd( Jy(i2+1,j2,k2 ,j_ind), Jy2(2))
#ifndef twoD
          stat=atomicAdd(Jy(i2  ,j2,k2+1,j_ind), Jy2(3))
          stat=atomicAdd(Jy(i2+1,j2,k2+1,j_ind), Jy2(4))
#endif


          Jz1(1)=Fz1 * (1.0-Wx1)*(1.0-Wy1)
          Jz1(2)=Fz1 *  Wx1    *(1.0-Wy1)
          Jz1(3)=Fz1 * (1.0-Wx1)* Wy1
          Jz1(4)=Fz1 *  Wx1    * Wy1
		  

          stat=atomicAdd(Jz(i1  ,j1  ,k1,j_ind), Jz1(1))
          stat=atomicAdd(Jz(i1+1,j1  ,k1,j_ind), Jz1(2))
          stat=atomicAdd(Jz(i1  ,j1+1,k1,j_ind), Jz1(3))
          stat=atomicAdd(Jz(i1+1,j1+1,k1,j_ind), Jz1(4))		  

          Jz2(1)=Fz2 * (1.0-Wx2)*(1.0-Wy2)
          Jz2(2)=Fz2 *  Wx2    *(1.0-Wy2)
          Jz2(3)=Fz2 * (1.0-Wx2)* Wy2
          Jz2(4)=Fz2 *  Wx2    * Wy2


          stat=atomicAdd(Jz(i2  ,j2  ,k2,j_ind), Jz2(1))
          stat=atomicAdd(Jz(i2+1,j2  ,k2,j_ind), Jz2(2))
          stat=atomicAdd(Jz(i2  ,j2+1,k2,j_ind), Jz2(3))
          stat=atomicAdd(Jz(i2+1,j2+1,k2,j_ind), Jz2(4))

		 
		 
	 end subroutine 

end module deposit_gpu