module movdep_gpu 
	use parameters 
	use vars 
	use var_gpu 
	use cudafor
contains 
	subroutine MoveDepositPrtlGPU
		integer :: kc, indi,indf
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		       call MoveDepositPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,qp_gpu,max_np_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
		                     Ex_gpu,Ey_gpu,Ez_gpu,Bx_gpu,By_gpu,Bz_gpu,Jx_gpu,Jy_gpu,Jz_gpu,Bx_ext0,By_ext0,Bz_ext0,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,indi,indf)
		end do 
	end subroutine MoveDepositPrtlGPU
	
	attributes(global) subroutine MoveDepositPrtlKernel(xp,yp,zp,up,vp,wp,flvp,qp,psize,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Bx_ext0,By_ext0,Bz_ext0,x1,x2,y1,y2,z1,z2,indi,indf)
		real, dimension(:) ::xp,yp,zp,up,vp,wp,qp
		integer, dimension(:) :: flvp
		real, dimension(:)  ::flvrqm
		integer, value :: psize 
		real, value :: c,sqc,cinv,qi
		integer :: x1,x2,y1,y2,z1,z2
		integer, value :: indi,indf
#ifndef twoD  		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z1+2)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#else 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,1:1)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#endif			
		real, value :: Bx_ext0,By_ext0,Bz_ext0
		
		integer :: n 
		integer :: stat
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g
		real :: pEx,pEy,pEz,pBx,pBy,pBz
		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: FldVec,wt 
#else 
        real,dimension(8) :: FldVec,wt
#endif 	
		!variables used in depositing current 
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
		
			
			
			
			
			
			  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
			  if(n.gt.indf) return 
			  
			  x0=xp(n)
	          y0=yp(n)
	          z0=zp(n)

			 ! load q/m for each particles into a vector
     	    if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
			  qm=0
		    else
  		      qm=flvrqm(flvp(n))*0.5_psn !for optimisation 0.5 is multiplied here itself
		    end if

			!interpolation
		  i=x0
		  j=y0
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=z0
          dz=z0-k
	      wt(1)=(1.0_psn-dx)*(1.0_psn-dy)*(1.0_psn-dz)
	      wt(2)=dx          *(1.0_psn-dy)*(1.0_psn-dz)
	      wt(3)=(1._psn-dx) *dy          *(1.0_psn-dz)
	      wt(4)=dx          *dy          *(1.0_psn-dz)
	      wt(5)=(1.0_psn-dx)*(1.0_psn-dy)*dz
	      wt(6)=dx          *(1.0_psn-dy)*dz
	      wt(7)=(1.0_psn-dx)*dy          *dz
	      wt(8)=dx          *dy          *dz
#else
          k=1
		  wt(1)=(1.0_psn-dx)*(1.0_psn-dy)
		  wt(2)=dx          *(1.0_psn-dy)
		  wt(3)=(1._psn-dx) *dy
		  wt(4)=dx          *dy
#endif
		  FldVec(1)=Ex(i-1,j,k)+Ex(i,j,k)
		  FldVec(2)=Ex(i,j,k)+Ex(i+1,j,k)
		  FldVec(3)=Ex(i-1,j+1,k)+Ex(i,j+1,k)
		  FldVec(4)=Ex(i,j+1,k)+Ex(i+1,j+1,k)
#ifndef twoD
		  FldVec(5)=Ex(i-1,j,k+1)+Ex(i,j,k+1)
		  FldVec(6)=Ex(i,j,k+1)+Ex(i+1,j,k+1)
		  FldVec(7)=Ex(i-1,j+1,k+1)+Ex(i,j+1,k+1)
		  FldVec(8)=Ex(i,j+1,k+1)+Ex(i+1,j+1,k+1)
#endif

#ifdef twoD
          pEx=0.5_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#else
          pEx=0.5_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#endif

		  FldVec(1)=Ey(i,j-1,k)+Ey(i,j,k)
		  FldVec(2)=Ey(i+1,j-1,k)+Ey(i+1,j,k)
		  FldVec(3)=Ey(i,j,k)+Ey(i,j+1,k)
		  FldVec(4)=Ey(i+1,j,k)+Ey(i+1,j+1,k)
#ifndef twoD
		  FldVec(5)=Ey(i,j-1,k+1)+Ey(i,j,k+1)
		  FldVec(6)=Ey(i+1,j-1,k+1)+Ey(i+1,j,k+1)
		  FldVec(7)=Ey(i,j,k+1)+Ey(i,j+1,k+1)
		  FldVec(8)=Ey(i+1,j,k+1)+Ey(i+1,j+1,k+1)
#endif

#ifdef twoD
          pEy=0.5_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#else
          pEy=0.5_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#endif

#ifndef twoD
		  FldVec(1)=Ez(i,j,k-1)+Ez(i,j,k)
		  FldVec(2)=Ez(i+1,j,k-1)+Ez(i+1,j,k)
		  FldVec(3)=Ez(i,j+1,k-1)+Ez(i,j+1,k)
		  FldVec(4)=Ez(i+1,j+1,k-1)+Ez(i+1,j+1,k)
		  FldVec(5)=Ez(i,j,k)+Ez(i,j,k+1)
		  FldVec(6)=Ez(i+1,j,k)+Ez(i+1,j,k+1)
		  FldVec(7)=Ez(i,j+1,k)+Ez(i,j+1,k+1)
		  FldVec(8)=Ez(i+1,j+1,k)+Ez(i+1,j+1,k+1)
          pEz=0.5_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#else
			FldVec(1)=Ez(i,j,k)
			FldVec(2)=Ez(i+1,j,k)
			FldVec(3)=Ez(i,j+1,k)
			FldVec(4)=Ez(i+1,j+1,k)
            pEz=(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#endif


#ifndef twoD
		  FldVec(1)=Bx(i,j-1,k-1)+Bx(i,j,k-1)+Bx(i,j-1,k)+Bx(i,j,k)
		  FldVec(2)=Bx(i+1,j-1,k-1)+Bx(i+1,j,k-1)+Bx(i+1,j-1,k)+Bx(i+1,j,k)
		  FldVec(3)=Bx(i,j,k-1)+Bx(i,j+1,k-1)+Bx(i,j,k)+Bx(i,j+1,k)
		  FldVec(4)=Bx(i+1,j,k-1)+Bx(i+1,j+1,k-1)+Bx(i+1,j,k)+Bx(i+1,j+1,k)
		  FldVec(5)=Bx(i,j-1,k)+Bx(i,j,k)+Bx(i,j-1,k+1)+Bx(i,j,k+1)
		  FldVec(6)=Bx(i+1,j-1,k)+Bx(i+1,j,k)+Bx(i+1,j-1,k+1)+Bx(i+1,j,k+1)
		  FldVec(7)=Bx(i,j,k)+Bx(i,j+1,k)+Bx(i,j,k+1)+Bx(i,j+1,k+1)
		  FldVec(8)=Bx(i+1,j,k)+Bx(i+1,j+1,k)+Bx(i+1,j,k+1)+Bx(i+1,j+1,k+1)
          pBx=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))+Bx_ext0
#else
			FldVec(1)=Bx(i,j-1,k)+Bx(i,j,k)
			FldVec(2)=Bx(i+1,j-1,k)+Bx(i+1,j,k)
			FldVec(3)=Bx(i,j,k)+Bx(i,j+1,k)
			FldVec(4)=Bx(i+1,j,k)+Bx(i+1,j+1,k)
	        pBx=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))+Bx_ext0
#endif

#ifndef twoD
		  FldVec(1)=By(i-1,j,k-1)+By(i-1,j,k)+By(i,j,k-1)+By(i,j,k)
		  FldVec(2)=By(i,j,k-1)+By(i,j,k)+By(i+1,j,k-1)+By(i+1,j,k)
		  FldVec(3)=By(i-1,j+1,k-1)+By(i-1,j+1,k)+By(i,j+1,k-1)+By(i,j+1,k)
		  FldVec(4)=By(i,j+1,k-1)+By(i,j+1,k)+By(i+1,j+1,k-1)+By(i+1,j+1,k)
		  FldVec(5)=By(i-1,j,k)+By(i-1,j,k+1)+By(i,j,k)+By(i,j,k+1)
		  FldVec(6)=By(i,j,k)+By(i,j,k+1)+By(i+1,j,k)+By(i+1,j,k+1)
		  FldVec(7)=By(i-1,j+1,k)+By(i-1,j+1,k+1)+By(i,j+1,k)+By(i,j+1,k+1)
		  FldVec(8)=By(i,j+1,k)+By(i,j+1,k+1)+By(i+1,j+1,k)+By(i+1,j+1,k+1)
          pBy=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))+By_ext0
#else
			FldVec(1)=By(i-1,j,k)+By(i,j,k)
			FldVec(2)=By(i,j,k)+By(i+1,j,k)
			FldVec(3)=By(i-1,j+1,k)+By(i,j+1,k)
			FldVec(4)=By(i,j+1,k)+By(i+1,j+1,k)
	        pBz=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))+By_ext0
#endif

		  FldVec(1)=Bz(i-1,j-1,k)+Bz(i-1,j,k)+Bz(i,j-1,k)+Bz(i,j,k)
		  FldVec(2)=Bz(i,j-1,k)+Bz(i,j,k)+Bz(i+1,j-1,k)+Bz(i+1,j,k)
		  FldVec(3)=Bz(i-1,j,k)+Bz(i-1,j+1,k)+Bz(i,j,k)+Bz(i,j+1,k)
		  FldVec(4)=Bz(i,j,k)+Bz(i,j+1,k)+Bz(i+1,j,k)+Bz(i+1,j+1,k)
#ifndef twoD
		  FldVec(5)=Bz(i-1,j-1,k+1)+Bz(i-1,j,k+1)+Bz(i,j-1,k+1)+Bz(i,j,k+1)
		  FldVec(6)=Bz(i,j-1,k+1)+Bz(i,j,k+1)+Bz(i+1,j-1,k+1)+Bz(i+1,j,k+1)
		  FldVec(7)=Bz(i-1,j,k+1)+Bz(i-1,j+1,k+1)+Bz(i,j,k+1)+Bz(i,j+1,k+1)
		  FldVec(8)=Bz(i,j,k+1)+Bz(i,j+1,k+1)+Bz(i+1,j,k+1)+Bz(i+1,j+1,k+1)
#endif

#ifdef twoD
        pBz=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))+Bz_ext0
#else
        pBz=0.25_psn*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))+Bz_ext0
#endif




            pEx=pEx*qm
            pEy=pEy*qm
            pEz=pEz*qm
            pBx=pBx*qm
            pBy=pBy*qm
            pBz=pBz*qm

			   !Boris Pusher
               u0=c*up(n)+pEx
               v0=c*vp(n)+pEy
               w0=c*wp(n)+pEz

               g=1.0_psn/sqrt(sqc+u0**2+v0**2+w0**2)   ! 1/c*gamma
               pBx=g*pBx
               pBy=g*pBy
               pBz=g*pBz

               f=2.0_psn/(1.0_psn+pBx*pBx+pBy*pBy+pBz*pBz)
               u1=(u0+v0*pBz-w0*pBy)*f
               v1=(v0+w0*pBx-u0*pBz)*f
               w1=(w0+u0*pBy-v0*pBx)*f

               u0=u0+v1*pBz-w1*pBy+pEx
               v0=v0+w1*pBx-u1*pBz+pEy
               w0=w0+u1*pBy-v1*pBx+pEz

               up(n)=u0*cinv
               vp(n)=v0*cinv
               wp(n)=w0*cinv

               g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g
               yp(n)=yp(n) + vp(n)*g
               zp(n)=zp(n) + wp(n)*g

!  		do n=1,VecBlockSizeThis
!  	        call DepositCurrentPIC(x0(n),y0(n),z0(n),xp(n+off),yp(n+off),zp(n+off),qp(n+off)) !deposit charges
!  		end do

          qthis=qp(n)*qi   ! q = particle's weight X sign of the charge

          i1=floor(x0)
          i2=floor(xp(n))
          j1=floor(y0)
          j2=floor(yp(n))
          k1=floor(z0)
          k2=floor(zp(n))
#ifdef twoD
        k1=1
        k2=1
#endif
               xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),0.5_psn*(x0+xp(n))))
               yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),0.5_psn*(y0+yp(n))))
               zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),0.5_psn*(z0+zp(n))))
          Fx1=qthis*(xr-x0)
          Fy1=qthis*(yr-y0)
          Fz1=qthis*(zr-z0)

          Wx1=0.5_psn*(x0+xr)-i1
          Wy1=0.5_psn*(y0+yr)-j1

          Wx2=0.5_psn*(xp(n)+xr)-i2
          Wy2=0.5_psn*(yp(n)+yr)-j2

#ifdef twoD
          Wz1=0.0_psn
          Wz2=0.0_psn
#else
          Wz1=0.5_psn*(z0+zr)-k1
          Wz2=0.5_psn*(zp(n)+zr)-k2
#endif


          Fx2=qthis*(xp(n)-xr)
          Fy2=qthis*(yp(n)-yr)
          Fz2=qthis*(zp(n)-zr)


          Jx1(1)=Fx1 * (1.0_psn-Wy1)*(1.0_psn-Wz1)
          Jx1(2)=Fx1 *  Wy1    *(1.0_psn-Wz1)
#ifndef twoD
              Jx1(3)= Fx1 * (1.0_psn-Wy1) * Wz1
              Jx1(4)= Fx1 *  Wy1    * Wz1
#endif

          Jx2(1)=Fx2 * (1.0_psn-Wy2)*(1.0_psn-Wz2)
          Jx2(2)=Fx2 *  Wy2    *(1.0_psn-Wz2)
#ifndef twoD
               Jx2(3)=Fx2 * (1.0_psn-Wy2)* Wz2
               Jx2(4)=Fx2 *  Wy2    * Wz2
#endif


          Jy1(1)=Fy1 * (1.0_psn-Wx1)*(1.0_psn-Wz1)
          Jy1(2)=Fy1 *  Wx1    *(1.0_psn-Wz1)
#ifndef twoD
               Jy1(3)=Fy1 * (1.0_psn-Wx1)* Wz1
               Jy1(4)=Fy1 *  Wx1    * Wz1
#endif
          Jy2(1)=Fy2 * (1.0_psn-Wx2)*(1.0_psn-Wz2)
          Jy2(2)=Fy2 *  Wx2    *(1.0_psn-Wz2)
#ifndef twoD
             Jy1(3)= Fy2 * (1.0_psn-Wx2)* Wz2
             Jy1(4)= Fy2 *  Wx2    * Wz2
#endif


          Jz1(1)=Fz1 * (1.0_psn-Wx1)*(1.0_psn-Wy1)
          Jz1(2)=Fz1 *  Wx1    *(1.0_psn-Wy1)
          Jz1(3)=Fz1 * (1.0_psn-Wx1)* Wy1
          Jz1(4)=Fz1 *  Wx1    * Wy1

          Jz2(1)=Fz2 * (1.0_psn-Wx2)*(1.0_psn-Wy2)
          Jz2(2)=Fz2 *  Wx2    *(1.0_psn-Wy2)
          Jz2(3)=Fz2 * (1.0_psn-Wx2)* Wy2
          Jz2(4)=Fz2 *  Wx2    * Wy2








          stat=atomicAdd(Jx(i1,j1,  k1  ), Jx1(1) )
          stat=atomicAdd(Jx(i1,j1+1,k1  ), Jx1(2) )
#ifndef twoD
               stat=atomicAdd(Jx(i1,j1,  k1+1), Jx1(3))
               stat=atomicAdd(Jx(i1,j1+1,k1+1), Jx1(4))
#endif

          stat=atomicAdd(Jx(i2,j2,  k2  ), Jx2(1) )
          stat=atomicAdd(Jx(i2,j2+1,k2  ), Jx2(1) )
#ifndef twoD
               stat= atomicAdd(Jx(i2,j2,  k2+1), Jx2(3))
               stat= atomicAdd(Jx(i2,j2+1,k2+1), Jx2(4))
#endif


          stat=atomicAdd(Jy(i1  ,j1,k1  ), Jy1(1))
          stat=atomicAdd(Jy(i1+1,j1,k1  ), Jy1(2))
#ifndef twoD
               stat=atomicAdd(Jy(i1  ,j1,k1+1), Jy1(3))
               stat=atomicAdd(Jy(i1+1,j1,k1+1), Jy1(4))
#endif
          stat=atomicAdd(Jy(i2  ,j2,k2  ), Jy2(1))
          stat=atomicAdd( Jy(i2+1,j2,k2  ), Jy2(2))
#ifndef twoD
               stat=atomicAdd(Jy(i2  ,j2,k2+1), Jy2(3))
               stat=atomicAdd(Jy(i2+1,j2,k2+1), Jy2(4))
#endif


          stat=atomicAdd(Jz(i1  ,j1  ,k1), Jz1(1))
          stat=atomicAdd(Jz(i1+1,j1  ,k1), Jz1(2))
          stat=atomicAdd(Jz(i1  ,j1+1,k1), Jz1(3))
          stat=atomicAdd(Jz(i1+1,j1+1,k1), Jz1(4))

          stat=atomicAdd(Jz(i2  ,j2  ,k2), Jz2(1))
          stat=atomicAdd(Jz(i2+1,j2  ,k2), Jz2(2))
          stat=atomicAdd(Jz(i2  ,j2+1,k2), Jz2(3))
          stat=atomicAdd(Jz(i2+1,j2+1,k2), Jz2(4))

!           if(xp(n).ge.1027.0_psn) then
!             xp(n)=xp(n)-1024.0_psn
!           else if(xp(n).lt.3.0_psn) then
!             xp(n)=1024.0_psn+xp(n)
!           end if
!
!
!           if(yp(n).ge.1027.0_psn) then
!             yp(n)=yp(n)-1024.0_psn
!           else if(yp(n).lt.3.0_psn) then
!             yp(n)=1024.0_psn+yp(n)
!           end if
              
! Periodic Boundary conditions 
#ifdef twoD
             if(zp(n).ge.2.0_psn) then
               zp(n)=zp(n)-1.0_psn
             else if(zp(n).lt.1.0_psn) then
               zp(n)=1.0_psn+zp(n)
             end if 
#endif       

end subroutine MoveDepositPrtlKernel
	
end module movdep_gpu