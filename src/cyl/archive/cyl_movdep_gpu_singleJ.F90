module movdep_gpu 
	use parameters 
	use vars 
	use var_gpu 
	use fields_gpu !Fitlering subroutines are in fields_gpu
	use interp_gpu
	use cyl_bc
	use cudafor
	use cudadevice
	real, texture, pointer, dimension(:,:,:) :: tEx,tEy,tEz,tBx,tBy,tBz
contains 
	subroutine MoveDepositPrtlGPU
		integer :: kc, indi,indf
			
		call UpdateFldAxis
		 
		call SetTexFlds
		
		tEx=>TexEx_gpu
		tEy=>TexEy_gpu
		tEz=>TexEz_gpu
		tBx=>TexBx_gpu
		tBy=>TexBy_gpu
		tBz=>TexBz_gpu
		

		
		!if((inc_axis.eq..true.).and.(procxind.eq.0)) then
		if((inc_axis.eq..true.).and.(procxind.le.0)) then
			do kc=1,Nchunk_prtl_gpu
				indi=(kc-1)*chunk_size_prtl_gpu+1
				indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			       call MoveDepositAxisPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(prtl_arr_size,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,qp_gpu,tagp_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
			                     TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,Jx_gpu,Jy_gpu,Jz_gpu,mx_gpu,my_gpu,mz_gpu,indi,indf,rborders(procxind)-xmin,inv_dtheta,ymin,ymax,ylen,ny)
		 
			end do 
		else 

			do kc=1,Nchunk_prtl_gpu
				indi=(kc-1)*chunk_size_prtl_gpu+1
				indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			       call MoveDepositPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(prtl_arr_size,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,qp_gpu,tagp_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
			                     TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,Jx_gpu,Jy_gpu,Jz_gpu,mx_gpu,my_gpu,mz_gpu,indi,indf,rborders(procxind)-xmin,inv_dtheta)
		 
			end do 
		
		end if 
			
		
		nullify(tEx)
		nullify(tEy)
		nullify(tEz)
		nullify(tBx)
		nullify(tBy)
		nullify(tBz)
	end subroutine MoveDepositPrtlGPU
	
	
	subroutine MoveTestPrtlGPU

	end subroutine MoveTestPrtlGPU
	
	
	attributes(global) subroutine MoveDepositPrtlKernel(size,xp,yp,zp,up,vp,wp,flvp,qp,tagp,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,mx,my,mz,indi,indf,r_proc,inv_dtheta)
	    integer, value :: size
		real, dimension(size) ::xp,yp,zp,up,vp,wp,qp!,var1p
		integer, dimension(size) :: flvp,tagp
		real, dimension(*)  ::flvrqm
		integer             :: mx,my,mz
		real, value :: c,sqc,cinv,qi,r_proc,inv_dtheta 
		integer, value :: indi,indf		
		real, dimension(mx,my,mz)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
		
		integer :: n 
		integer :: stat
		!integer :: InsertAt ! used in prtl outliers
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g,x1,y1,r1
		real :: pEx,pEy,pEz,pBx,pBy,pBz
		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: wt 
#else 
        real,dimension(8) :: wt
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
		  if(qp(n).eq.0) return
		  !if(xp(n)+r_proc.le.1.0) return
		  !tagp(n)=0

		  x0=xp(n)
          y0=yp(n)
          z0=zp(n)

		  !load q/m for each particles into a vector
 	      if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
		    qm=0
	      else
	        qm=flvrqm(flvp(n))*0.5 !for optimisation 0.5 is multiplied here itself
	      end if

		  !interpolation
		  i=floor(x0)
		  j=floor(y0)
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=floor(z0)
          dz=z0-k
	      wt(1)=(1.0-dx)*(1.0-dy)*(1.0-dz)
	      wt(2)=dx          *(1.0-dy)*(1.0-dz)
	      wt(3)=(1.0-dx) *dy          *(1.0-dz)
	      wt(4)=dx          *dy          *(1.0-dz)
	      wt(5)=(1.0-dx)*(1.0-dy)*dz
	      wt(6)=dx          *(1.0-dy)*dz
	      wt(7)=(1.0-dx)*dy          *dz
	      wt(8)=dx          *dy          *dz
#else
          k=1
		  wt(1)=(1.0-dx)*(1.0-dy)
		  wt(2)=dx          *(1.0-dy)
		  wt(3)=(1.0-dx) *dy
		  wt(4)=dx          *dy
#endif


#ifdef twoD
          pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)
		  pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)
		  pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)
          pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)
		  pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)
		  pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)
#else
          pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)+wt(5)*Ex(i,j,k+1)+wt(6)*Ex(i+1,j,k+1)+wt(7)*Ex(i,j+1,k+1)+wt(8)*Ex(i+1,j+1,k+1)
          pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)+wt(5)*Ey(i,j,k+1)+wt(6)*Ey(i+1,j,k+1)+wt(7)*Ey(i,j+1,k+1)+wt(8)*Ey(i+1,j+1,k+1)
          pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)+wt(5)*Ez(i,j,k+1)+wt(6)*Ez(i+1,j,k+1)+wt(7)*Ez(i,j+1,k+1)+wt(8)*Ez(i+1,j+1,k+1)
          pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)+wt(5)*Bx(i,j,k+1)+wt(6)*Bx(i+1,j,k+1)+wt(7)*Bx(i,j+1,k+1)+wt(8)*Bx(i+1,j+1,k+1)
          pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)+wt(5)*By(i,j,k+1)+wt(6)*By(i+1,j,k+1)+wt(7)*By(i,j+1,k+1)+wt(8)*By(i+1,j+1,k+1)
          pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)+wt(5)*Bz(i,j,k+1)+wt(6)*Bz(i+1,j,k+1)+wt(7)*Bz(i,j+1,k+1)+wt(8)*Bz(i+1,j+1,k+1)
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

               g=1.0/sqrt(sqc+u0**2+v0**2+w0**2)   ! 1/c*gamma
               pBx=g*pBx
               pBy=g*pBy
               pBz=g*pBz

               f=2.0/(1.0+pBx*pBx+pBy*pBy+pBz*pBz)
               u1=(u0+v0*pBz-w0*pBy)*f
               v1=(v0+w0*pBx-u0*pBz)*f
               w1=(w0+u0*pBy-v0*pBx)*f

               u0=u0+v1*pBz-w1*pBy+pEx
               v0=v0+w1*pBx-u1*pBz+pEy
               w0=w0+u1*pBy-v1*pBx+pEz

              			 
               u1=u0*cinv
               v1=v0*cinv
               wp(n)=w0*cinv

               g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

               x1=xp(n) + u1*g + r_proc
			   y1=v1*g
			   r1=sqrt(x1**2+y1**2)
			   xp(n)=r1 - r_proc
               yp(n)=yp(n) + atan(y1/x1)*inv_dtheta
               zp(n)=zp(n) + wp(n)*g
			   
			   r1=1.0/r1
			   up(n)=(x1*u1+y1*v1)*r1
	           vp(n)=(-y1*u1+x1*v1)*r1
			   
			   		   

          qthis=qp(n)*qi   ! q = particle's weight X sign of the charge


!current deposition
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
               xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),0.5*(x0+xp(n))))
               yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),0.5*(y0+yp(n))))
               zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),0.5*(z0+zp(n))))
          Fx1=qthis*(xr-x0)
          Fy1=qthis*(yr-y0)
          Fz1=qthis*(zr-z0)

          Wx1=0.5*(x0+xr)-i1
          Wy1=0.5*(y0+yr)-j1

          Wx2=0.5*(xp(n)+xr)-i2
          Wy2=0.5*(yp(n)+yr)-j2

#ifdef twoD
          Wz1=0.0
          Wz2=0.0
#else
          Wz1=0.5*(z0+zr)-k1
          Wz2=0.5*(zp(n)+zr)-k2
#endif


          Fx2=qthis*(xp(n)-xr)
          Fy2=qthis*(yp(n)-yr)
          Fz2=qthis*(zp(n)-zr)


          Jx1(1)=Fx1 * (1.0-Wy1)*(1.0-Wz1)
          Jx1(2)=Fx1 *  Wy1    *(1.0-Wz1)
#ifndef twoD
              Jx1(3)= Fx1 * (1.0-Wy1) * Wz1
              Jx1(4)= Fx1 *  Wy1    * Wz1
#endif

          Jx2(1)=Fx2 * (1.0-Wy2)*(1.0-Wz2)
          Jx2(2)=Fx2 *  Wy2    *(1.0-Wz2)
#ifndef twoD
               Jx2(3)=Fx2 * (1.0-Wy2)* Wz2
               Jx2(4)=Fx2 *  Wy2    * Wz2
#endif


          Jy1(1)=Fy1 * (1.0-Wx1)*(1.0-Wz1)
          Jy1(2)=Fy1 *  Wx1    *(1.0-Wz1)
#ifndef twoD
               Jy1(3)=Fy1 * (1.0-Wx1)* Wz1
               Jy1(4)=Fy1 *  Wx1    * Wz1
#endif

          Jy2(1)=Fy2 * (1.0-Wx2)*(1.0-Wz2)
          Jy2(2)=Fy2 *  Wx2    *(1.0-Wz2)
#ifndef twoD
             Jy2(3)= Fy2 * (1.0-Wx2)* Wz2
             Jy2(4)= Fy2 *  Wx2    * Wz2
#endif


          Jz1(1)=Fz1 * (1.0-Wx1)*(1.0-Wy1)
          Jz1(2)=Fz1 *  Wx1    *(1.0-Wy1)
          Jz1(3)=Fz1 * (1.0-Wx1)* Wy1
          Jz1(4)=Fz1 *  Wx1    * Wy1

          Jz2(1)=Fz2 * (1.0-Wx2)*(1.0-Wy2)
          Jz2(2)=Fz2 *  Wx2    *(1.0-Wy2)
          Jz2(3)=Fz2 * (1.0-Wx2)* Wy2
          Jz2(4)=Fz2 *  Wx2    * Wy2



          stat=atomicAdd(Jx(i1,j1,  k1  ), Jx1(1) )
          stat=atomicAdd(Jx(i1,j1+1,k1  ), Jx1(2) )
#ifndef twoD
               stat=atomicAdd(Jx(i1,j1,  k1+1), Jx1(3))
               stat=atomicAdd(Jx(i1,j1+1,k1+1), Jx1(4))
#endif

          stat=atomicAdd(Jx(i2,j2,  k2  ), Jx2(1) )
          stat=atomicAdd(Jx(i2,j2+1,k2  ), Jx2(2) )
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
 
              
! Periodic Boundary conditions 
#ifdef twoD
             if(zp(n).ge.2.0) then
               zp(n)=zp(n)-1.0
             else if(zp(n).lt.1.0) then
               zp(n)=1.0+zp(n)
             end if
#endif
	  

end subroutine MoveDepositPrtlKernel


	attributes(global) subroutine MoveDepositAxisPrtlKernel(size,xp,yp,zp,up,vp,wp,flvp,qp,tagp,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,mx,my,mz,indi,indf,r_proc,inv_dtheta,ymin,ymax,ylen,ny)
	    integer, value :: size
		real, dimension(size) ::xp,yp,zp,up,vp,wp,qp!,var1p
		integer, dimension(size) :: flvp,tagp
		real, dimension(*)  ::flvrqm
		integer             :: mx,my,mz
		real, value :: c,sqc,cinv,qi,r_proc,inv_dtheta, ymin, ymax, ylen
		integer, value :: indi,indf, ny		
		real, dimension(mx,my,mz)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
		
		integer :: n 
		integer :: stat
		!integer :: InsertAt ! used in prtl outliers
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g,x1,y1,r1, t, dt, delt
		real :: pEx,pEy,pEz,pBx,pBy,pBz
		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: wt 
#else 
        real,dimension(8) :: wt
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
		real :: alpha, beta, dt_temp
								
			
		  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		  if(n.gt.indf) return
		  if(qp(n).eq.0) return
		  !if(xp(n)+r_proc.gt.1.0) return
 		  t=0.0
 		  dt=1.0
          !tagp(n)=0
          do while(t.lt.1.0)
			  !tagp(n)=tagp(n)+1 
		  x0=xp(n)
          y0=yp(n)
          z0=zp(n)

		  !load q/m for each particles into a vector
 	      if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
		    qm=0
	      else
	        qm=flvrqm(flvp(n))*0.5 !for optimisation 0.5 is multiplied here itself
	      end if

		  !interpolation
		  i=floor(x0)
		  j=floor(y0)
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=floor(z0)
          dz=z0-k
	      wt(1)=(1.0-dx)*(1.0-dy)*(1.0-dz)
	      wt(2)=dx          *(1.0-dy)*(1.0-dz)
	      wt(3)=(1.0-dx) *dy          *(1.0-dz)
	      wt(4)=dx          *dy          *(1.0-dz)
	      wt(5)=(1.0-dx)*(1.0-dy)*dz
	      wt(6)=dx          *(1.0-dy)*dz
	      wt(7)=(1.0-dx)*dy          *dz
	      wt(8)=dx          *dy          *dz
#else
          k=1
		  wt(1)=(1.0-dx)*(1.0-dy)
		  wt(2)=dx          *(1.0-dy)
		  wt(3)=(1.0-dx) *dy
		  wt(4)=dx          *dy
#endif


#ifdef twoD
          pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)
		  pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)
		  pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)
          pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)
		  pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)
		  pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)
#else
          pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)+wt(5)*Ex(i,j,k+1)+wt(6)*Ex(i+1,j,k+1)+wt(7)*Ex(i,j+1,k+1)+wt(8)*Ex(i+1,j+1,k+1)
          pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)+wt(5)*Ey(i,j,k+1)+wt(6)*Ey(i+1,j,k+1)+wt(7)*Ey(i,j+1,k+1)+wt(8)*Ey(i+1,j+1,k+1)
          pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)+wt(5)*Ez(i,j,k+1)+wt(6)*Ez(i+1,j,k+1)+wt(7)*Ez(i,j+1,k+1)+wt(8)*Ez(i+1,j+1,k+1)
          pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)+wt(5)*Bx(i,j,k+1)+wt(6)*Bx(i+1,j,k+1)+wt(7)*Bx(i,j+1,k+1)+wt(8)*Bx(i+1,j+1,k+1)
          pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)+wt(5)*By(i,j,k+1)+wt(6)*By(i+1,j,k+1)+wt(7)*By(i,j+1,k+1)+wt(8)*By(i+1,j+1,k+1)
          pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)+wt(5)*Bz(i,j,k+1)+wt(6)*Bz(i+1,j,k+1)+wt(7)*Bz(i,j+1,k+1)+wt(8)*Bz(i+1,j+1,k+1)
#endif

            delt=dt
			if(t+dt.gt.1.0) delt=1.0-t 
			qm=qm*delt
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

               g=1.0/sqrt(sqc+u0**2+v0**2+w0**2)   ! 1/c*gamma
               pBx=g*pBx
               pBy=g*pBy
               pBz=g*pBz

               f=2.0/(1.0+pBx*pBx+pBy*pBy+pBz*pBz)
               u1=(u0+v0*pBz-w0*pBy)*f
               v1=(v0+w0*pBx-u0*pBz)*f
               w1=(w0+u0*pBy-v0*pBx)*f

               u0=u0+v1*pBz-w1*pBy+pEx
               v0=v0+w1*pBx-u1*pBz+pEy
               w0=w0+u1*pBy-v1*pBx+pEz

              			 
	           u1=u0*cinv
	           v1=v0*cinv

	           g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

	         x1=xp(n) + u1*g*delt - 3.5
	  	     y1=v1*g*delt
	  	     r1=sqrt(x1**2+y1**2)
		 
	  		 if(x1.gt.0) then 
	  		     dy= atan(y1/x1)*inv_dtheta
	  		 else 
	  			 dy= 0.0
	  		 end if 
	  		 !print*,'n',n,'dt',dt,'dy',dy,'u1',u1,'v1',v1
	  		 if(abs(dy).gt.0.45) then 
				 !dt=min(dt*0.5,(0.4/dy)*dt)
				 
				 alpha=tan(0.4/inv_dtheta)
				 if(v1.lt.0) alpha=-alpha

				 beta=abs(g*(v1-alpha*u1))
				 if(beta.ne.0) then
				     dt_temp=alpha*(xp(n)-3.5)/beta
				 else
					 dt_temp=1.0
				 end if
				 dt=min(abs(dt_temp),dt*0.5)
				 
				 cycle 
			 end if 
			 t=t+delt
			 dt=1.0
		 		
	  		 if(x1.gt.0) then 	 
	  	         xp(n)=r1 + 3.5
	  	         up(n)=(x1*u1+y1*v1)/r1
	             vp(n)=(-y1*u1+x1*v1)/r1	
	  		 else
	  			 xp(n)=x1 + 3.5
	  	         up(n)=u1
	  			 vp(n)=v1
	  		 end if  
	         yp(n)=yp(n)+dy
	
	  		 wp(n)=w0*cinv
	  		 zp(n)=zp(n) + wp(n)*g*delt
			   
			   		   

          qthis=qp(n)*qi   ! q = particle's weight X sign of the charge


!current deposition
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
               xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),0.5*(x0+xp(n))))
               yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),0.5*(y0+yp(n))))
               zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),0.5*(z0+zp(n))))
          Fx1=qthis*(xr-x0)
          Fy1=qthis*(yr-y0)
          Fz1=qthis*(zr-z0)

          Wx1=0.5*(x0+xr)-i1
          Wy1=0.5*(y0+yr)-j1

          Wx2=0.5*(xp(n)+xr)-i2
          Wy2=0.5*(yp(n)+yr)-j2

#ifdef twoD
          Wz1=0.0
          Wz2=0.0
#else
          Wz1=0.5*(z0+zr)-k1
          Wz2=0.5*(zp(n)+zr)-k2
#endif


          Fx2=qthis*(xp(n)-xr)
          Fy2=qthis*(yp(n)-yr)
          Fz2=qthis*(zp(n)-zr)


          Jx1(1)=Fx1 * (1.0-Wy1)*(1.0-Wz1)
          Jx1(2)=Fx1 *  Wy1    *(1.0-Wz1)
#ifndef twoD
              Jx1(3)= Fx1 * (1.0-Wy1) * Wz1
              Jx1(4)= Fx1 *  Wy1    * Wz1
#endif

          Jx2(1)=Fx2 * (1.0-Wy2)*(1.0-Wz2)
          Jx2(2)=Fx2 *  Wy2    *(1.0-Wz2)
#ifndef twoD
               Jx2(3)=Fx2 * (1.0-Wy2)* Wz2
               Jx2(4)=Fx2 *  Wy2    * Wz2
#endif


          Jy1(1)=Fy1 * (1.0-Wx1)*(1.0-Wz1)
          Jy1(2)=Fy1 *  Wx1    *(1.0-Wz1)
#ifndef twoD
               Jy1(3)=Fy1 * (1.0-Wx1)* Wz1
               Jy1(4)=Fy1 *  Wx1    * Wz1
#endif

          Jy2(1)=Fy2 * (1.0-Wx2)*(1.0-Wz2)
          Jy2(2)=Fy2 *  Wx2    *(1.0-Wz2)
#ifndef twoD
             Jy2(3)= Fy2 * (1.0-Wx2)* Wz2
             Jy2(4)= Fy2 *  Wx2    * Wz2
#endif


          Jz1(1)=Fz1 * (1.0-Wx1)*(1.0-Wy1)
          Jz1(2)=Fz1 *  Wx1    *(1.0-Wy1)
          Jz1(3)=Fz1 * (1.0-Wx1)* Wy1
          Jz1(4)=Fz1 *  Wx1    * Wy1

          Jz2(1)=Fz2 * (1.0-Wx2)*(1.0-Wy2)
          Jz2(2)=Fz2 *  Wx2    *(1.0-Wy2)
          Jz2(3)=Fz2 * (1.0-Wx2)* Wy2
          Jz2(4)=Fz2 *  Wx2    * Wy2



          stat=atomicAdd(Jx(i1,j1,  k1  ), Jx1(1) )
          stat=atomicAdd(Jx(i1,j1+1,k1  ), Jx1(2) )
#ifndef twoD
               stat=atomicAdd(Jx(i1,j1,  k1+1), Jx1(3))
               stat=atomicAdd(Jx(i1,j1+1,k1+1), Jx1(4))
#endif

          stat=atomicAdd(Jx(i2,j2,  k2  ), Jx2(1) )
          stat=atomicAdd(Jx(i2,j2+1,k2  ), Jx2(2) )
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
 
              
! Periodic Boundary conditions 
#ifdef twoD
             if(zp(n).ge.2.0) then
               zp(n)=zp(n)-1.0
             else if(zp(n).lt.1.0) then
               zp(n)=1.0+zp(n)
             end if
#endif


			if(xp(n).le.3.5) then
			  xp(n)= 3.5 + (3.5-xp(n))
			  yp(n)= yp(n) + ny/2.0
			  up(n) = - up(n)
			  vp(n) = - vp(n)   
			end if 
			if(yp(n).gt.ymax) yp(n)=yp(n)-ylen
			if(yp(n).lt.ymin) yp(n)=yp(n)+ylen
      
	  end do! end of the while loop  

end subroutine MoveDepositAxisPrtlKernel



end module movdep_gpu 
	