module movdep_gpu 
	use parameters 
	use vars 
	use var_gpu 
	use fields_gpu !Fitlering subroutines are in fields_gpu
	use cudafor
	use cudadevice
	real, texture, pointer, dimension(:,:,:) :: tEx,tEy,tEz,tBx,tBy,tBz

	
contains 
	subroutine MoveDepositPrtlGPU
		integer :: kc, indi,indf
		!integer :: np_send_gpu_temp
		

        call ResetVecCurrentGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(VecJ_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)			
		call InitTexExGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexEx_gpu,Ex_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)
		call InitTexEyGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexEy_gpu,Ey_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)
		call InitTexEzGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexEz_gpu,Ez_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu)
		call InitTexBxGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexBx_gpu,Bx_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,Bx_ext0)
		call InitTexByGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexBy_gpu,By_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,By_ext0)
		call InitTexBzGPUKernel<<<tGrid_gpu_global,tBlock_gpu_global>>>(TexBz_gpu,Bz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,Bz_ext0)
		
		tEx=>TexEx_gpu
		tEy=>TexEy_gpu
		tEz=>TexEz_gpu
		tBx=>TexBx_gpu
		tBy=>TexBy_gpu
		tBz=>TexBz_gpu
		
				
		
		!np_send_gpu_temp=0 
		!np_send_gpu=np_send_gpu_temp
		
		do kc=1,Nchunk_prtl_gpu
			indi=(kc-1)*chunk_size_prtl_gpu+1
			indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
		       call MoveDepositPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,qp_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
		                     TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,Jx_gpu,Jy_gpu,Jz_gpu,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,indi,indf)
			  ! call MoveDepositPrtlKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,qp_gpu,var1p_gpu,tagp_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
			  !		         TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,Jx_gpu,Jy_gpu,Jz_gpu,Bx_ext0,By_ext0,Bz_ext0,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,indi,indf,&
			!				 qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu,flvp_send_gpu,tagp_send_gpu)			 
		end do 
		!np_recv_host=np_send_gpu 
		
		
		
		nullify(tEx)
		nullify(tEy)
		nullify(tEz)
		nullify(tBx)
		nullify(tBy)
		nullify(tBz)
	end subroutine MoveDepositPrtlGPU
	
	
	subroutine MoveTestPrtlGPU
		integer :: kc, indi,indf
#ifdef GPU_EXCLUSIVE		
		call SetFilteredEfldGPU_Exclusive
#endif		
		do kc=1,Nchunk_test_prtl_gpu
			indi=(kc-1)*chunk_size_test_prtl_gpu+1
			indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
	       call MoveTestPrtlKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,flvtp_gpu,qtp_gpu,flvrqm_gpu,c,sqc,cinv,qi,&
	                     FilteredEx_gpu,FilteredEy_gpu,FilteredEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,Bx_ext0,By_ext0,Bz_ext0,xmin1_gpu,xmax1_gpu,ymin1_gpu,ymax1_gpu,zmin1_gpu,zmax1_gpu,indi,indf)			
		end do 
	end subroutine MoveTestPrtlGPU
	
	
#ifndef GPU_USE_INTRINSICS	
	
	!attributes(global) subroutine  MoveDepositPrtlKernel(xp,yp,zp,up,vp,wp,flvp,qp,var1p,tagp,psize,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,Bx_ext0,By_ext0,Bz_ext0,x1,x2,y1,y2,z1,z2,indi,indf,&
		!qout,xout,yout,zout,uout,vout,wout,var1out,flvout,tagout)
	attributes(global) subroutine MoveDepositPrtlKernel(xp,yp,zp,up,vp,wp,flvp,qp,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,x1,x2,y1,y2,z1,z2,indi,indf)
		real, dimension(:) ::xp,yp,zp,up,vp,wp,qp!,var1p
		integer, dimension(:) :: flvp!,tagp
		real, dimension(:)  ::flvrqm
		real, value :: c,sqc,cinv,qi !should put then on gpu as constants
		integer :: x1,x2,y1,y2,z1,z2
		integer, value :: indi,indf
#ifndef twoD  		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#else 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,1:1)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#endif			
		
		!real, dimension(:) ::xout,yout,zout,uout,vout,wout,qout,var1out
		!integer, dimension(:) :: flvout,tagout
		
		integer :: n 
		integer :: stat
		!integer :: InsertAt ! used in prtl outliers
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g
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

			  x0=xp(n)
	          y0=yp(n)
	          z0=zp(n)

			 ! load q/m for each particles into a vector
     	    if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
			  qm=0
		    else
  		      qm=flvrqm(flvp(n))*0.5 !for optimisation 0.5 is multiplied here itself
		    end if

			!interpolation
		  i=x0
		  j=y0
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=z0
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

               up(n)=u0*cinv
               vp(n)=v0*cinv
               wp(n)=w0*cinv

               g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g
               yp(n)=yp(n) + vp(n)*g
               zp(n)=zp(n) + wp(n)*g


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

! !separate the particles that are leaving the GPU domain [led to slow down, possible due to atomic operation]
! #ifndef twoD
!          if((xp(n).lt.x1).or.(xp(n).gt.x2).or.(yp(n).lt.y1).or.(yp(n).gt.y2).or.(zp(n).lt.z1).or.(zp(n).gt.z2)) then
! #else
!          if((xp(n).lt.x1).or.(xp(n).gt.x2).or.(yp(n).lt.y1).or.(yp(n).gt.y2)) then
! #endif
! 						 InsertAt=atomicinc(np_send_gpu,1147483647)
! 						 InsertAt=InsertAt+1
! 						 qout(InsertAt)=qp(n)
! 					     xout(InsertAt)=xp(n)
! 						 yout(InsertAt)=yp(n)
! 						 zout(InsertAt)=zp(n)
! 						 uout(InsertAt)=up(n)
! 						 vout(InsertAt)=vp(n)
! 						 wout(InsertAt)=wp(n)
! 						 var1out(InsertAt)=var1p(n)
! 						 flvout(InsertAt)=flvp(n)
! 						 tagout(InsertAt)=tagp(n)
!
! 						 ! delete the particle
! 						 qp(n)=0
! 						 xp(n)=x1+0.5
! 						 yp(n)=y1+0.5
! 						 zp(n)=z1+0.5
! 						 up(n)=0
! 						 vp(n)=0
! 						 wp(n)=0
! 						 var1p(n)=0
! 						 tagp(n)=0
! 						 flvp(n)=0
! 		 end if			  

end subroutine MoveDepositPrtlKernel

#endif


!---------------------------------------------------------------------------
! Test Particle Mover 
!---------------------------------------------------------------------------
attributes(global) subroutine MoveTestPrtlKernel(xp,yp,zp,up,vp,wp,flvp,qp,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Bx_ext0,By_ext0,Bz_ext0,x1,x2,y1,y2,z1,z2,indi,indf)
		real, dimension(:) ::xp,yp,zp,up,vp,wp,qp
		integer, dimension(:) :: flvp
		real, dimension(:)  ::flvrqm
		real, value :: c,sqc,cinv,qi
		integer :: x1,x2,y1,y2,z1,z2
		integer, value :: indi,indf
#ifndef twoD  		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2)  :: Ex,Ey,Ez,Bx,By,Bz
#else 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,1:1)  :: Ex,Ey,Ez,Bx,By,Bz
#endif			
		real, value :: Bx_ext0,By_ext0,Bz_ext0
		
		integer :: n 
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g
		real :: pEx,pEy,pEz,pBx,pBy,pBz
		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: wt 
#else 
        real,dimension(8) :: wt
#endif 	


			
			  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
			  if(n.gt.indf) return 
			  
			  x0=xp(n)
	          y0=yp(n)
	          z0=zp(n)

			 ! load q/m for each particles into a vector
     	    if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
			  qm=0
		    else
  		      qm=flvrqm(flvp(n))*0.5 !for optimisation 0.5 is multiplied here itself
		    end if

			!interpolation
		  i=x0
		  j=y0
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=z0
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

               up(n)=u0*cinv
               vp(n)=v0*cinv
               wp(n)=w0*cinv

               g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g
               yp(n)=yp(n) + vp(n)*g
               zp(n)=zp(n) + wp(n)*g
			   
#ifdef twoD
            if(zp(n).ge.2.0) then
              zp(n)=zp(n)-1.0
            else if(zp(n).lt.1.0) then
              zp(n)=1.0+zp(n)
            end if 
#endif 			   
	
end subroutine MoveTestPrtlKernel




!---------------------------------------------------------------------------
! Optimized Particle Mover 
!---------------------------------------------------------------------------

#ifdef GPU_USE_INTRINSICS	
	attributes(global) subroutine MoveDepositPrtlKernel(xp,yp,zp,up,vp,wp,flvp,qp,flvrqm,c,sqc,cinv,qi,Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,x1,x2,y1,y2,z1,z2,indi,indf)
		real, dimension(:) ::xp,yp,zp,up,vp,wp,qp!,var1p
		integer, dimension(:) :: flvp!,tagp
		real, dimension(:)  ::flvrqm
		real, value :: c,sqc,cinv,qi
		integer :: x1,x2,y1,y2,z1,z2
		integer, value :: indi,indf
#ifndef twoD  		
		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#else 		
		real, dimension(x1-2:x2+2,y1-2:y2+2,1:1)  :: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
#endif			
		
		!real, dimension(:) ::xout,yout,zout,uout,vout,wout,qout,var1out
		!integer, dimension(:) :: flvout,tagout
		
		integer :: n 
		integer :: stat
		!integer :: InsertAt ! used in prtl outliers
		!varaibles used in moving particles 
		real :: x0,y0,z0,qm,u0,v0,w0,u1,v1,w1,f,g
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
        real,dimension(4) :: Jz1,Jz2,recvJ
		
		integer :: cell_index,other_cell_index
		integer :: lane
		integer :: unclaimed
		logical :: is_peer
		integer :: peers
		integer :: first_lane
		integer :: rel_pos
		integer :: next
		integer :: done
		logical :: not_done
		
								
			
			  n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
			  !if(n.gt.indf) return 
			  if(n.gt.indf) return 
			  
			  
			  
			  x0=xp(n)
	          y0=yp(n)
	          z0=zp(n)

			 ! load q/m for each particles into a vector
     	    if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0
			  qm=0
		    else
  		      qm=flvrqm(flvp(n))*0.5 !for optimisation 0.5 is multiplied here itself
		    end if

			!interpolation
		  i=x0
		  j=y0
		  dx=x0-i
		  dy=y0-j
#ifndef twoD
          k=z0
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


! #ifdef twoD
!           pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)
! 		  pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)
! 		  pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)
!           pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)
! 		  pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)
! 		  pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)
! #else
!           pEx=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)+wt(5)*Ex(i,j,k+1)+wt(6)*Ex(i+1,j,k+1)+wt(7)*Ex(i,j+1,k+1)+wt(8)*Ex(i+1,j+1,k+1)
!           pEy=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)+wt(5)*Ey(i,j,k+1)+wt(6)*Ey(i+1,j,k+1)+wt(7)*Ey(i,j+1,k+1)+wt(8)*Ey(i+1,j+1,k+1)
!           pEz=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)+wt(5)*Ez(i,j,k+1)+wt(6)*Ez(i+1,j,k+1)+wt(7)*Ez(i,j+1,k+1)+wt(8)*Ez(i+1,j+1,k+1)
!           pBx=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)+wt(5)*Bx(i,j,k+1)+wt(6)*Bx(i+1,j,k+1)+wt(7)*Bx(i,j+1,k+1)+wt(8)*Bx(i+1,j+1,k+1)
!           pBy=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)+wt(5)*By(i,j,k+1)+wt(6)*By(i+1,j,k+1)+wt(7)*By(i,j+1,k+1)+wt(8)*By(i+1,j+1,k+1)
!           pBz=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)+wt(5)*Bz(i,j,k+1)+wt(6)*Bz(i+1,j,k+1)+wt(7)*Bz(i,j+1,k+1)+wt(8)*Bz(i+1,j+1,k+1)
! #endif



#ifdef twoD
          pEx=wt(1)*tEx(i,j,k)+wt(2)*tEx(i+1,j,k)+wt(3)*tEx(i,j+1,k)+wt(4)*tEx(i+1,j+1,k)
		  pEy=wt(1)*tEy(i,j,k)+wt(2)*tEy(i+1,j,k)+wt(3)*tEy(i,j+1,k)+wt(4)*tEy(i+1,j+1,k)
		  pEz=wt(1)*tEz(i,j,k)+wt(2)*tEz(i+1,j,k)+wt(3)*tEz(i,j+1,k)+wt(4)*tEz(i+1,j+1,k)
          pBx=wt(1)*tBx(i,j,k)+wt(2)*tBx(i+1,j,k)+wt(3)*tBx(i,j+1,k)+wt(4)*tBx(i+1,j+1,k)
		  pBy=wt(1)*tBy(i,j,k)+wt(2)*tBy(i+1,j,k)+wt(3)*tBy(i,j+1,k)+wt(4)*tBy(i+1,j+1,k)
		  pBz=wt(1)*tBz(i,j,k)+wt(2)*tBz(i+1,j,k)+wt(3)*tBz(i,j+1,k)+wt(4)*tBz(i+1,j+1,k)
#else
          pEx=wt(1)*tEx(i,j,k)+wt(2)*tEx(i+1,j,k)+wt(3)*tEx(i,j+1,k)+wt(4)*tEx(i+1,j+1,k)+wt(5)*tEx(i,j,k+1)+wt(6)*tEx(i+1,j,k+1)+wt(7)*tEx(i,j+1,k+1)+wt(8)*tEx(i+1,j+1,k+1)
          pEy=wt(1)*tEy(i,j,k)+wt(2)*tEy(i+1,j,k)+wt(3)*tEy(i,j+1,k)+wt(4)*tEy(i+1,j+1,k)+wt(5)*tEy(i,j,k+1)+wt(6)*tEy(i+1,j,k+1)+wt(7)*tEy(i,j+1,k+1)+wt(8)*tEy(i+1,j+1,k+1)
          pEz=wt(1)*tEz(i,j,k)+wt(2)*tEz(i+1,j,k)+wt(3)*tEz(i,j+1,k)+wt(4)*tEz(i+1,j+1,k)+wt(5)*tEz(i,j,k+1)+wt(6)*tEz(i+1,j,k+1)+wt(7)*tEz(i,j+1,k+1)+wt(8)*tEz(i+1,j+1,k+1)
          pBx=wt(1)*tBx(i,j,k)+wt(2)*tBx(i+1,j,k)+wt(3)*tBx(i,j+1,k)+wt(4)*tBx(i+1,j+1,k)+wt(5)*tBx(i,j,k+1)+wt(6)*tBx(i+1,j,k+1)+wt(7)*tBx(i,j+1,k+1)+wt(8)*tBx(i+1,j+1,k+1)
          pBy=wt(1)*tBy(i,j,k)+wt(2)*tBy(i+1,j,k)+wt(3)*tBy(i,j+1,k)+wt(4)*tBy(i+1,j+1,k)+wt(5)*tBy(i,j,k+1)+wt(6)*tBy(i+1,j,k+1)+wt(7)*tBy(i,j+1,k+1)+wt(8)*tBy(i+1,j+1,k+1)
          pBz=wt(1)*tBz(i,j,k)+wt(2)*tBz(i+1,j,k)+wt(3)*tBz(i,j+1,k)+wt(4)*tBz(i+1,j+1,k)+wt(5)*tBz(i,j,k+1)+wt(6)*tBz(i+1,j,k+1)+wt(7)*tBz(i,j+1,k+1)+wt(8)*tBz(i+1,j+1,k+1)	  
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

               up(n)=u0*cinv
               vp(n)=v0*cinv
               wp(n)=w0*cinv

               g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g
               yp(n)=yp(n) + vp(n)*g
               zp(n)=zp(n) + wp(n)*g


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



		  !------------------------------------------------
		  ! Frist warp level deposition
		  !-----------------------------------------------
#ifdef twoD
          cell_index=j1*(x2-x1+5)+i1
#else
          cell_index=k1*(y2-y1+5)*(x2-x1+5) + j1*(x2-x1+5) + i1
#endif

          lane=mod((threadIdx%x-1), 32) + 1 !WARP_SIZE is assumed to be 32
          unclaimed= B'11111111111111111111111111111111'
		  is_peer=.false.
		  do while(.not.is_peer)
              other_cell_index=__SHFL(cell_index,__FFS(unclaimed))
		      is_peer=cell_index.eq.other_cell_index
			  peers=ballot(is_peer)
			  unclaimed=XOR(unclaimed,peers)
	      end do

		  first_lane=__FFS(peers)
		  rel_pos=__popc(IBITS(peers,0,lane))-1
          peers=IAND(peers,ISHFT(B'11111111111111111111111111111110',lane-1))
		  do while(anythread((peers.ne.0)))

			  next=__FFS(peers)

			  !Copy data from other threads in the warp
			  recvJ(1)=__SHFL(Jx1(1),next)
			  recvJ(2)=__SHFL(Jx1(2),next)
#ifndef twoD
			  recvJ(3)=__SHFL(Jx1(3),next)
			  recvJ(4)=__SHFL(Jx1(4),next)
#endif
			  if(next.gt.0) then
				  Jx1(1)=Jx1(1)+recvJ(1)
				  Jx1(2)=Jx1(2)+recvJ(2)
#ifndef twoD
				  Jx1(3)=Jx1(3)+recvJ(3)
				  Jx1(4)=Jx1(4)+recvJ(4)
#endif
			  end if

			  recvJ(1)=__SHFL(Jy1(1),next)
			  recvJ(2)=__SHFL(Jy1(2),next)
#ifndef twoD
			  recvJ(3)=__SHFL(Jy1(3),next)
			  recvJ(4)=__SHFL(Jy1(4),next)
#endif
			  if(next.gt.0) then
				  Jy1(1)=Jy1(1)+recvJ(1)
				  Jy1(2)=Jy1(2)+recvJ(2)
#ifndef twoD
				  Jy1(3)=Jy1(3)+recvJ(3)
				  Jy1(4)=Jy1(4)+recvJ(4)
#endif
			  end if

			  recvJ(1)=__SHFL(Jz1(1),next)
			  recvJ(2)=__SHFL(Jz1(2),next)
			  recvJ(3)=__SHFL(Jz1(3),next)
			  recvJ(4)=__SHFL(Jz1(4),next)
			  if(next.gt.0) then
				  Jz1(1)=Jz1(1)+recvJ(1)
				  Jz1(2)=Jz1(2)+recvJ(2)
				  Jz1(3)=Jz1(3)+recvJ(3)
				  Jz1(4)=Jz1(4)+recvJ(4)
			  end if
			  !----------------------------

			  done=IAND(rel_pos,1)
			  if(done.eq.1) then
				   not_done=.false.
			  else
				   not_done=.true.
			  end if
			  peers=IAND(peers,ballot(not_done))
			  !peers=IAND(peers,ballot(NOT(done)))
			  rel_pos=ISHFT(rel_pos,-1)
		  end do




		  if(lane.eq.first_lane) then
 			  !if(lane.eq.1) then


          stat=atomicAdd(Jx(i1,j1,  k1  ), Jx1(1) )
          stat=atomicAdd(Jx(i1,j1+1,k1  ), Jx1(2) )
#ifndef twoD
               stat=atomicAdd(Jx(i1,j1,  k1+1), Jx1(3))
               stat=atomicAdd(Jx(i1,j1+1,k1+1), Jx1(4))
#endif


          stat=atomicAdd(Jy(i1  ,j1,k1  ), Jy1(1))
          stat=atomicAdd(Jy(i1+1,j1,k1  ), Jy1(2))
#ifndef twoD
               stat=atomicAdd(Jy(i1  ,j1,k1+1), Jy1(3))
               stat=atomicAdd(Jy(i1+1,j1,k1+1), Jy1(4))
#endif



          stat=atomicAdd(Jz(i1  ,j1  ,k1), Jz1(1))
          stat=atomicAdd(Jz(i1+1,j1  ,k1), Jz1(2))
          stat=atomicAdd(Jz(i1  ,j1+1,k1), Jz1(3))
          stat=atomicAdd(Jz(i1+1,j1+1,k1), Jz1(4))

 	     end if
		  






!------------------------------------------------
! warp level deposition for the second  time
!-----------------------------------------------

#ifdef twoD
       cell_index=j2*(x2-x1+5)+i2
#else
       cell_index=k2*(y2-y1+5)*(x2-x1+5)+j2*(x2-x1+5)+i2
#endif

       lane=mod((threadIdx%x-1), 32) + 1 !WARP_SIZE is assumed to be 32
       unclaimed= B'11111111111111111111111111111111'
	   is_peer=.false.
	  do while(.not.is_peer)
          other_cell_index=__SHFL(cell_index,__FFS(unclaimed))
	      is_peer=cell_index.eq.other_cell_index
		  peers=ballot(is_peer)
		  unclaimed=XOR(unclaimed,peers)
      end do

	  first_lane=__FFS(peers)
	  rel_pos=__popc(IBITS(peers,0,lane))-1
      peers=IAND(peers,ISHFT(B'11111111111111111111111111111110',lane-1))
	  do while(anythread((peers.ne.0)))


		  next=__FFS(peers)

		  !Copy data from other threads in the warp
			  recvJ(1)=__SHFL(Jx2(1),next)
			  recvJ(2)=__SHFL(Jx2(2),next)
#ifndef twoD
			  recvJ(3)=__SHFL(Jx2(3),next)
			  recvJ(4)=__SHFL(Jx2(4),next)
#endif
			  if(next.gt.0) then
				  Jx2(1)=Jx2(1)+recvJ(1)
				  Jx2(2)=Jx2(2)+recvJ(2)
#ifndef twoD
				  Jx2(3)=Jx2(3)+recvJ(3)
				  Jx2(4)=Jx2(4)+recvJ(4)
#endif
			  end if

			  recvJ(1)=__SHFL(Jy2(1),next)
			  recvJ(2)=__SHFL(Jy2(2),next)
#ifndef twoD
			  recvJ(3)=__SHFL(Jy2(3),next)
			  recvJ(4)=__SHFL(Jy2(4),next)
#endif
			  if(next.gt.0) then
				  Jy2(1)=Jy2(1)+recvJ(1)
				  Jy2(2)=Jy2(2)+recvJ(2)
#ifndef twoD
				  Jy2(3)=Jy2(3)+recvJ(3)
				  Jy2(4)=Jy2(4)+recvJ(4)
#endif
			  end if

			  recvJ(1)=__SHFL(Jz2(1),next)
			  recvJ(2)=__SHFL(Jz2(2),next)
			  recvJ(3)=__SHFL(Jz2(3),next)
			  recvJ(4)=__SHFL(Jz2(4),next)
			  if(next.gt.0) then
				  Jz2(1)=Jz2(1)+recvJ(1)
				  Jz2(2)=Jz2(2)+recvJ(2)
				  Jz2(3)=Jz2(3)+recvJ(3)
				  Jz2(4)=Jz2(4)+recvJ(4)
			  end if

			  !----------------------------
			  done=IAND(rel_pos,1)
			  if(done.eq.1) then
				   not_done=.false.
			  else
				   not_done=.true.
			  end if
			  peers=IAND(peers,ballot(not_done))
			  !peers=IAND(peers,ballot(NOT(done)))
			  rel_pos=ISHFT(rel_pos,-1)
	  end do







	      if(lane.eq.first_lane) then
		  !   if(lane.eq.1) then


          stat=atomicAdd(Jx(i2,j2,  k2  ), Jx2(1) )
          stat=atomicAdd(Jx(i2,j2+1,k2  ), Jx2(2) )
#ifndef twoD
               stat= atomicAdd(Jx(i2,j2,  k2+1), Jx2(3))
               stat= atomicAdd(Jx(i2,j2+1,k2+1), Jx2(4))
#endif

          stat=atomicAdd(Jy(i2  ,j2,k2  ), Jy2(1))
          stat=atomicAdd( Jy(i2+1,j2,k2  ), Jy2(2))
#ifndef twoD
               stat=atomicAdd(Jy(i2  ,j2,k2+1), Jy2(3))
               stat=atomicAdd(Jy(i2+1,j2,k2+1), Jy2(4))
#endif


          stat=atomicAdd(Jz(i2  ,j2  ,k2), Jz2(1))
          stat=atomicAdd(Jz(i2+1,j2  ,k2), Jz2(2))
          stat=atomicAdd(Jz(i2  ,j2+1,k2), Jz2(3))
          stat=atomicAdd(Jz(i2+1,j2+1,k2), Jz2(4))

	     end if

 
 
 
 
              
! Periodic Boundary conditions 
#ifdef twoD
             if(zp(n).ge.2.0) then
               zp(n)=zp(n)-1.0
             else if(zp(n).lt.1.0) then
               zp(n)=1.0+zp(n)
             end if 
#endif   
	  

end subroutine MoveDepositPrtlKernel


#endif 

















 	attributes(global) subroutine ResetVecCurrentGPUKernel(Fld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD
 		real, dimension(12,x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: Fld
#else 
        real, dimension(8,x1-2:x2+2,y1-2:y2+2,1:1) :: Fld
#endif  		
 		integer :: i,j,k,nn
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif  
         if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) then 
#ifdef twoD			 
			 do nn=1,8
#else 
             do nn=1,12
#endif 				 
			     Fld(nn,i,j,k)=0.0
		     end do 
		 end if 
 	end subroutine ResetVecCurrentGPUKernel
	
	
	attributes(global) subroutine InitTexExGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-2
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif 		
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.5*(Fld(i-1,j,k)+Fld(i,j,k))
	end subroutine InitTexExGPUKernel
	
	attributes(global) subroutine InitTexEyGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-2
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif 		
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.5*(Fld(i,j-1,k)+Fld(i,j,k))
	end subroutine InitTexEyGPUKernel
	
	attributes(global) subroutine InitTexEzGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-2
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.5*(Fld(i,j,k-1)+Fld(i,j,k))
#else
        k=1
		if((i.le.x2+2).and.(j.le.y2+2)) TexFld(i,j,k)=Fld(i,j,k)
#endif 		
	end subroutine InitTexEzGPUKernel
	
	attributes(global) subroutine InitTexBxGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2,extB)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  		
        real, value :: extB 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-3
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-2
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-2
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.25*(Fld(i,j-1,k-1)+Fld(i,j,k-1)+Fld(i,j-1,k)+Fld(i,j,k)) + extB
#else
        k=1
		if((i.le.x2+2).and.(j.le.y2+2)) TexFld(i,j,k)=0.5*(Fld(i,j-1,k)+Fld(i,j,k)) +extB
#endif 		
	end subroutine InitTexBxGPUKernel
	attributes(global) subroutine InitTexByGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2,extB)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  
        real, value :: extB 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-2
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-3
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-2
		if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.25*(Fld(i-1,j,k-1)+Fld(i-1,j,k)+Fld(i,j,k-1)+Fld(i,j,k)) +extB
#else
        k=1
		if((i.le.x2+2).and.(j.le.y2+2)) TexFld(i,j,k)=0.5*(Fld(i-1,j,k)+Fld(i,j,k)) +extB
#endif 		
	end subroutine InitTexByGPUKernel
	attributes(global) subroutine InitTexBzGPUKernel(TexFld,Fld,x1,x2,y1,y2,z1,z2,extB)
 		integer :: x1,x2,y1,y2,z1,z2
#ifndef twoD 		
 		real, dimension(x1-2:x2+2,y1-2:y2+2,z1-2:z2+2) :: TexFld,Fld
#else 
        real, dimension(x1-2:x2+2,y1-2:y2+2,1) :: TexFld,Fld
#endif  
        real, value :: extB		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +x1-2
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +y1-2
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +z1-3
#else
        k=1
#endif 		
        if((i.le.x2+2).and.(j.le.y2+2).and.(k.le.z2+2)) TexFld(i,j,k)=0.25*(Fld(i-1,j-1,k)+Fld(i-1,j,k)+Fld(i,j-1,k)+Fld(i,j,k)) +extB
	end subroutine InitTexBzGPUKernel
	
end module movdep_gpu