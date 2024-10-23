module movdep
     use parameters
     use vars 
	 use cyl_common
	 use movdep_routines
     use deposit 
     use interpolation
	 
     implicit none

contains 

	 
     subroutine MoveDepositPrtl
		  logical :: axis_proc
          integer :: n,jc,kc,off,nn,off_max,thread_prtl_block_size
		  real(psn), dimension(VecBlockSize) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,ui,vi,wi,u1,v1,w1,f,g
          real(psn), dimension(VecBlockSize) :: x0,y0,z0
		  real(psn), dimension(VecBlockSize) :: x1,y1,r1
		  integer, dimension(VecBlockSize):: ip0,jp0,kp0
		  integer :: ip,jp,kp
		  real(psn)    :: xr,yr,zr
		  real(psn)    :: Fx1,Fx2,Fy1,Fy2,Fz1,Fz2,Wx1,Wx2,Wy1,Wy2,Wz1,Wz2
		  real(psn), dimension(VecBlockSize) :: qthis
		  real(psn) :: r_proc
          integer :: i,j,k
		  real(psn), dimension(VecBlockSize) :: qm
		  real(psn) :: dx,dy,dz
		  integer :: VecBlockSizeThis
		  integer :: joff,koff,short_arr_ind
#ifdef twoD 
          real(psn), dimension(4) :: FldVec,wt
#else
          real(psn), dimension(8) :: FldVec,wt 
#endif  


 		  call GatherVecEfld(VecEx,VecEy,VecEz,Ex,Ey,Ez)

 		  call GatherVecBfld

		  VecJ=0.0_psn ! reset the value of current to zero


!-------------------------------------------------------------------------------------		
! First Initialize All the fields needed in the mover 
!-------------------------------------------------------------------------------------
		  qm=0.0_psn 
		  pEx=0.0_psn 
		  pEy=0.0_psn 
		  pEz=0.0_psn
		  pBx=0.0_psn 
		  pBy=0.0_psn 
		  pBz=0.0_psn 
		  qthis=0.0_psn

		  r_proc=xborders(procxind)+rshift-xmin
		  axis_proc=.false.
		  if(inc_axis) then
			  if(procxind.eq.0) then 
				  axis_proc=.true. 		
			  end if 
		  end if 

!call UpdateFldAxis
 



!off=ThreadID*thread_prtl_block_size 
!off_max=min(off+thread_prtl_block_size,used_prtl_arr_size)
!-------------------------------------------------------------------------------------		
! Now Move the Particles 
!-------------------------------------------------------------------------------------
	
#ifndef OPEN_MP
    ThreadID=0	
#else
!$OMP PARALLEL DO PRIVATE(off,VecBlockSizeThis,n,x0,y0,z0,qm,pEx,pEy,pEz,pBx,pBy,pBz,ip0,jp0,kp0,dx,dy,dz,nn,wt,short_arr_ind,u0,v0,w0,ui,vi,wi,x1,y1,r1,g,f,u1,v1,w1,ip,jp,kp,xr,yr,zr,qthis,Fx1,Fy1,Fz1,Wx1,Wy1,Wz1,Fx2,Fy2,Fz2,Wx2,Wy2,Wz2,ThreadID) 
#endif 
	 do off=0,used_prtl_arr_size-1,VecBlockSize 
		 VecBlockSizeThis=min(VecBlockSize,used_prtl_arr_size-off)
#ifdef OPEN_MP		 
		 ThreadID=OMP_GET_THREAD_NUM()
#endif 		 
		!initial position of particles  
		do n=1,VecBlockSizeThis 
			  x0(n)=xp(n+off) !local r
	          y0(n)=yp(n+off) !local theta
	          z0(n)=zp(n+off)
	          qthis(n)=qp(n+off)*qi
		end do 
		if(axis_proc) call CopyInitMom(off,VecBlockSizeThis,ui,vi,wi)

		! load q/m for each particles into a vector
		do n=1,VecBlockSizeThis 
		    if(flvp(n+off).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0 
			  qm(n)=0
		    else 			
  		      qm(n)=flvrqm(flvp(n+off))*0.5_psn !for optimisation 0.5 is multiplied here itself 
		    end if 
		end do 
        !interpolation 
	    pEx=0.0_psn
	    pEy=0.0_psn
	    pEz=0.0_psn
	    pBx=0.0_psn
	    pBy=0.0_psn
	    pBz=0.0_psn				
		do n=1,VecBlockSizeThis 
		  ip0(n)=x0(n)
		  jp0(n)=y0(n)
		  dx=x0(n)-ip0(n)
		  dy=y0(n)-jp0(n) 
#ifndef twoD 
          kp0(n)=z0(n)
          dz=z0(n)-kp0(n) 	
!$OMP SIMD
		  do nn=1,8
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
		  end do

#else
          kp0(n)=1
!$OMP SIMD
		  do nn=1,4
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
		  end do
#endif


#ifdef twoD
          short_arr_ind= (jp0(n)-1)*mx + ip0(n)
#else
          short_arr_ind= (kp0(n)-1)*mx*my  + (jp0(n)-1)*mx + ip0(n)
#endif


!$OMP SIMD
#ifdef twoD
		  do nn=1,4
#else
          do nn=1,8
#endif
			  pEx(n)=pEx(n)+wt(nn)*VecEx(nn,short_arr_ind)
			  pEy(n)=pEy(n)+wt(nn)*VecEy(nn,short_arr_ind)
			  pEz(n)=pEz(n)+wt(nn)*VecEz(nn,short_arr_ind)
			  pBx(n)=pBx(n)+wt(nn)*VecBx(nn,short_arr_ind)
			  pBy(n)=pBy(n)+wt(nn)*VecBy(nn,short_arr_ind)
			  pBz(n)=pBz(n)+wt(nn)*VecBz(nn,short_arr_ind)

		  end do
           !call InterpEMFldGridPointYeeMesh(xp(n+off),yp(n+off),zp(n+off),pEx(n),pEy(n),pEz(n),pBx(n),pBy(n),pBz(n),Ex,Ey,Ez,Bx,By,Bz)  
		
		end do
		
		
            pEx=pEx*qm
            pEy=pEy*qm
            pEz=pEz*qm
            pBx=pBx*qm
            pBy=pBy*qm
            pBz=pBz*qm      

!$OMP SIMD
       do n=1,VecBlockSizeThis    
              !Boris Pusher 
               u0(n)=c*up(n+off)+pEx(n)
               v0(n)=c*vp(n+off)+pEy(n)
               w0(n)=c*wp(n+off)+pEz(n)

               g(n)=1.0_psn/sqrt(sqc+u0(n)**2+v0(n)**2+w0(n)**2)   ! 1/c*gamma
               pBx(n)=g(n)*pBx(n)
               pBy(n)=g(n)*pBy(n)
               pBz(n)=g(n)*pBz(n)

               f(n)=2.0_psn/(1.0_psn+pBx(n)*pBx(n)+pBy(n)*pBy(n)+pBz(n)*pBz(n))
               u1(n)=(u0(n)+v0(n)*pBz(n)-w0(n)*pBy(n))*f(n)
               v1(n)=(v0(n)+w0(n)*pBx(n)-u0(n)*pBz(n))*f(n)
               w1(n)=(w0(n)+u0(n)*pBy(n)-v0(n)*pBx(n))*f(n)

               u0(n)=u0(n)+v1(n)*pBz(n)-w1(n)*pBy(n)+pEx(n)
               v0(n)=v0(n)+w1(n)*pBx(n)-u1(n)*pBz(n)+pEy(n)
               w0(n)=w0(n)+u1(n)*pBy(n)-v1(n)*pBx(n)+pEz(n)

               u1(n)=u0(n)*cinv
               v1(n)=v0(n)*cinv
               wp(n+off)=w0(n)*cinv

               g(n)=sqc/sqrt(sqc+u0(n)**2+v0(n)**2+w0(n)**2)

               x1(n)=(xp(n+off) + r_proc)*grid_dx + u1(n)*g(n)
			   y1(n)=v1(n)*g(n)
			   r1(n)=sqrt(x1(n)**2+y1(n)**2)
			   xp(n+off)=r1(n)*grid_inv_dx - r_proc
               yp(n+off)=yp(n+off) + atan(y1(n)/x1(n))*inv_dtheta
               zp(n+off)=zp(n+off) + wp(n+off)*g(n)*grid_inv_dz
			   
			   r1(n)=1.0_psn/r1(n)
			   up(n+off)=(x1(n)*u1(n)+y1(n)*v1(n))*r1(n)
	           vp(n+off)=(-y1(n)*u1(n)+x1(n)*v1(n))*r1(n)
        end do     

		if(axis_proc) call MovDepAxis(off,VecBlockSizeThis,qthis,x0,y0,z0,ui,vi,wi)

        do n=1,VecBlockSizeThis

 	          ip=xp(n+off)
 	          jp=yp(n+off)
 	          kp=zp(n+off)
#ifdef twoD
              kp=1
#endif
                xr=min(real(min(ip0(n),ip)+1,psn),max(real(max(ip0(n),ip),psn),0.5_psn*(x0(n)+xp(n+off))))
                yr=min(real(min(jp0(n),jp)+1,psn),max(real(max(jp0(n),jp),psn),0.5_psn*(y0(n)+yp(n+off))))
                zr=min(real(min(kp0(n),kp)+1,psn),max(real(max(kp0(n),kp),psn),0.5_psn*(z0(n)+zp(n+off))))

 				Fx1=qthis(n)*(xr-x0(n))
 		        Fy1=qthis(n)*(yr-y0(n))
				Fz1=qthis(n)*(zr-z0(n))
 		        Wx1=0.5_psn*(x0(n)+xr)-ip0(n)
 		        Wy1=0.5_psn*(y0(n)+yr)-jp0(n)
#ifdef twoD
		        Wz1=0.0_psn
				short_arr_ind= (jp0(n)-1)*mx + ip0(n)
#else
		        Wz1=0.5_psn*(z0(n)+zr)-kp0(n)
				short_arr_ind= (kp0(n)-1)*mx*my  + (jp0(n)-1)*mx + ip0(n)
#endif



!$OMP SIMD
#ifdef twoD
                do nn=1,8
#else
                do nn=1,12
#endif
                       VecJ(nn,short_arr_ind,ThreadID+1)=VecJ(nn,short_arr_ind,ThreadID+1)+(wtFx(nn)*Fx1+wtFy(nn)*Fy1+wtFz(nn)*Fz1)*(Jwtx1(nn)+Jwtx2(nn)*Wx1)*(Jwty1(nn)+Jwty2(nn)*Wy1)*(Jwtz1(nn)+Jwtz2(nn)*Wz1)
                end do

		        Fx2=qthis(n)*(xp(n+off)-xr)
		        Fy2=qthis(n)*(yp(n+off)-yr)
		        Fz2=qthis(n)*(zp(n+off)-zr)
		        Wx2=0.5_psn*(xp(n+off)+xr)-ip
		        Wy2=0.5_psn*(yp(n+off)+yr)-jp
#ifdef twoD
		        Wz2=0.0_psn
				short_arr_ind= (jp0(n)-1)*mx + ip0(n)
#else
		        Wz2=0.5_psn*(zp(n+off)+zr)-kp
				short_arr_ind= (kp0(n)-1)*mx*my  + (jp0(n)-1)*mx + ip0(n)
#endif


!$OMP SIMD
#ifdef twoD
                do nn=1,8
#else
                do nn=1,12
#endif

					   VecJ(nn,short_arr_ind,ThreadID+1)=VecJ(nn,short_arr_ind,ThreadID+1)+(wtFx(nn)*Fx2+wtFy(nn)*Fy2+wtFz(nn)*Fz2)*(Jwtx1(nn)+Jwtx2(nn)*Wx2)*(Jwty1(nn)+Jwty2(nn)*Wy2)*(Jwtz1(nn)+Jwtz2(nn)*Wz2)
				end do


  		end do


! Periodic Boundary conditions
#ifdef twoD
        do n=1,VecBlockSizeThis
             if(zp(n+off).gt.zmax) then
               zp(n+off)=zp(n+off)-zlen
             else if(zp(n+off).lt.zmin) then
               zp(n+off)=zlen+zp(n+off)
             end if
		end do
#endif

    end do ! end of the main do loop


     call ReduceCurrentMatrix

	 call AxisCurrentBC   
	 
	 
	 end subroutine MoveDepositPrtl
	 
	 subroutine CopyInitMom(off,pcount,u,v,w) 
		integer :: n,off,pcount
		real(psn), dimension(VecBlockSize) :: u,v,w
		do n=1,pcount 
 			  u(n)=up(n+off)
 	          v(n)=vp(n+off) 
 	          w(n)=wp(n+off)
 	    end do 
	 end subroutine CopyInitMom 
	 
	 subroutine MovDepAxis(off,pcount,q,x,y,z,u,v,w)
		 integer :: off, pcount
		 real(psn), dimension(VecBlockSize) ::q,x,y,z,u,v,w
		 integer :: n
		 
		 do n=1,pcount !first rectify the angular cordinate 
			 if(x(n).lt.4.5_psn) then ! if the particle is more than 1 unit away from the axis
				 call SubStepMovDep(n+off,x(n),y(n),z(n),u(n),v(n),w(n))
				 q(n)=0.0_psn
			 end if  
		 end do
		 
	 end subroutine MovDepAxis
	 
	 subroutine SubStepMovDep(n,x0,y0,z0,u0,v0,w0)
		 integer   :: n
		 real(psn) :: x0,y0,z0,u0,v0,w0
		 real(psn) :: t,dt,delt
		 integer   :: res
		 
		 if(flvp(n).eq.0) return
		 call ResetPosMom(n,x0,y0,z0,u0,v0,w0)
		 t=0.0_psn
		 dt=1.0_psn
		 do while(t.lt.1.0) ! advance along the straight line 
			 !determine dt
             res=0 
			 !print*,'n ',n,'t ',t,' dt ',dt,'r1',r1,'r2',r2
			 !if(r1.lt.0.5) print*,'dt',dt,'r1',r1,'r2',r2,'V0',v0,'u0',u0
			 delt=dt
			 if(t+dt.gt.1.0_psn) delt=1.0_psn-t
			 !print*,'n',n,'dt',dt
			 call MovePrtlAxis(n,x0,y0,z0,u0,v0,w0,delt,res)
			 if(res.eq.1) then 
				 t=t+delt
				 call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
				 call PrtlAxisBC(xp(n),yp(n),zp(n),up(n),vp(n),wp(n))
				 dt=dt*1.1_psn
				 x0=xp(n)
				 y0=yp(n)
				 z0=zp(n)
			 end if 
			 if(res.eq.0) then 
				 dt=dt*0.5_psn
			 end if 


		 end do 
	 end subroutine SubStepMovDep
	 subroutine ResetPosMom(n,x,y,z,u,v,w)
		 integer :: n
		 real(psn) :: x,y,z,u,v,w
		 xp(n)=x
		 yp(n)=y
		 zp(n)=z
		 up(n)=u
		 vp(n)=v
		 wp(n)=w
	 end subroutine ResetPosMom
	 
	 subroutine PrtlAxisBC(x,y,z,u,v,w)
		  real(psn) :: x,y,z,u,v,w
		  if(x.lt.3.5_psn) then
			  x= 3.5_psn + (3.5_psn-x)
			  y= y + ny/2.0_psn
			  u = - u
			  v = - v   
		  end if 
		  if(y.gt.ymax) y=y-ylen
		  if(y.lt.ymin) y=y+ylen
	 end subroutine PrtlAxisBC
	 
	 subroutine MovePrtlAxis(n,x,y,z,u,v,w,dt,res)
		 integer   :: n, res
		 real(psn) :: x,y,z,u,v,w,dt
		 real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz
		 real(psn) :: qm
		 real(psn) :: u0,v0,w0,u1,v1,w1,f,g,x1,y1,r1,dy
		 
		 call Interp_VecEM_GridPoints(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
		 !call InterpEMFldGridPointYeeMesh(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)  
		 
		 if(flvp(n).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0 
			  qm=0
		 else 			
  		      qm=flvrqm(flvp(n))*0.5_psn !for optimisation 0.5 is multiplied here
		 end if 
         pEx=pEx*qm*dt
         pEy=pEy*qm*dt
         pEz=pEz*qm*dt
         pBx=pBx*qm*dt
         pBy=pBy*qm*dt
         pBz=pBz*qm*dt      
		 
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

         u1=u0*cinv
         v1=v0*cinv

         g=sqc/sqrt(sqc+u0**2+v0**2+w0**2)

         x1= (x - 3.5_psn )*grid_dx + u1*g*dt 
	     y1=v1*g*dt
	     r1=sqrt(x1**2+y1**2)
		 
		 if(x1.gt.0) then 
		     dy= atan(y1/x1)*inv_dtheta
		 else 
			 dy= 0.0_psn
		 end if 
		 !print*,'n',n,'dt',dt,'dy',dy,'u1',u1,'v1',v1
		 if(abs(dy).gt.0.45_psn) return 
		 		
		 if(x1.gt.0) then 	 
	         xp(n)=r1*grid_inv_dx + 3.5_psn
	         up(n)=(x1*u1+y1*v1)/r1
             vp(n)=(-y1*u1+x1*v1)/r1	
		 else
			 xp(n)=x1*grid_inv_dx + 3.5_psn
	         up(n)=u1
			 vp(n)=v1
		 end if  
         yp(n)=yp(n)+dy
		 
		 wp(n)=w0*cinv
		 zp(n)=zp(n) + wp(n)*g*dt *grid_inv_dz
		 	   
		 res=1 
		 		 
	 end subroutine MovePrtlAxis
	 
	 
	 !Current "boundary condition" at the axis. current at negative radius is substracted because particle's effective charge ->0 as r->0
     subroutine AxisCurrentBC
 		integer :: j,k,jp,jm,jpp
		
 		if(.not.inc_axis) return
 		if(procxind.ne.0) return

           
 	    Jy(4,:,:)=Jy(4,:,:)-Jy(3,:,:)
        Jz(4,:,:)=Jz(4,:,:)-Jz(3,:,:)

 		Jx(3,:,:) = 0.0_psn
 		Jy(2,:,:) = 0.0_psn
 		Jz(2,:,:) = 0.0_psn

 	end subroutine AxisCurrentBC 
	 
	 

	     
     
end module movdep