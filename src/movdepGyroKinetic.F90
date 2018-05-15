module movdep
     use parameters
     use vars 
     use interpolation 
     use deposit 
     implicit none
     real  :: OmegaCdtMAX=pi/6
contains 
     subroutine MoveParticles
          integer :: n,i,j,k,l
          real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,u1,v1,w1,f,f1,f2,g,h,Bmag,rho
          real(psn) :: x0,y0,z0
          real(psn) :: xg,yg,zg !guiding center 
          real(psn) :: Exa,Eya,Eza !electric field along the magnetic field 
           real(psn) :: EXB_x,EXB_y,EXB_z !Drift velocity 
          real(psn) :: ua,va,wa,un,vn,wn
     do n=1,prtl_arr_size

        if(qp(n).eq.0) cycle

          x0=xp(n)
        y0=yp(n)
          z0=zp(n)
          !Choose the interpolation shceme (see interpolation.F90)
          !call InterpFields2D(x0,y0,z0,qp(n),pEx,pEy,pEz,pBx,pBy,pBz)
          call InterpFields_GridPoints_2D(x0,y0,z0,qp(n),pEx,pEy,pEz,pBx,pBy,pBz)
          
          
          !print*,'pEx is',pEx,qi     
          Bmag=sqrt(pBx*pBx+pBy*pBy+pBz*pBz)
          g=sqrt(1+up(n)**2+vp(n)**2+wp(n)**2) 
          
         !if(Bmag/(c*g).gt.OmegaCdtMAX) then !try to get rid of factor c here, for optimisation 
          if(.false.) then 
!--------------------------------------------------------------------
!MAGNETIC ROTATION+DRIFT : USED FOR MAGNETISED PARTICLES 
!--------------------------------------------------------------------
            !Normalise the magnetic field compoents 
               pBx=pBx/Bmag
               pBy=pBy/Bmag
               pBz=pBz/Bmag 
               !Drift : EXB/c
               EXB_x=(pEy*pBz-pEz*pBy)/Bmag
               EXB_y=(pEz*pBx-pEx*pBz)/Bmag
               EXB_z=(pEx*pBy-pEy*pBx)/Bmag
               !print*,'Drift EXB  is',EXB_x,EXB_y,EXB_z
               !Lorent factor of the particle
               g=sqrt(1+up(n)**2+vp(n)**2+wp(n)**2) 
               !print*,'particle gamma is ',g
               
               !get the component of gamma*v/c along the magnetic field 
               h=up(n)*pBx+vp(n)*pBy+wp(n)*pBz
               ua=pBx*h
               va=pBy*h 
               wa=pBz*h
               !print*,'parallel gv/c is',ua,va,wa
            !parallel component of the electric field 
               h=2*(pEx*pBx+pEy*pBy+pEz*pBz)
               Exa=pBx*h
               Eya=pBy*h
               Eza=pBz*h
               !print*,'parallel compoent of electric field is',Exa,Eya,Eza
               !get the perpendicular compoent of v/c normal to the mangetic field 
               un=(up(n)-ua)/g
               vn=(vp(n)-va)/g
               wn=(wp(n)-wa)/g
               !print*,'normal v/c is',un,vn,wn
               !change the perpenduclar component of v/c the gyrating part only = v_perp/c
               !h=(1-un*EXB-vn*EXB_y-wn*EXB_z)
               un=(un-EXB_x)!/h
               vn=(vn-EXB_y)!/h
               wn=(wn-EXB_z)!/h
               !compute the larmor radius 
               h=sqrt(un*un+vn*vn+wn*wn)
               rho=(g*sqc)/(2*Bmag) !rho is Gyroradius/sqrt(un*un+vn*vn+wn*wn)
               !print*,'rho is',rho,h
               !first find the current location of Guiding center
               xg=xp(n)+rho*(vn*pBz-wn*pBy)
               yg=yp(n)+rho*(wn*pBx-un*pBz)
               zg=zp(n)+rho*(un*pBy-vn*pBx)
               !print*,'guiding center is at',xg,yg,zg

            !acceleration due to the electric field   
               ua=ua+Exa*cinv
               va=va+Eya*cinv
               wa=wa+Eza*cinv
               !rotate the normal compoent of velicity vector 
               u0=un
               v0=vn
               w0=wn
               u1=(vn*pBz-wn*pBy)
               v1=(wn*pBx-un*pBz)
               w1=(un*pBy-vn*pBx)
               
               !print*,'angle of rotation',qp(n), 2*Bmag/g 
               f=2*Bmag/(c*g) !angle of rotation: omega*dt
               f1=cos(f)  
               f2=sin(f) 
               un=(u0*f1+u1*f2)
               vn=(v0*f1+v1*f2)
               wn=(w0*f1+w1*f2)
               !add the drift term to the velocity 
               un=un+EXB_x
               vn=vn+EXB_y
               wn=wn+EXB_z
               !update the four velocity of particle  
               up(n)=ua+un*g
               vp(n)=va+vn*g
               wp(n)=wa+wn*g
               
               !rotate the position vector w.r.t. guding center 
               xp(n)=xg-rho*(-u0*f2+u1*f1)
               yp(n)=yg-rho*(-v0*f2+v1*f1)
               zp(n)=zg-rho*(-w0*f2+w1*f1)
               
               !add the drift term 
               xp(n)=xp(n) + EXB_x*c
               yp(n)=yp(n) + EXB_y*c 
               zp(n)=zp(n) + EXB_z*c
               !move particle along the magentic field 
               xp(n)=xp(n) + (ua*c)/g
               yp(n)=yp(n) + (va*c)/g 
               zp(n)=zp(n) + (wa*c)/g

!--------------------------------------------------------------------
!END OF MAGENTIC ROTATION+DRIFT
!-------------------------------------------------------------------- 
        else 
!---------------------------------------------------------------------
! BORISH PUSHER : USED FOR PARTICLES WITH LARGE LARMOR RADIUS 
!---------------------------------------------------------------------
               u0=c*up(n)+pEx
               v0=c*vp(n)+pEy
               w0=c*wp(n)+pEz
               !First half magnetic rotation, with Lorentz gamma
               g=1/sqrt(sqc+u0**2+v0**2+w0**2)   ! 1/c*gamma 
               !print*,'Gamma is:',qp(n), 1/(g*c)
               
               pBx=g*pBx
               pBy=g*pBy
               pBz=g*pBz
               
               f=2.0/(1.0+pBx*pBx+pBy*pBy+pBz*pBz)
               u1=(u0+v0*pBz-w0*pBy)*f
               v1=(v0+w0*pBx-u0*pBz)*f
               w1=(w0+u0*pBy-v0*pBx)*f
               ! another half magnatic rotation and elelctric acceleration
               u0=u0+v1*pBz-w1*pBy+pEx
               v0=v0+w1*pBx-u1*pBz+pEy 
               w0=w0+u1*pBy-v1*pBx+pEz 
               !4-velocity,Normalised
               up(n)=u0*cinv
               vp(n)=v0*cinv
               wp(n)=w0*cinv
          
              !move particle 
               g=c/sqrt(c**2+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g*c
               yp(n)=yp(n) + vp(n)*g*c 
               zp(n)=zp(n) + wp(n)*g*c
               
               !print*,'Change in Y is:',vp(n)*g*c,'and X is:',up(n)*g*c
               
!---------------------------------------------------------------------------
!   END OF BORISH PUSHER
!---------------------------------------------------------------------------               
       end if


#ifdef mulflvr
           if(p(n)%tp.ne.1) call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n)) !it is a test particle  
#else 
           call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n)) 
#endif
             
! Periodic Boundary conditions 

             if(zp(n).gt.zmax) then
               zp(n)=zp(n)-zlen
             else if(zp(n).lt.zmin) then
               zp(n)=zlen+zp(n)
           end if      

    end do 
                    
     end subroutine MoveParticles 
     
!      subroutine InterpFields2D(x,y,z,q,pEx,pEy,pEz,pBx,pBy,pBz)
!           !Note: this subroutine returns field*0.5, the 0.5 factor is multiplied here for optimisation purpose
!           real(psn) :: qthis,x,y,z,pEx,pEy,pEz,pBx,pBy,pBz
!           integer :: i,j
!           real(psn):: dx,dy
!           real :: q
!
! !                 if(xp(n).lt.2) then
! !                      print*,'caught here',xp(n),tagp(n),up(n),proc
! !                 end if
!                i=anint(x)
!                j=y
!                dx=x-i+0.5
!                dy=y-j
!                call Interp2D(i-1,j,Ex,dx,dy,pEx)
!
!                i=x
!                j=anint(y)
!                dx=x-i
!                dy=y-j+0.5
!                call Interp2D(i,j-1,Ey,dx,dy,pEy)
!
!                i=x
!                j=y
!                dx=x-i
!                dy=y-j
!                call Interp2D(i,j,Ez,dx,dy,pEz)
!
! ! now interpolate the magnetic field
!                i=x
!                j=aint(y-0.5)
!                dx=x-i
!                dy=y-j-0.5
!                call Interp2D(i,j,Bx,dx,dy,pBx)
!
!                i=aint(x-0.5)
!                j=y
!                dx=x-i-0.5
!                dy=y-j
!                call Interp2D(i,j,By,dx,dy,pBy)
!
!                i=aint(x-0.5)
!                j=aint(y-0.5)
!                dx=x-i-0.5
!                dy=y-j-0.5
!                call Interp2D(i,j,Bz,dx,dy,pBz)
!
! #ifndef mulflvr
!                if(q.gt.0) then
!                     qthis=qmi*0.5
!                end if
!                if(q.lt.0) then
!                     qthis=qme*0.5
!                end if
! #else
! !!!!! This is a temporary fix
!                if(q.gt.0) then
!                     qthis=qmi*0.5
!                end if
!                if(q.lt.0) then
!                     qthis=abs(q)*qme*0.5
!                end if
!
! #endif
!                     pEx=pEx*qthis
!                     pEy=pEy*qthis
!                     pEz=pEz*qthis
! #ifndef Bext0
!                     pBx=pBx*qthis
!                     pBy=pBy*qthis
!                     pBz=pBz*qthis
! #else
!                     pBx=(pBx+Bx_ext0)*qthis
!                     pBy=(pBy+By_ext0)*qthis
!                     pBz=(pBz+Bz_ext0)*qthis
! #endif
!
!      end subroutine InterpFields2D
!      subroutine Interp3D(i,j,k,F,dx,dy,dz,res)
!           integer :: i,j,k
!           real(psn) :: dx,dy,dz,res
!           real(psn),dimension(mx,my,mz) :: F
!           res=(1-dx)*(1-dy)*(1-dz)*F(i,  j  ,k)+&
!               dx    *(1-dy)*(1-dz)*F(i+1,j  ,k)+&
!                (1-dx)*dy    *(1-dz)*F(i  ,j+1,k)+&
!                (1-dx)*(1-dy)*dz    *F(i,  j,  k+1)+&
!                (1-dx)*dy    *dz    *F(i,  j+1,k+1)+&
!                dx    *(1-dy)*dz    *F(i+1,j,  k+1)+&
!                dx    *dy    *(1-dz)*F(i+1,j+1,k)+&
!                dx    *dy    *dz    *F(i+1,j+1,k+1)
!      end subroutine Interp3D
!      subroutine Interp2D(i,j,F,dx,dy,res)
!           integer :: i,j
!           real(psn) :: dx,dy,res
!           real(psn),dimension(mx,my,1) :: F
!           res=(1-dx)*(1-dy)*F(i,  j  ,1)+&
!               dx    *(1-dy)*F(i+1,j  ,1)+&
!                (1-dx)*dy    *F(i  ,j+1,1)+&
!                dx    *dy    *F(i+1,j+1,1)
! !                if(isnan(res)) then
! !                     print*,i,j,dx,dy,F(i,j,1),F(i+1,j+1,1),F(i+1,j,1),F(i,j+1,1),mx,my
! !                end if
!      end subroutine Interp2D
     
!      subroutine DepositCurrentPIC(x0,y0,z0,x,y,z,q)
!           implicit none
!
!           ! local variables
!           real(psn) ::q
!           real(xpsn)::x
!           real(ypsn)::y
!           real(zpsn)::z
!           real(psn) :: xr,yr,zr,x0,y0,z0
!           real(psn) ::qthis
!
!           real(psn) ::Fx1, Fx2, Fy1, Fy2, Fz1, Fz2
!           real(psn) ::Wx1, Wx2, Wy1, Wy2, Wz1, Wz2
!           integer :: i1, i2, j1, j2, k1, k2
!
!         qthis=q*qi   ! qi=-qe
!
!           i1=aint(x0)
!           i2=aint(x)
!           j1=aint(y0)
!           j2=aint(y)
!           k1=aint(z0)
!           k2=aint(z)
! #ifdef twoD
!         k1=1
!         k2=1
! #endif
!                xr=min(real(min(i1,i2)+1),max(real(max(i1,i2)),.5*(x0+x)))
!                yr=min(real(min(j1,j2)+1),max(real(max(j1,j2)),.5*(y0+y)))
!                zr=min(real(min(k1,k2)+1),max(real(max(k1,k2)),.5*(z0+z)))
!           Fx1=qthis*(xr-x0)
!           Fy1=qthis*(yr-y0)
!           Fz1=qthis*(zr-z0)
!
!           Wx1=.5*(x0+xr)-i1
!           Wy1=.5*(y0+yr)-j1
!
!           Wx2=.5*(x+xr)-i2
!           Wy2=.5*(y+yr)-j2
!
! #ifdef twoD
!           Wz1=0
!           Wz2=0
! #else
!           Wz1=.5*(z0+zr)-k1
!           Wz2=.5*(z+zr)-k2
! #endif
!
!
!           Fx2=qthis*(x-xr)
!           Fy2=qthis*(y-yr)
!           Fz2=qthis*(z-zr)
!
!           Jx(i1,j1,  k1  )= Jx(i1,j1,  k1  )+Fx1 * (1.-Wy1)*(1.-Wz1)
!           Jx(i1,j1+1,k1  )= Jx(i1,j1+1,k1  )+Fx1 *  Wy1    *(1.-Wz1)
! #ifndef twoD
!                Jx(i1,j1,  k1+1)= Jx(i1,j1,  k1+1)+Fx1 * (1-Wy1) * Wz1
!                Jx(i1,j1+1,k1+1)= Jx(i1,j1+1,k1+1)+Fx1 *  Wy1    * Wz1
! #endif
!
!           Jx(i2,j2,  k2  )= Jx(i2,j2,  k2  )+Fx2 * (1.-Wy2)*(1.-Wz2)
!           Jx(i2,j2+1,k2  )= Jx(i2,j2+1,k2  )+Fx2 *  Wy2    *(1.-Wz2)
! #ifndef twoD
!                Jx(i2,j2,  k2+1)= Jx(i2,j2,  k2+1)+Fx2 * (1.-Wy2)* Wz2
!                Jx(i2,j2+1,k2+1)= Jx(i2,j2+1,k2+1)+Fx2 *  Wy2    * Wz2
! #endif
!
!
!           Jy(i1  ,j1,k1  )= Jy(i1  ,j1,k1  )+Fy1 * (1.-Wx1)*(1.-Wz1)
!           Jy(i1+1,j1,k1  )= Jy(i1+1,j1,k1  )+Fy1 *  Wx1    *(1.-Wz1)
! #ifndef twoD
!                Jy(i1  ,j1,k1+1)= Jy(i1  ,j1,k1+1)+Fy1 * (1.-Wx1)* Wz1
!                Jy(i1+1,j1,k1+1)= Jy(i1+1,j1,k1+1)+Fy1 *  Wx1    * Wz1
! #endif
!           Jy(i2  ,j2,k2  )= Jy(i2  ,j2,k2  )+Fy2 * (1.-Wx2)*(1.-Wz2)
!           Jy(i2+1,j2,k2  )= Jy(i2+1,j2,k2  )+Fy2 *  Wx2    *(1.-Wz2)
! #ifndef twoD
!                Jy(i2  ,j2,k2+1)= Jy(i2  ,j2,k2+1)+Fy2 * (1.-Wx2)* Wz2
!                Jy(i2+1,j2,k2+1)= Jy(i2+1,j2,k2+1)+Fy2 *  Wx2    * Wz2
! #endif
!
!
!           Jz(i1  ,j1  ,k1)= Jz(i1  ,j1  ,k1)+Fz1 * (1.-Wx1)*(1.-Wy1)
!           Jz(i1+1,j1  ,k1)= Jz(i1+1,j1  ,k1)+Fz1 *  Wx1    *(1.-Wy1)
!           Jz(i1  ,j1+1,k1)= Jz(i1  ,j1+1,k1)+Fz1 * (1.-Wx1)* Wy1
!           Jz(i1+1,j1+1,k1)= Jz(i1+1,j1+1,k1)+Fz1 *  Wx1    * Wy1
!
!           Jz(i2  ,j2  ,k2)= Jz(i2  ,j2  ,k2)+Fz2 * (1.-Wx2)*(1.-Wy2)
!           Jz(i2+1,j2  ,k2)= Jz(i2+1,j2  ,k2)+Fz2 *  Wx2    *(1.-Wy2)
!           Jz(i2  ,j2+1,k2)= Jz(i2  ,j2+1,k2)+Fz2 * (1.-Wx2)* Wy2
!           Jz(i2+1,j2+1,k2)= Jz(i2+1,j2+1,k2)+Fz2 *  Wx2    * Wy2
!
!      end subroutine DepositCurrentPIC
     
end module movdep