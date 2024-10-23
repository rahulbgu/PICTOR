module deposit 
     use parameters
     use vars
     implicit none 
contains 
     
!----------------------------------------------------------------------------------------------     
!DepositCurrentPIC : Charge conserving deposition method, used in the normal PIC simualtion 
!----------------------------------------------------------------------------------------------
     
     subroutine DepositCurrentPIC(x0,y0,z0,x,y,z,q)
          implicit none

          ! local variables
          real(psn), intent(IN)::x0,y0,z0,x,y,z,q
          real(psn) :: xr,yr,zr
          real(psn) ::qthis 
          
          real(psn) ::Fx1, Fx2, Fy1, Fy2, Fz1, Fz2
          real(psn) ::Wx1, Wx2, Wy1, Wy2, Wz1, Wz2
          integer :: i1, i2, j1, j2, k1, k2
               
          qthis=q*qi   ! q = particle's weight X sign of the charge
     
          i1=floor(x0)
          i2=floor(x)
          j1=floor(y0)
          j2=floor(y)
          k1=floor(z0)
          k2=floor(z)
! #ifdef twoD
!         k1=1
!         k2=1
! #endif
          xr=min(real(min(i1,i2)+1,psn),max(real(max(i1,i2),psn),0.5_psn*(x0+x)))
          yr=min(real(min(j1,j2)+1,psn),max(real(max(j1,j2),psn),0.5_psn*(y0+y)))
          zr=min(real(min(k1,k2)+1,psn),max(real(max(k1,k2),psn),0.5_psn*(z0+z)))
			   
          Fx1=qthis*(xr-x0)
          Fy1=qthis*(yr-y0)
          Fz1=qthis*(zr-z0)
     
          Wx1=0.5_psn*(x0+xr)-i1
          Wy1=0.5_psn*(y0+yr)-j1

          Wx2=0.5_psn*(x+xr)-i2
          Wy2=0.5_psn*(y+yr)-j2
          
#ifdef twoD
		  k1=1
		  k2=1
          Wz1=0.0_psn
          Wz2=0.0_psn
#else
          Wz1=0.5_psn*(z0+zr)-k1
          Wz2=0.5_psn*(z+zr)-k2
#endif

     
          Fx2=qthis*(x-xr)
          Fy2=qthis*(y-yr)
          Fz2=qthis*(z-zr)
		  
          Jx(i1,j1,  k1  )= Jx(i1,j1,  k1  )+Fx1 * (1.0_psn-Wy1)*(1.0_psn-Wz1)
          Jx(i1,j1+1,k1  )= Jx(i1,j1+1,k1  )+Fx1 *  Wy1    *(1.0_psn-Wz1)
#ifndef twoD
               Jx(i1,j1,  k1+1)= Jx(i1,j1,  k1+1)+Fx1 * (1.0_psn-Wy1) * Wz1
               Jx(i1,j1+1,k1+1)= Jx(i1,j1+1,k1+1)+Fx1 *  Wy1    * Wz1
#endif
     
          Jx(i2,j2,  k2  )= Jx(i2,j2,  k2  )+Fx2 * (1.0_psn-Wy2)*(1.0_psn-Wz2)
          Jx(i2,j2+1,k2  )= Jx(i2,j2+1,k2  )+Fx2 *  Wy2    *(1.0_psn-Wz2)
#ifndef twoD
               Jx(i2,j2,  k2+1)= Jx(i2,j2,  k2+1)+Fx2 * (1.0_psn-Wy2)* Wz2
               Jx(i2,j2+1,k2+1)= Jx(i2,j2+1,k2+1)+Fx2 *  Wy2    * Wz2
#endif
     

          Jy(i1  ,j1,k1  )= Jy(i1  ,j1,k1  )+Fy1 * (1.0_psn-Wx1)*(1.0_psn-Wz1)
          Jy(i1+1,j1,k1  )= Jy(i1+1,j1,k1  )+Fy1 *  Wx1    *(1.0_psn-Wz1) 
#ifndef twoD
               Jy(i1  ,j1,k1+1)= Jy(i1  ,j1,k1+1)+Fy1 * (1.0_psn-Wx1)* Wz1
               Jy(i1+1,j1,k1+1)= Jy(i1+1,j1,k1+1)+Fy1 *  Wx1    * Wz1
#endif
          Jy(i2  ,j2,k2  )= Jy(i2  ,j2,k2  )+Fy2 * (1.0_psn-Wx2)*(1.0_psn-Wz2)
          Jy(i2+1,j2,k2  )= Jy(i2+1,j2,k2  )+Fy2 *  Wx2    *(1.0_psn-Wz2) 
#ifndef twoD
               Jy(i2  ,j2,k2+1)= Jy(i2  ,j2,k2+1)+Fy2 * (1.0_psn-Wx2)* Wz2
               Jy(i2+1,j2,k2+1)= Jy(i2+1,j2,k2+1)+Fy2 *  Wx2    * Wz2
#endif
          

          Jz(i1  ,j1  ,k1)= Jz(i1  ,j1  ,k1)+Fz1 * (1.0_psn-Wx1)*(1.0_psn-Wy1)
          Jz(i1+1,j1  ,k1)= Jz(i1+1,j1  ,k1)+Fz1 *  Wx1    *(1.0_psn-Wy1)
          Jz(i1  ,j1+1,k1)= Jz(i1  ,j1+1,k1)+Fz1 * (1.0_psn-Wx1)* Wy1
          Jz(i1+1,j1+1,k1)= Jz(i1+1,j1+1,k1)+Fz1 *  Wx1    * Wy1
     
          Jz(i2  ,j2  ,k2)= Jz(i2  ,j2  ,k2)+Fz2 * (1.0_psn-Wx2)*(1.0_psn-Wy2)
          Jz(i2+1,j2  ,k2)= Jz(i2+1,j2  ,k2)+Fz2 *  Wx2    *(1.0_psn-Wy2)
          Jz(i2  ,j2+1,k2)= Jz(i2  ,j2+1,k2)+Fz2 * (1.0_psn-Wx2)* Wy2
          Jz(i2+1,j2+1,k2)= Jz(i2+1,j2+1,k2)+Fz2 *  Wx2    * Wy2
          
     end subroutine DepositCurrentPIC
     
!----------------------------------------------------------------------------------------------
! The following subroutine store particles sum(Weight*velocity) in J's and sum(Weight) in F0
! It is useful to compute flow velocity of particles
! Important: Make sure that J's and F0 are avaialble and have been reset to 0
!----------------------------------------------------------------------------------------------

subroutine DepCurrDenCellCenter(q,x,y,z,u,v,w)
     real(psn) ::q
     real(psn)::x
     real(psn)::y
     real(psn)::z
     real(psn) ::u,v,w
     real(psn) ::qthis
     real(psn) ::vx,vy,vz
     real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp,invg
     integer :: i,j,k

     qthis=abs(q)

     invg=1.0_psn/sqrt(1.0_psn+u*u+v*v+w*w)
     vx=u*invg
     vy=v*invg
     vz=w*invg


     i=floor(x-0.5_psn)
     j=floor(y-0.5_psn)
#ifndef twoD
     k=floor(z-0.5_psn)
#endif
     Wxp=x-i-0.5_psn
     Wx =1.0_psn-Wxp
     Wyp=y-j-0.5_psn
     Wy =1.0_psn-Wyp

#ifndef twoD
     Wzp=z-k-0.5_psn
     Wz =1-Wzp
#else
    k=1
    Wz=1.0_psn
#endif

                    Jx(i   ,j   ,k  )=Jx(i   ,j   ,k  )+ Wx *Wy *Wz*vx*qthis
                    Jx(i+1 ,j   ,k  )=Jx(i+1 ,j   ,k  )+ Wxp*Wy *Wz*vx*qthis
                    Jx(i   ,j+1 ,k  )=Jx(i   ,j+1 ,k  )+ Wx *Wyp*Wz*vx*qthis
                    Jx(i+1 ,j+1 ,k  )=Jx(i+1 ,j+1 ,k  )+ Wxp*Wyp*Wz*vx*qthis
#ifndef twoD
                    Jx(i  ,j  ,k+1 )=Jx(i  ,j  ,k+1 )+ Wx *Wy *Wzp*vx*qthis
                    Jx(i+1,j  ,k+1 )=Jx(i+1,j  ,k+1 )+ Wxp*Wy *Wzp*vx*qthis
                    Jx(i  ,j+1,k+1 )=Jx(i  ,j+1,k+1 )+ Wx *Wyp*Wzp*vx*qthis
                    Jx(i+1,j+1,k+1 )=Jx(i+1,j+1,k+1 )+ Wxp*Wyp*Wzp*vx*qthis
#endif


                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qthis
                    Jy(i+1,j  ,k  )=Jy(i+1,j  ,k  )+ Wxp*Wy *Wz*vy*qthis
                    Jy(i  ,j+1,k  )=Jy(i  ,j+1,k  )+ Wx *Wyp*Wz*vy*qthis
                    Jy(i+1,j+1,k  )=Jy(i+1,j+1,k  )+ Wxp*Wyp*Wz*vy*qthis
#ifndef twoD
                    Jy(i  ,j  ,k+1 )=Jy(i  ,j  ,k+1 )+ Wx *Wy *Wzp*vy*qthis
                    Jy(i+1,j  ,k+1 )=Jy(i+1,j  ,k+1 )+ Wxp*Wy *Wzp*vy*qthis
                    Jy(i  ,j+1,k+1 )=Jy(i  ,j+1,k+1 )+ Wx *Wyp*Wzp*vy*qthis
                    Jy(i+1,j+1,k+1 )=Jy(i+1,j+1,k+1 )+ Wxp*Wyp*Wzp*vy*qthis
#endif

                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qthis
                    Jz(i+1,j  ,k  )=Jz(i+1,j  ,k  )+ Wxp*Wy *Wz*vz*qthis
                    Jz(i  ,j+1,k  )=Jz(i  ,j+1,k  )+ Wx *Wyp*Wz*vz*qthis
                    Jz(i+1,j+1,k  )=Jz(i+1,j+1,k  )+ Wxp*Wyp*Wz*vz*qthis
#ifndef twoD
                    Jz(i  ,j  ,k+1 )=Jz(i  ,j  ,k+1 )+ Wx *Wy *Wzp*vz*qthis
                    Jz(i+1,j  ,k+1 )=Jz(i+1,j  ,k+1 )+ Wxp*Wy *Wzp*vz*qthis
                    Jz(i  ,j+1,k+1 )=Jz(i  ,j+1,k+1 )+ Wx *Wyp*Wzp*vz*qthis
                    Jz(i+1,j+1,k+1 )=Jz(i+1,j+1,k+1 )+ Wxp*Wyp*Wzp*vz*qthis
#endif

                    F0(i   ,j   ,k  )=F0(i   ,j   ,k  )+ Wx *Wy *Wz*qthis
                    F0(i+1 ,j   ,k  )=F0(i+1 ,j   ,k  )+ Wxp*Wy *Wz*qthis
                    F0(i   ,j+1 ,k  )=F0(i   ,j+1 ,k  )+ Wx *Wyp*Wz*qthis
                    F0(i+1 ,j+1 ,k  )=F0(i+1 ,j+1 ,k  )+ Wxp*Wyp*Wz*qthis
#ifndef twoD
                    F0(i  ,j  ,k+1 )=F0(i  ,j  ,k+1 )+ Wx *Wy *Wzp*qthis
                    F0(i+1,j  ,k+1 )=F0(i+1,j  ,k+1 )+ Wxp*Wy *Wzp*qthis
                    F0(i  ,j+1,k+1 )=F0(i  ,j+1,k+1 )+ Wx *Wyp*Wzp*qthis
                    F0(i+1,j+1,k+1 )=F0(i+1,j+1,k+1 )+ Wxp*Wyp*Wzp*qthis
#endif


end subroutine DepCurrDenCellCenter
     
     
     
end module deposit 