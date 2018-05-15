module movdep
     use parameters
     use vars 
     use deposit 
     use interpolation
     implicit none
contains 
	 subroutine InitMoveDeposit
	 end subroutine InitMoveDeposit
	 subroutine ReshapeShortMoverFldArr
	 end subroutine ReshapeShortMoverFldArr
	 
     subroutine MoveDepositPrtl
          integer :: n
          real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,u1,v1,w1,f,g
          real(psn) :: x0,y0,z0
          real(psn) :: qm
     do n=1,used_prtl_arr_size
          x0=xp(n)
          y0=yp(n)
          z0=zp(n)
          !determine the charge to mass ratio of the particle           
		  qm=flvrqm(flvp(n))*0.5_psn ! for optimisation 0.5 is multiplied here itself 
		   
          call InterpEMFldGridPointAvg(xp(n),yp(n),zp(n),pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)     
          
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
          
               g=c/sqrt(c**2+u0**2+v0**2+w0**2)

               xp(n)=xp(n) + up(n)*g*c
               yp(n)=yp(n) + vp(n)*g*c 
               zp(n)=zp(n) + wp(n)*g*c
               

               call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n)) 
 
              
! Periodic Boundary conditions 
#ifdef twoD
             if(zp(n).gt.zmax) then
               zp(n)=zp(n)-zlen
             else if(zp(n).lt.zmin) then
               zp(n)=zlen+zp(n)
             end if 
#endif                  

    end do 
                    
     end subroutine MoveDepositPrtl 
	 
     subroutine MoveTestPrtl
          integer :: n
          real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,u1,v1,w1,f,g
          real(psn) :: qm
     do n=1,used_test_prtl_arr_size
          !determine the charge to mass ratio of the particle           
		  qm=flvrqm(flvp(n))*0.5_psn ! for optimisation 0.5 is multiplied here itself 
		   
          call InterpEMFldGridPointAvg(xtp(n),ytp(n),ztp(n),pEx,pEy,pEz,pBx,pBy,pBz,FilteredEx,FilteredEy,FilteredEz,Bx,By,Bz)     
          
          pEx=pEx*qm
          pEy=pEy*qm
          pEz=pEz*qm
          pBx=pBx*qm
          pBy=pBy*qm
          pBz=pBz*qm
          
              !Boris Pusher 
               u0=c*utp(n)+pEx
               v0=c*vtp(n)+pEy
               w0=c*wtp(n)+pEz

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

               utp(n)=u0*cinv
               vtp(n)=v0*cinv
               wtp(n)=w0*cinv
          
               g=c/sqrt(c**2+u0**2+v0**2+w0**2)

               xtp(n)=xtp(n) + utp(n)*g*c
               ytp(n)=ytp(n) + vtp(n)*g*c 
               ztp(n)=ztp(n) + wtp(n)*g*c
              
! Periodic Boundary conditions 
#ifdef twoD
             if(ztp(n).gt.zmax) then
               ztp(n)=ztp(n)-zlen
             else if(ztp(n).lt.zmin) then
               ztp(n)=zlen+ztp(n)
             end if 
#endif                  

    end do 
                    
     end subroutine MoveTestPrtl
     

     
end module movdep