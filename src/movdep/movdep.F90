module movdep
     use parameters
     use vars 
	 use movdep_routines
#ifdef OPEN_MP 
	 use omp_lib
#endif	 
     implicit none

contains 
	 
     subroutine MoveDepositPrtl
          integer :: n,jc,kc,off,nn,off_max,thread_prtl_block_size
		  real(psn), dimension(VecBlockSize) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,u1,v1,w1,f,g
          real(psn), dimension(VecBlockSize) :: x0,y0,z0
		  integer, dimension(VecBlockSize):: ip0,jp0,kp0
		  integer :: ip,jp,kp
		  real(psn)    :: xr,yr,zr
		  real(psn)    :: Fx1,Fx2,Fy1,Fy2,Fz1,Fz2,Wx1,Wx2,Wy1,Wy2,Wz1,Wz2
		  real(psn) :: qthis
          integer :: i,j,k
		  real(psn), dimension(VecBlockSize) :: qm
		  real(psn) :: dx,dy,dz
		  integer :: VecBlockSizeThis
		  integer :: i1,j1,k1
		  integer :: joff,koff,short_arr_ind
#ifdef twoD 
          real(psn), dimension(4) :: FldVec,wt
#else
          real(psn), dimension(8) :: FldVec,wt 
#endif  

		
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

           call GatherVecEfld(VecEx,VecEy,VecEz,Ex,Ey,Ez) 
			
		   call GatherVecBfld

		   VecJ=0.0_psn ! reset the value of current to zero

!off=ThreadID*thread_prtl_block_size 
!off_max=min(off+thread_prtl_block_size,used_prtl_arr_size)
!-------------------------------------------------------------------------------------		
! Now move particles
!-------------------------------------------------------------------------------------
 				 
#ifndef OPEN_MP
    ThreadID=0	
#else 
!$OMP PARALLEL DO  PRIVATE(off,VecBlockSizeThis,n,x0,y0,z0,qm,pEx,pEy,pEz,pBx,pBy,pBz,ip0,jp0,kp0,dx,dy,dz,nn,wt,short_arr_ind,u0,v0,w0,g,f,u1,v1,w1,ip,jp,kp,xr,yr,zr,qthis,Fx1,Fy1,Fz1,Wx1,Wy1,Wz1,Fx2,Fy2,Fz2,Wx2,Wy2,Wz2,ThreadID) 
#endif 
	 do off=0, used_prtl_arr_size-1,VecBlockSize 
		 VecBlockSizeThis=min(VecBlockSize,used_prtl_arr_size-off)
#ifdef OPEN_MP		 
		 ThreadID=OMP_GET_THREAD_NUM()
#endif 		 
	    !initial position of particles  
		 do n=1,VecBlockSizeThis 
			  x0(n)=xp(n+off)
	          y0(n)=yp(n+off)
	          z0(n)=zp(n+off)
	    end do 
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
		  kp0(n)=z0(n)
		  
		  dx=x0(n)-ip0(n)
		  dy=y0(n)-jp0(n) 
#ifndef twoD      
          dz=z0(n)-kp0(n) 	
!$OMP SIMD
		  do nn=1,8
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
		  end do 	  

#else
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

               up(n+off)=u0(n)*cinv
               vp(n+off)=v0(n)*cinv
               wp(n+off)=w0(n)*cinv

               g(n)=sqc/sqrt(sqc+u0(n)**2+v0(n)**2+w0(n)**2)

               xp(n+off)=xp(n+off) + up(n+off)*g(n)
               yp(n+off)=yp(n+off) + vp(n+off)*g(n)
               zp(n+off)=zp(n+off) + wp(n+off)*g(n)
	   
        end do     



         do n=1,VecBlockSizeThis

 	            ip=xp(n+off)
 	            jp=yp(n+off)
 	            kp=zp(n+off)

                xr=min(real(min(ip0(n),ip)+1,psn),max(real(max(ip0(n),ip),psn),0.5_psn*(x0(n)+xp(n+off))))
                yr=min(real(min(jp0(n),jp)+1,psn),max(real(max(jp0(n),jp),psn),0.5_psn*(y0(n)+yp(n+off))))
                zr=min(real(min(kp0(n),kp)+1,psn),max(real(max(kp0(n),kp),psn),0.5_psn*(z0(n)+zp(n+off))))

                qthis=qp(n+off)*qi   ! q = particle's weight x sign of the charge

 				Fx1=qthis*(xr-x0(n))
 		        Fy1=qthis*(yr-y0(n))
				Fz1=qthis*(zr-z0(n))
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
				

		        Fx2=qthis*(xp(n+off)-xr)
		        Fy2=qthis*(yp(n+off)-yr)
		        Fz2=qthis*(zp(n+off)-zr)
		        Wx2=0.5_psn*(xp(n+off)+xr)-ip
		        Wy2=0.5_psn*(yp(n+off)+yr)-jp
#ifdef twoD
		        Wz2=0.0_psn
				short_arr_ind= (jp-1)*mx + ip
#else
		        Wz2=0.5_psn*(zp(n+off)+zr)-kp
				short_arr_ind= (kp-1)*mx*my  + (jp-1)*mx + ip
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
		
              
! 2D: periodic Boundary conditions in the z-direction 
#ifdef twoD
        do n=1,VecBlockSizeThis
             if(zp(n+off).gt.zmax) then
               zp(n+off)=zp(n+off)-zlen
             else if(zp(n+off).lt.zmin) then
               zp(n+off)=zlen+zp(n+off)
             end if 
		end do 
#endif         
         
    end do ! end of the main loop
	


	call ReduceCurrentMatrix

	 	 
	 
	 end subroutine MoveDepositPrtl
	 

	 
	      
end module movdep
