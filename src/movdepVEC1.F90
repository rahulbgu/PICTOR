! -------------------------------------------------------------
! The interpolation is done in a differet way in this routine  
!The foolwing subrroutine is currently not being used
!It is not clear whether it is significantly faster than the other method in use 
! if used, var.F90,Makefile and loadbalance.F90 needs to be changed accordignly  
! --------------------------------------------------------



module movdep
     use parameters
     use vars 
     use deposit 
     use interpolation
     implicit none
	 integer, parameter :: VecBlockSize=32
contains 
	 subroutine InitMoveDeposit
		  ShortFldArrSizeY=my-1
#ifdef twoD	
          ShortFldArrSizeZ=1
          allocate(VecEM(6,mx*my))
          allocate(VecJ(8,mx*my))
#else
          ShortFldArrSizeZ=mz-1
          allocate(VecEM(6,mx*my*mz))
		  allocate(VecJ(16,mx*my*mz))
		  	
#endif

VecEM=0
	 end subroutine InitMoveDeposit 
	 
     subroutine MoveDepositPrtl
          integer :: n,jc,kc,off,nn,mm
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
		  integer :: i1,j1,k1,ShortFldArrWidth2
		  integer :: joff,koff,short_arr_ind
		  real(psn), dimension(6) :: tempEM
#ifdef twoD 
          real(psn), dimension(4) :: FldVec,wt
		  integer, dimension(4) :: arr_ind
#else
          real(psn), dimension(8) :: FldVec,wt 
		  integer, dimension(8) :: arr_ind
#endif  

#ifdef twoD
          real(psn), dimension(4), parameter :: wtx1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
          real(psn), dimension(4), parameter :: wtx2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
          real(psn), dimension(4), parameter :: wty1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
          real(psn), dimension(4), parameter :: wty2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
		  
		real(psn), dimension(8), parameter :: wtFx= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(8), parameter :: wtFy= (/0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(8), parameter :: wtFz= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
		
		real(psn), dimension(8), parameter :: Jwtx1= (/1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
		real(psn), dimension(8), parameter :: Jwtx2= (/0.0_psn, 0.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
		real(psn), dimension(8), parameter :: Jwty1= (/1.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(8), parameter :: Jwty2= (/-1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
		real(psn), dimension(8), parameter :: Jwtz1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
		real(psn), dimension(8), parameter :: Jwtz2= (/-1.0_psn, -1.0_psn, -1.0_psn, -1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
#else 
		  real(psn), dimension(8), parameter :: wtx1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
		  real(psn), dimension(8), parameter :: wtx2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
		  real(psn), dimension(8), parameter :: wty1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn  /)
		  real(psn), dimension(8), parameter :: wty2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
		  real(psn), dimension(8), parameter :: wtz1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn  /)
		  real(psn), dimension(8), parameter :: wtz2= (/-1.0_psn, -1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn  /)
		  
		real(psn), dimension(12), parameter :: wtFx= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(12), parameter :: wtFy= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(12), parameter :: wtFz= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
		
		real(psn), dimension(12), parameter :: Jwtx1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
		real(psn), dimension(12), parameter :: Jwtx2= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
		real(psn), dimension(12), parameter :: Jwty1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn,  1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
		real(psn), dimension(12), parameter :: Jwty2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn/)
		real(psn), dimension(12), parameter :: Jwtz1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
		real(psn), dimension(12), parameter :: Jwtz2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
		  
#endif 

		
!-------------------------------------------------------------------------------------		
! First Initialize All the fields needed in the mover 
!-------------------------------------------------------------------------------------

#ifndef twoD 				
				do k1=2,mz
#else 
                do k1=1,1
#endif 					
					do j1=2,my
						do i1=2,mx
						  short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
				  		  VecEM(1,short_arr_ind)=(Ex(i1-1,j1,k1)+Ex(i1,j1,k1))*0.5_psn
   		                  VecEM(2,short_arr_ind)=(Ey(i1,j1-1,k1)+Ey(i1,j1,k1))*0.5_psn
#ifdef twoD 	
                          VecEM(3,short_arr_ind)=Ez(i1,j1,k1)
						  VecEM(4,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.5_psn+Bx_ext0
					  	  VecEM(5,short_arr_ind)=(By(i1-1,j1,k1)+By(i1,j1,k1))*0.5_psn+By_ext0						  							  
#else 					  
						  VecEM(3,short_arr_ind)=(Ez(i1,j1,k1-1)+Ez(i1,j1,k1))*0.5_psn
						  VecEM(4,short_arr_ind)=(Bx(i1,j1-1,k1-1)+Bx(i1,j1,k1-1)+Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.25_psn	+Bx_ext0						  
					  	  VecEM(5,short_arr_ind)=(By(i1-1,j1,k1-1)+By(i1-1,j1,k1)+By(i1,j1,k1-1)+By(i1,j1,k1))*0.25_psn +By_ext0
#endif 				
                          VecEM(6,short_arr_ind)=(Bz(i1-1,j1-1,k1)+Bz(i1-1,j1,k1)+Bz(i1,j1-1,k1)+Bz(i1,j1,k1))*0.25_psn	+Bz_ext0	 
 
						end do 
					end do    
				end do	
												
! #ifndef twoD
! 				do k1=2,ShortFldArrSizeZ
! #else
!                 do k1=1,1
! #endif
! 					do j1=2,ShortFldArrSizeY
! 						do i1=2,mx-1
! 					     short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
!   		                 VecEy(1,short_arr_ind)=(Ey(i1,j1-1,k1)+Ey(i1,j1,k1))*0.5_psn
!   		                 VecEy(2,short_arr_ind)=(Ey(i1+1,j1-1,k1)+Ey(i1+1,j1,k1))*0.5_psn
!   		                 VecEy(3,short_arr_ind)=(Ey(i1,j1,k1)+Ey(i1,j1+1,k1))*0.5_psn
!   		                 VecEy(4,short_arr_ind)=(Ey(i1+1,j1,k1)+Ey(i1+1,j1+1,k1))*0.5_psn
! #ifndef twoD
!   		                 VecEy(5,short_arr_ind)=(Ey(i1,j1-1,k1+1)+Ey(i1,j1,k1+1))*0.5_psn
!   		                 VecEy(6,short_arr_ind)=(Ey(i1+1,j1-1,k1+1)+Ey(i1+1,j1,k1+1))*0.5_psn
!   		                 VecEy(7,short_arr_ind)=(Ey(i1,j1,k1+1)+Ey(i1,j1+1,k1+1))*0.5_psn
!   		                 VecEy(8,short_arr_ind)=(Ey(i1+1,j1,k1+1)+Ey(i1+1,j1+1,k1+1))*0.5_psn
! #endif
! 						end do
! 					end do
! 				end do
!
!
! #ifndef twoD
! 				do k1=2,ShortFldArrSizeZ
! #else
!                 do k1=1,1
! #endif
! 					do j1=2,ShortFldArrSizeY
! 						do i1=2,mx-1
!    					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
! #ifndef twoD
! 						  VecEz(1,short_arr_ind)=(Ez(i1,j1,k1-1)+Ez(i1,j1,k1))*0.5_psn
! 						  VecEz(2,short_arr_ind)=(Ez(i1+1,j1,k1-1)+Ez(i1+1,j1,k1))*0.5_psn
! 						  VecEz(3,short_arr_ind)=(Ez(i1,j1+1,k1-1)+Ez(i1,j1+1,k1))*0.5_psn
! 						  VecEz(4,short_arr_ind)=(Ez(i1+1,j1+1,k1-1)+Ez(i1+1,j1+1,k1))*0.5_psn
! 						  VecEz(5,short_arr_ind)=(Ez(i1,j1,k1)+Ez(i1,j1,k1+1))*0.5_psn
! 						  VecEz(6,short_arr_ind)=(Ez(i1+1,j1,k1)+Ez(i1+1,j1,k1+1))*0.5_psn
! 						  VecEz(7,short_arr_ind)=(Ez(i1,j1+1,k1)+Ez(i1,j1+1,k1+1))*0.5_psn
! 						  VecEz(8,short_arr_ind)=(Ez(i1+1,j1+1,k1)+Ez(i1+1,j1+1,k1+1))*0.5_psn
! #else
! 						  VecEz(1,short_arr_ind)=Ez(i1,j1,k1)
! 						  VecEz(2,short_arr_ind)=Ez(i1+1,j1,k1)
! 						  VecEz(3,short_arr_ind)=Ez(i1,j1+1,k1)
! 						  VecEz(4,short_arr_ind)=Ez(i1+1,j1+1,k1)
! #endif
!
! 						end do
! 					end do
! 				end do
!
!
! #ifndef twoD
! 				do k1=2,ShortFldArrSizeZ
! #else
!                 do k1=1,1
! #endif
!                 do j1=2,ShortFldArrSizeY
!
! 			       do i1=2,mx-1
! 					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
!
! #ifndef twoD
! 						  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1-1)+Bx(i1,j1,k1-1)+Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.25_psn
! 						  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1-1)+Bx(i1+1,j1,k1-1)+Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.25_psn
! 						  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1-1)+Bx(i1,j1+1,k1-1)+Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.25_psn
! 						  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1-1)+Bx(i1+1,j1+1,k1-1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.25_psn
! 						  VecBx(5,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1)+Bx(i1,j1-1,k1+1)+Bx(i1,j1,k1+1))*0.25_psn
! 						  VecBx(6,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1-1,k1+1)+Bx(i1+1,j1,k1+1))*0.25_psn
! 						  VecBx(7,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1)+Bx(i1,j1,k1+1)+Bx(i1,j1+1,k1+1))*0.25_psn
! 						  VecBx(8,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1)+Bx(i1+1,j1,k1+1)+Bx(i1+1,j1+1,k1+1))*0.25_psn
! #else
! 						  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.5_psn
! 						  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.5_psn
! 						  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.5_psn
! 						  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.5_psn
! #endif
!                      end do
! 					 end do
! 				end do
!
! #ifndef twoD
! 				do k1=2,ShortFldArrSizeZ
! #else
!                 do k1=1,1
! #endif
!                 do j1=2,ShortFldArrSizeY
! 			       do i1=2,mx-1
! 					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
! #ifndef twoD
! 						  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1-1)+By(i1-1,j1,k1)+By(i1,j1,k1-1)+By(i1,j1,k1))*0.25_psn
! 						  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1-1)+By(i1,j1,k1)+By(i1+1,j1,k1-1)+By(i1+1,j1,k1))*0.25_psn
! 						  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1-1)+By(i1-1,j1+1,k1)+By(i1,j1+1,k1-1)+By(i1,j1+1,k1))*0.25_psn
! 						  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1-1)+By(i1,j1+1,k1)+By(i1+1,j1+1,k1-1)+By(i1+1,j1+1,k1))*0.25_psn
! 						  	  VecBy(5,short_arr_ind)=(By(i1-1,j1,k1)+By(i1-1,j1,k1+1)+By(i1,j1,k1)+By(i1,j1,k1+1))*0.25_psn
! 						  	  VecBy(6,short_arr_ind)=(By(i1,j1,k1)+By(i1,j1,k1+1)+By(i1+1,j1,k1)+By(i1+1,j1,k1+1))*0.25_psn
! 						  	  VecBy(7,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1-1,j1+1,k1+1)+By(i1,j1+1,k1)+By(i1,j1+1,k1+1))*0.25_psn
! 						  	  VecBy(8,short_arr_ind)=(By(i1,j1+1,k1)+By(i1,j1+1,k1+1)+By(i1+1,j1+1,k1)+By(i1+1,j1+1,k1+1))*0.25_psn
! #else
! 						  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1)+By(i1,j1,k1))*0.5_psn
! 						  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1)+By(i1+1,j1,k1))*0.5_psn
! 						  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1,j1+1,k1))*0.5_psn
! 						  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1)+By(i1+1,j1+1,k1))*0.5_psn
! #endif
!
!                      end do
! 					 end do
! 				end do
!
! #ifndef twoD
! 				do k1=2,ShortFldArrSizeZ
! #else
!                 do k1=1,1
! #endif
!                 do j1=2,ShortFldArrSizeY
! 			       do i1=2,mx-1
! 					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
!
! 						  VecBz(1,short_arr_ind)=(Bz(i1-1,j1-1,k1)+Bz(i1-1,j1,k1)+Bz(i1,j1-1,k1)+Bz(i1,j1,k1))*0.25_psn
! 						  VecBz(2,short_arr_ind)=(Bz(i1,j1-1,k1)+Bz(i1,j1,k1)+Bz(i1+1,j1-1,k1)+Bz(i1+1,j1,k1))*0.25_psn
! 						  VecBz(3,short_arr_ind)=(Bz(i1-1,j1,k1)+Bz(i1-1,j1+1,k1)+Bz(i1,j1,k1)+Bz(i1,j1+1,k1))*0.25_psn
! 						  VecBz(4,short_arr_ind)=(Bz(i1,j1,k1)+Bz(i1,j1+1,k1)+Bz(i1+1,j1,k1)+Bz(i1+1,j1+1,k1))*0.25_psn
! #ifndef twoD
! 						  VecBz(5,short_arr_ind)=(Bz(i1-1,j1-1,k1+1)+Bz(i1-1,j1,k1+1)+Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1))*0.25_psn
! 						  VecBz(6,short_arr_ind)=(Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1+1,j1-1,k1+1)+Bz(i1+1,j1,k1+1))*0.25_psn
! 						  VecBz(7,short_arr_ind)=(Bz(i1-1,j1,k1+1)+Bz(i1-1,j1+1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1))*0.25_psn
! 						  VecBz(8,short_arr_ind)=(Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1)+Bz(i1+1,j1,k1+1)+Bz(i1+1,j1+1,k1+1))*0.25_psn
! #endif
!                      end do
! 					 end do
! 				end do				
			
				VecJ=0.0_psn ! reset the value of current to zero




!-------------------------------------------------------------------------------------		
! First move all the particles in the Sorted part of the prtl array using optimized algorithm 
!-------------------------------------------------------------------------------------
		
		
off=0	  
do kc=1,mz 
	do jc=1,my 	 	
		!print*,'off',off,'Sorted',SortedPrtlCountYZ(jc,kc)
		
	    !First Step is to Manange Short Arrays  	
        koff=kc-1
		joff=jc-1
	     
	!-------------------------------------------------------------------------------------		
	! Now Move the Particles 
	!-------------------------------------------------------------------------------------
	 do while(off.lt.SortedPrtlCountYZ(jc,kc)) 
		 VecBlockSizeThis=min(VecBlockSize,SortedPrtlCountYZ(jc,kc)-off)
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
! 	  pEx=0.0_psn
! 	  pEy=0.0_psn
! 	  pEz=0.0_psn
! 	  pBx=0.0_psn
! 	  pBy=0.0_psn
! 	  pBz=0.0_psn				
		do n=1,VecBlockSizeThis 
		  ip0(n)=x0(n)
		  jp0(n)=y0(n)
		  dx=x0(n)-ip0(n)
		  dy=y0(n)-jp0(n) 
#ifndef twoD 
          kp0(n)=z0(n)
          dz=z0(n)-kp0(n) 	
		  do nn=1,8
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
		  end do 	  
! 	      wt(1)=(1.0_psn-dx)*(1.0_psn-dy)*(1.0_psn-dz)
! 	      wt(2)=dx          *(1.0_psn-dy)*(1.0_psn-dz)
! 	      wt(3)=(1.0_psn-dx) *dy          *(1.0_psn-dz)
! 	      wt(4)=dx          *dy          *(1.0_psn-dz)
! 	      wt(5)=(1.0_psn-dx)*(1.0_psn-dy)*dz
! 	      wt(6)=dx          *(1.0_psn-dy)*dz
! 	      wt(7)=(1.0_psn-dx)*dy          *dz
! 	      wt(8)=dx          *dy          *dz
#else
          kp0(n)=1 
		  do nn=1,4 
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
		  end do     
#endif   

          arr_ind(1)= (kp0(n)-1)*mx*my  + (jp0(n)-1)*mx + ip0(n) 
		  arr_ind(2)=arr_ind(1)+1 
		  arr_ind(3)=arr_ind(1)+mx
		  arr_ind(4)=arr_ind(1)+mx+1 
          tempEM=0.0_psn

!$OMP SIMD
#ifdef twoD
		  do nn=1,4
#else
          do nn=1,8
#endif
			        tempEM=tempEM+wt(nn)*VecEM(1:6,arr_ind(nn))
		  end do

		  pEx(n)=tempEM(1)
		  pEy(n)=tempEM(2)
		  pEz(n)=tempEM(3)
		  pBx(n)=tempEM(4)
		  pBy(n)=tempEM(5)
		  pBz(n)=tempEM(6)
		  


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
 	        call DepositCurrentPIC(x0(n),y0(n),z0(n),xp(n+off),yp(n+off),zp(n+off),qp(n+off)) !deposit charges
 		end do


!          do n=1,VecBlockSizeThis
!
!  	          ip=xp(n+off)
!  	          jp=yp(n+off)
!  	          kp=zp(n+off)
! #ifdef twoD
!               kp=1
! #endif
!                 xr=min(real(min(ip0(n),ip)+1),max(real(max(ip0(n),ip)),0.5_psn*(x0(n)+xp(n+off))))
!                 yr=min(real(min(jp0(n),jp)+1),max(real(max(jp0(n),jp)),0.5_psn*(y0(n)+yp(n+off))))
!                 zr=min(real(min(kp0(n),kp)+1),max(real(max(kp0(n),kp)),0.5_psn*(z0(n)+zp(n+off))))
!
!                 qthis=qp(n+off)*qi   ! q = particle's weight X sign of the charge
!
!  				Fx1=qthis*(xr-x0(n))
!  		        Fy1=qthis*(yr-y0(n))
! 				Fz1=qthis*(zr-z0(n))
!  		        Wx1=0.5_psn*(x0(n)+xr)-ip0(n)
!  		        Wy1=0.5_psn*(y0(n)+yr)-jp0(n)
! #ifdef twoD
! 		        Wz1=0.0_psn
! #else
! 		        Wz1=0.5_psn*(z0(n)+zr)-kp0(n)
! #endif
!                 short_arr_ind= (kp0(n)-1)*mx*my  + (jp0(n)-1)*mx + ip0(n)
!
! !$OMP SIMD
! #ifdef twoD
!                 do nn=1,8
! #else
!                 do nn=1,12
! #endif
!                        VecJ(nn,short_arr_ind)=VecJ(nn,short_arr_ind)+(wtFx(nn)*Fx1+wtFy(nn)*Fy1+wtFz(nn)*Fz1)*(Jwtx1(nn)+Jwtx2(nn)*Wx1)*(Jwty1(nn)+Jwty2(nn)*Wy1)*(Jwtz1(nn)+Jwtz2(nn)*Wz1)
!                 end do
!
!
! 		        Fx2=qthis*(xp(n+off)-xr)
! 		        Fy2=qthis*(yp(n+off)-yr)
! 		        Fz2=qthis*(zp(n+off)-zr)
! 		        Wx2=0.5_psn*(xp(n+off)+xr)-ip
! 		        Wy2=0.5_psn*(yp(n+off)+yr)-jp
! #ifdef twoD
! 		        Wz2=0.0_psn
! #else
! 		        Wz2=0.5_psn*(zp(n+off)+zr)-kp
! #endif
!                 short_arr_ind= (kp-1)*mx*my  + (jp-1)*mx + ip
! !$OMP SIMD
! #ifdef twoD
!                 do nn=1,8
! #else
!                 do nn=1,12
! #endif
!
! !if(proc.eq.0) print*,'F',(Jwtx1(nn)+Jwtx2(nn)*Wx2),(Jwty1(nn)+Jwty2(nn)*Wy2),(Jwtz1(nn)+Jwtz2(nn)*Wz2),nn
! 					   VecJ(nn,short_arr_ind)=VecJ(nn,short_arr_ind)+(wtFx(nn)*Fx2+wtFy(nn)*Fy2+wtFz(nn)*Fz2)*(Jwtx1(nn)+Jwtx2(nn)*Wx2)*(Jwty1(nn)+Jwty2(nn)*Wy2)*(Jwtz1(nn)+Jwtz2(nn)*Wz2)
! 				end do
!
!
!
!
!
!
!
!
!
!
!
!
! ! 		          Wx2(n)=0.5_psn*(xp(n+off)+xr(n))-ip(n)
! ! 		          Wy2(n)=0.5_psn*(yp(n+off)+yr(n))-jp(n)
! ! #ifdef twoD
! ! 		          Wz1(n)=0.0_psn
! ! 		          Wz2(n)=0.0_psn
! ! #else
! ! 		          Wz1(n)=0.5_psn*(z0(n)+zr(n))-kp0(n)
! ! 		          Wz2(n)=0.5_psn*(zp(n+off)+zr(n))-kp(n)
! ! #endif
! ! 		          Fx2(n)=qthis*(xp(n+off)-xr(n))
! ! 		          Fy2(n)=qthis*(yp(n+off)-yr(n))
! ! 		          Fz2(n)=qthis*(zp(n+off)-zr(n))
! ! 	     end do
! !
! ! 		 do n=1,VecBlockSizeThis
! !
! !              Jx(ip0(n),jp0(n),  kp0(n)  )= Jx(ip0(n),jp0(n),  kp0(n)  )+Fx1(n) * (1.0_psn-Wy1(n))*(1.0_psn-Wz1(n))
! !              Jx(ip0(n),jp0(n)+1,kp0(n)  )= Jx(ip0(n),jp0(n)+1,kp0(n)  )+Fx1(n) *  Wy1(n)    *(1.0_psn-Wz1(n))
! !              Jx(ip0(n),jp0(n),  kp0(n)+1)= Jx(ip0(n),jp0(n),  kp0(n)+1)+Fx1(n) * (1.0_psn-Wy1(n)) * Wz1(n)
! !              Jx(ip0(n),jp0(n)+1,kp0(n)+1)= Jx(ip0(n),jp0(n)+1,kp0(n)+1)+Fx1(n) *  Wy1(n)    * Wz1(n)
! !
!
! !
! !
! !              Jy(ip0(n)  ,jp0(n),kp0(n)  )= Jy(ip0(n)  ,jp0(n),kp0(n)  )+Fy1(n) * (1.0_psn-Wx1(n))*(1.0_psn-Wz1(n))
! !              Jy(ip0(n)+1,jp0(n),kp0(n)  )= Jy(ip0(n)+1,jp0(n),kp0(n)  )+Fy1(n) *  Wx1(n)    *(1.0_psn-Wz1(n))
! !              Jy(ip0(n)  ,jp0(n),kp0(n)+1)= Jy(ip0(n)  ,jp0(n),kp0(n)+1)+Fy1(n) * (1.0_psn-Wx1(n))* Wz1(n)
! !              Jy(ip0(n)+1,jp0(n),kp0(n)+1)= Jy(ip0(n)+1,jp0(n),kp0(n)+1)+Fy1(n) *  Wx1(n)    * Wz1(n)
!
! !
! !
! !              Jz(ip0(n)  ,jp0(n)  ,kp0(n))= Jz(ip0(n)  ,jp0(n)  ,kp0(n))+Fz1(n) * (1.0_psn-Wx1(n))*(1.0_psn-Wy1(n))
! !              Jz(ip0(n)+1,jp0(n)  ,kp0(n))= Jz(ip0(n)+1,jp0(n)  ,kp0(n))+Fz1(n) *  Wx1(n)    *(1.0_psn-Wy1(n))
! !              Jz(ip0(n)  ,jp0(n)+1,kp0(n))= Jz(ip0(n)  ,jp0(n)+1,kp0(n))+Fz1(n) * (1.0_psn-Wx1(n))* Wy1(n)
! !              Jz(ip0(n)+1,jp0(n)+1,kp0(n))= Jz(ip0(n)+1,jp0(n)+1,kp0(n))+Fz1(n) *  Wx1(n)    * Wy1(n)
! !
!
!   		end do
		
              
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
         
    off=off+VecBlockSizeThis
    end do ! end of the main while loop 
                    
end do 
end do



! ! Now deposit the current into the main array
! #ifndef twoD
! 				do k1=1 ,mz-1
! #else
!                 do k1=1,1
! #endif
!                    do j1=1,my-1
! 			          do i1=1,mx-1
! #ifndef twoD
! 					      short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
! 				          Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind)
! 						  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind)
! 						  Jx(i1,j1,k1+1)  =Jx(i1,j1,k1+1)  +VecJ(3,short_arr_ind)
! 						  Jx(i1,j1+1,k1+1)=Jx(i1,j1+1,k1+1)+VecJ(4,short_arr_ind)
!
! 				          Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(5,short_arr_ind)
! 						  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(6,short_arr_ind)
! 						  Jy(i1,j1,k1+1)  =Jy(i1,j1,k1+1)  +VecJ(7,short_arr_ind)
! 						  Jy(i1+1,j1,k1+1)=Jy(i1+1,j1,k1+1)+VecJ(8,short_arr_ind)
!
! 		                  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(9,short_arr_ind)
! 		                  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(10,short_arr_ind)
! 		                  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(11,short_arr_ind)
! 		                  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(12,short_arr_ind)
! #else
! 						  short_arr_ind=  (j1-1)*mx + i1
! 						  Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind)
! 						  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind)
!
! 						  Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(3,short_arr_ind)
! 						  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(4,short_arr_ind)
!
! 						  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(5,short_arr_ind)
! 						  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(6,short_arr_ind)
! 						  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(7,short_arr_ind)
! 						  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(8,short_arr_ind)
! #endif
!
! 					  end do
! 				   end do
! 			    end do

  

!-------------------------------------------------------------------------------------		
! Now move all the other particles that are not sorted yet 
!-------------------------------------------------------------------------------------
off=SortedPrtlCountYZ(my,mz)
     do while(off.lt.used_prtl_arr_size) 
		 VecBlockSizeThis=min(VecBlockSize,used_prtl_arr_size-off)
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
		do n=1,VecBlockSizeThis 
		  i=x0(n) 
		  if(i.eq.0) print*,'i',i,x0(n),qp(n+off),xp(n+off),yp(n+off),zp(n+off),up(n+off),vp(n+off),var1p(n+off)
		  j=y0(n)
		  dx=x0(n)-i
		  dy=y0(n)-j 
#ifndef twoD 
          k=z0(n)
          dz=z0(n)-k 		  
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
		  pEx(n)=sum(0.5_psn*wt*FldVec) 
		  
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
		  pEy(n)=sum(0.5_psn*wt*FldVec)

#ifndef twoD 		  
		  FldVec(1)=Ez(i,j,k-1)+Ez(i,j,k)
		  FldVec(2)=Ez(i+1,j,k-1)+Ez(i+1,j,k)
		  FldVec(3)=Ez(i,j+1,k-1)+Ez(i,j+1,k)
		  FldVec(4)=Ez(i+1,j+1,k-1)+Ez(i+1,j+1,k)		  
		  FldVec(5)=Ez(i,j,k)+Ez(i,j,k+1)
		  FldVec(6)=Ez(i+1,j,k)+Ez(i+1,j,k+1)
		  FldVec(7)=Ez(i,j+1,k)+Ez(i,j+1,k+1)
		  FldVec(8)=Ez(i+1,j+1,k)+Ez(i+1,j+1,k+1)
          pEz(n)=sum(0.5_psn*wt*FldVec)
#else
			FldVec(1)=Ez(i,j,k)
			FldVec(2)=Ez(i+1,j,k)
			FldVec(3)=Ez(i,j+1,k)
			FldVec(4)=Ez(i+1,j+1,k)		 
            pEz(n)=sum(wt*FldVec)			 
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
          pBx(n)=sum(0.25_psn*wt*FldVec) + Bx_ext0
#else
			FldVec(1)=Bx(i,j-1,k)+Bx(i,j,k)	
			FldVec(2)=Bx(i+1,j-1,k)+Bx(i+1,j,k)  	  
			FldVec(3)=Bx(i,j,k)+Bx(i,j+1,k)	
			FldVec(4)=Bx(i+1,j,k)+Bx(i+1,j+1,k)		
			pBx(n)=sum(0.5_psn*wt*FldVec) +Bx_ext0	    
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
          pBy(n)=sum(0.25_psn*wt*FldVec) +By_ext0
#else
			FldVec(1)=By(i-1,j,k)+By(i,j,k)
			FldVec(2)=By(i,j,k)+By(i+1,j,k)
			FldVec(3)=By(i-1,j+1,k)+By(i,j+1,k)
			FldVec(4)=By(i,j+1,k)+By(i+1,j+1,k)
			pBy(n)=sum(0.5_psn*wt*FldVec) +By_ext0
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
		  pBz(n)=sum(0.25*wt*FldVec) +Bz_ext0			  
						
		end do 
		
            pEx=pEx*qm
            pEy=pEy*qm
            pEz=pEz*qm
            pBx=pBx*qm
            pBy=pBy*qm
            pBz=pBz*qm      

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
          
               g(n)=c/sqrt(c**2+u0(n)**2+v0(n)**2+w0(n)**2)

               xp(n+off)=xp(n+off) + up(n+off)*g(n)*c
               yp(n+off)=yp(n+off) + vp(n+off)*g(n)*c 
               zp(n+off)=zp(n+off) + wp(n+off)*g(n)*c			   
        end do        

 		do n=1,VecBlockSizeThis
 	        call DepositCurrentPIC(x0(n),y0(n),z0(n),xp(n+off),yp(n+off),zp(n+off),qp(n+off)) !deposit charges
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
       
    off=off+VecBlockSizeThis
    end do ! end of the main while loop 

	 
	 
	 
	 end subroutine MoveDepositPrtl
	 
	 subroutine MoveTestPrtl
	 end subroutine MoveTestPrtl
     

     
!      subroutine DepositCurrent(x0,y0,z0,x,y,z,q)
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
!      end subroutine DepositCurrent
     
end module movdep