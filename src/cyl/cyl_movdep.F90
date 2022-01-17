module movdep
     use parameters
     use vars 
	 use cyl_vars
     use deposit 
     use interpolation
	 use cyl_bc
#ifdef OPEN_MP 
	 use omp_lib
#endif	 
     implicit none
	 integer, parameter :: VecBlockSize=32
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

contains 
	 subroutine InitMoveDeposit
		  ShortFldArrSizeY=my-1
#ifdef twoD	
          ShortFldArrSizeZ=1
          allocate(VecEx(4,mx*ShortFldArrSizeY),VecEy(4,mx*ShortFldArrSizeY),VecEz(4,mx*ShortFldArrSizeY))
          allocate(VecBx(4,mx*ShortFldArrSizeY),VecBy(4,mx*ShortFldArrSizeY),VecBz(4,mx*ShortFldArrSizeY))
          allocate(VecJ(8,mx*ShortFldArrSizeY,Nthreads))
#else
          ShortFldArrSizeZ=mz-1
          allocate(VecEx(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecEy(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecEz(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ))	
          allocate(VecBx(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecBy(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecBz(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ))	
          allocate(VecJ(12,mx*ShortFldArrSizeY*ShortFldArrSizeZ,Nthreads))
		  	
#endif
VecEx=0
VecEy=0
VecEz=0
VecBx=0
VecBy=0
VecBz=0
	 end subroutine InitMoveDeposit 
	 
	 subroutine ReshapeShortMoverFldArr
	     !Now update  arrays used in mover 
	      deallocate(VecEx,VecEy,VecEz,VecBx,VecBy,VecBz,VecJ)
		  ShortFldArrSizeY=my-1
#ifdef twoD	
	      ShortFldArrSizeZ=1
	      allocate(VecEx(4,mx*ShortFldArrSizeY),VecEy(4,mx*ShortFldArrSizeY),VecEz(4,mx*ShortFldArrSizeY))
	      allocate(VecBx(4,mx*ShortFldArrSizeY),VecBy(4,mx*ShortFldArrSizeY),VecBz(4,mx*ShortFldArrSizeY))
	      allocate(VecJ(8,mx*ShortFldArrSizeY,Nthreads))
#else
	      ShortFldArrSizeZ=mz-1
	      allocate(VecEx(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecEy(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecEz(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ))	
	      allocate(VecBx(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecBy(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ),VecBz(8,mx*ShortFldArrSizeY*ShortFldArrSizeZ))	
	      allocate(VecJ(12,mx*ShortFldArrSizeY*ShortFldArrSizeZ,Nthreads))		  	
#endif
	      VecEx=0
	      VecEy=0
	      VecEz=0
	      VecBx=0
	      VecBy=0
	      VecBz=0
	 end subroutine ReshapeShortMoverFldArr
	 
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
		  integer :: i1,j1,k1,ShortFldArrWidth2
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
qthis=0.0_psn
r_proc=rborders(procxind(proc))-xmin
axis_proc=.false.
if(inc_axis) then
	if(procxind(proc).eq.0) then 
		axis_proc=.true. 		
	end if 
end if 

!call UpdateFldAxis
 

#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 					
					do j1=2,ShortFldArrSizeY
						do i1=2,mx-1
						  short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
				  		  VecEx(1,short_arr_ind)=(Ex(i1-1,j1,k1)+Ex(i1,j1,k1))*0.5_psn
				  		  VecEx(2,short_arr_ind)=(Ex(i1,j1,k1)+Ex(i1+1,j1,k1))*0.5_psn
				  		  VecEx(3,short_arr_ind)=(Ex(i1-1,j1+1,k1)+Ex(i1,j1+1,k1))*0.5_psn
				  		  VecEx(4,short_arr_ind)=(Ex(i1,j1+1,k1)+Ex(i1+1,j1+1,k1))*0.5_psn
#ifndef twoD 		  
	                      VecEx(5,short_arr_ind)=(Ex(i1-1,j1,k1+1)+Ex(i1,j1,k1+1))*0.5_psn
	                      VecEx(6,short_arr_ind)=(Ex(i1,j1,k1+1)+Ex(i1+1,j1,k1+1))*0.5_psn
	                      VecEx(7,short_arr_ind)=(Ex(i1-1,j1+1,k1+1)+Ex(i1,j1+1,k1+1))*0.5_psn
	                      VecEx(8,short_arr_ind)=(Ex(i1,j1+1,k1+1)+Ex(i1+1,j1+1,k1+1))*0.5_psn
#endif 
						end do 
					end do    
				end do	
												
#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 					
					do j1=2,ShortFldArrSizeY
						do i1=2,mx-1
					     short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
  		                 VecEy(1,short_arr_ind)=(Ey(i1,j1-1,k1)+Ey(i1,j1,k1))*0.5_psn
  		                 VecEy(2,short_arr_ind)=(Ey(i1+1,j1-1,k1)+Ey(i1+1,j1,k1))*0.5_psn
  		                 VecEy(3,short_arr_ind)=(Ey(i1,j1,k1)+Ey(i1,j1+1,k1))*0.5_psn
  		                 VecEy(4,short_arr_ind)=(Ey(i1+1,j1,k1)+Ey(i1+1,j1+1,k1))*0.5_psn
#ifndef twoD 		  
  		                 VecEy(5,short_arr_ind)=(Ey(i1,j1-1,k1+1)+Ey(i1,j1,k1+1))*0.5_psn
  		                 VecEy(6,short_arr_ind)=(Ey(i1+1,j1-1,k1+1)+Ey(i1+1,j1,k1+1))*0.5_psn
  		                 VecEy(7,short_arr_ind)=(Ey(i1,j1,k1+1)+Ey(i1,j1+1,k1+1))*0.5_psn
  		                 VecEy(8,short_arr_ind)=(Ey(i1+1,j1,k1+1)+Ey(i1+1,j1+1,k1+1))*0.5_psn
#endif		
						end do 
					end do    
				end do	
				
	
#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 					
					do j1=2,ShortFldArrSizeY
						do i1=2,mx-1			
   					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
#ifndef twoD
						  VecEz(1,short_arr_ind)=(Ez(i1,j1,k1-1)+Ez(i1,j1,k1))*0.5_psn
						  VecEz(2,short_arr_ind)=(Ez(i1+1,j1,k1-1)+Ez(i1+1,j1,k1))*0.5_psn
						  VecEz(3,short_arr_ind)=(Ez(i1,j1+1,k1-1)+Ez(i1,j1+1,k1))*0.5_psn
						  VecEz(4,short_arr_ind)=(Ez(i1+1,j1+1,k1-1)+Ez(i1+1,j1+1,k1))*0.5_psn	  
						  VecEz(5,short_arr_ind)=(Ez(i1,j1,k1)+Ez(i1,j1,k1+1))*0.5_psn
						  VecEz(6,short_arr_ind)=(Ez(i1+1,j1,k1)+Ez(i1+1,j1,k1+1))*0.5_psn
						  VecEz(7,short_arr_ind)=(Ez(i1,j1+1,k1)+Ez(i1,j1+1,k1+1))*0.5_psn
						  VecEz(8,short_arr_ind)=(Ez(i1+1,j1+1,k1)+Ez(i1+1,j1+1,k1+1))*0.5_psn
#else
						  VecEz(1,short_arr_ind)=Ez(i1,j1,k1)
						  VecEz(2,short_arr_ind)=Ez(i1+1,j1,k1)
						  VecEz(3,short_arr_ind)=Ez(i1,j1+1,k1)
						  VecEz(4,short_arr_ind)=Ez(i1+1,j1+1,k1)		 
#endif 			
				
						end do 
					end do    
				end do	
				

#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 			
                do j1=2,ShortFldArrSizeY
		
			       do i1=2,mx-1	
					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1

#ifndef twoD 		  
						  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1-1)+Bx(i1,j1,k1-1)+Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.25_psn	+Bx_ext0
						  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1-1)+Bx(i1+1,j1,k1-1)+Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.25_psn+Bx_ext0  	  
						  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1-1)+Bx(i1,j1+1,k1-1)+Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.25_psn	+Bx_ext0
						  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1-1)+Bx(i1+1,j1+1,k1-1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.25_psn	+Bx_ext0
						  VecBx(5,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1)+Bx(i1,j1-1,k1+1)+Bx(i1,j1,k1+1))*0.25_psn	+Bx_ext0
						  VecBx(6,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1-1,k1+1)+Bx(i1+1,j1,k1+1))*0.25_psn +Bx_ext0 	  
						  VecBx(7,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1)+Bx(i1,j1,k1+1)+Bx(i1,j1+1,k1+1))*0.25_psn	+Bx_ext0
						  VecBx(8,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1)+Bx(i1+1,j1,k1+1)+Bx(i1+1,j1+1,k1+1))*0.25_psn	+Bx_ext0	
#else
						  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.5_psn	+Bx_ext0
						  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.5_psn +Bx_ext0 	  
						  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.5_psn	+Bx_ext0
						  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.5_psn+Bx_ext0		
#endif 
                     end do 
					 end do 
				end do 

#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 			
                do j1=2,ShortFldArrSizeY		
			       do i1=2,mx-1	
					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
#ifndef twoD 		  
						  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1-1)+By(i1-1,j1,k1)+By(i1,j1,k1-1)+By(i1,j1,k1))*0.25_psn +By_ext0
						  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1-1)+By(i1,j1,k1)+By(i1+1,j1,k1-1)+By(i1+1,j1,k1))*0.25_psn +By_ext0
						  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1-1)+By(i1-1,j1+1,k1)+By(i1,j1+1,k1-1)+By(i1,j1+1,k1))*0.25_psn +By_ext0
						  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1-1)+By(i1,j1+1,k1)+By(i1+1,j1+1,k1-1)+By(i1+1,j1+1,k1))*0.25_psn+By_ext0
						  	  VecBy(5,short_arr_ind)=(By(i1-1,j1,k1)+By(i1-1,j1,k1+1)+By(i1,j1,k1)+By(i1,j1,k1+1))*0.25_psn+By_ext0
						  	  VecBy(6,short_arr_ind)=(By(i1,j1,k1)+By(i1,j1,k1+1)+By(i1+1,j1,k1)+By(i1+1,j1,k1+1))*0.25_psn+By_ext0
						  	  VecBy(7,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1-1,j1+1,k1+1)+By(i1,j1+1,k1)+By(i1,j1+1,k1+1))*0.25_psn+By_ext0
						  	  VecBy(8,short_arr_ind)=(By(i1,j1+1,k1)+By(i1,j1+1,k1+1)+By(i1+1,j1+1,k1)+By(i1+1,j1+1,k1+1))*0.25_psn+By_ext0
#else
						  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1)+By(i1,j1,k1))*0.5_psn+By_ext0
						  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1)+By(i1+1,j1,k1))*0.5_psn+By_ext0
						  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1,j1+1,k1))*0.5_psn+By_ext0
						  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1)+By(i1+1,j1+1,k1))*0.5_psn+By_ext0
#endif 

                     end do 
					 end do 
				end do 

#ifndef twoD 				
				do k1=2,ShortFldArrSizeZ
#else 
                do k1=1,1
#endif 				
                do j1=2,ShortFldArrSizeY
			       do i1=2,mx-1	
					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
		  
						  VecBz(1,short_arr_ind)=(Bz(i1-1,j1-1,k1)+Bz(i1-1,j1,k1)+Bz(i1,j1-1,k1)+Bz(i1,j1,k1))*0.25_psn+Bz_ext0
						  VecBz(2,short_arr_ind)=(Bz(i1,j1-1,k1)+Bz(i1,j1,k1)+Bz(i1+1,j1-1,k1)+Bz(i1+1,j1,k1))*0.25_psn+Bz_ext0
						  VecBz(3,short_arr_ind)=(Bz(i1-1,j1,k1)+Bz(i1-1,j1+1,k1)+Bz(i1,j1,k1)+Bz(i1,j1+1,k1))*0.25_psn+Bz_ext0
						  VecBz(4,short_arr_ind)=(Bz(i1,j1,k1)+Bz(i1,j1+1,k1)+Bz(i1+1,j1,k1)+Bz(i1+1,j1+1,k1))*0.25_psn+Bz_ext0
#ifndef twoD 		  
						  VecBz(5,short_arr_ind)=(Bz(i1-1,j1-1,k1+1)+Bz(i1-1,j1,k1+1)+Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1))*0.25_psn+Bz_ext0
						  VecBz(6,short_arr_ind)=(Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1+1,j1-1,k1+1)+Bz(i1+1,j1,k1+1))*0.25_psn+Bz_ext0
						  VecBz(7,short_arr_ind)=(Bz(i1-1,j1,k1+1)+Bz(i1-1,j1+1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1))*0.25_psn+Bz_ext0
						  VecBz(8,short_arr_ind)=(Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1)+Bz(i1+1,j1,k1+1)+Bz(i1+1,j1+1,k1+1))*0.25_psn+Bz_ext0
#endif 
                     end do 
					 end do 
				end do 				
			
				VecJ=0.0_psn ! reset the value of current to zero




!-------------------------------------------------------------------------------------		
! First move all the particles in the Sorted part of the prtl array using optimized algorithm 
!-------------------------------------------------------------------------------------
		


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
		  do nn=1,4 
			  wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
		  end do     
#endif   

          short_arr_ind= (kp0(n)-1)*mx*ShortFldArrSizeY  + (jp0(n)-1)*mx + ip0(n) 
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

               u1(n)=u0(n)*cinv
               v1(n)=v0(n)*cinv
               wp(n+off)=w0(n)*cinv

               g(n)=sqc/sqrt(sqc+u0(n)**2+v0(n)**2+w0(n)**2)

               x1(n)=xp(n+off) + u1(n)*g(n) + r_proc
			   y1(n)=v1(n)*g(n)
			   r1(n)=sqrt(x1(n)**2+y1(n)**2)
			   xp(n+off)=r1(n) - r_proc
               yp(n+off)=yp(n+off) + atan(y1(n)/x1(n))*inv_dtheta
               zp(n+off)=zp(n+off) + wp(n+off)*g(n)
			   
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
#else
		        Wz1=0.5_psn*(z0(n)+zr)-kp0(n)
#endif
                short_arr_ind= (kp0(n)-1)*mx*ShortFldArrSizeY  + (jp0(n)-1)*mx + ip0(n)
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
#else
		        Wz2=0.5_psn*(zp(n+off)+zr)-kp
#endif
                short_arr_ind= (kp-1)*mx*ShortFldArrSizeY  + (jp-1)*mx + ip
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


!-------------Collect Current From All threads ----------!
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(j1,k1)
#endif
do j1=1,size(VecJ,2)
do k1=2,Nthreads
	VecJ(:,j1,1)=VecJ(:,j1,1)+VecJ(:,j1,k1)
end do 
end do 
!-------------------------------------------------------!



! Now deposit the current into the main array
#ifndef twoD
				do k1=1 ,mz-1
#else
                do k1=1,1
#endif
                   do j1=1,my-1
			          do i1=1,mx-1
#ifndef twoD 
					      short_arr_ind= (k1-1)*ShortFldArrSizeY*mx + (j1-1)*mx + i1
				          Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind,1)
						  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind,1)
						  Jx(i1,j1,k1+1)  =Jx(i1,j1,k1+1)  +VecJ(3,short_arr_ind,1)
						  Jx(i1,j1+1,k1+1)=Jx(i1,j1+1,k1+1)+VecJ(4,short_arr_ind,1)

				          Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(5,short_arr_ind,1)
						  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(6,short_arr_ind,1)
						  Jy(i1,j1,k1+1)  =Jy(i1,j1,k1+1)  +VecJ(7,short_arr_ind,1)
						  Jy(i1+1,j1,k1+1)=Jy(i1+1,j1,k1+1)+VecJ(8,short_arr_ind,1)

		                  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(9,short_arr_ind,1)
		                  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(10,short_arr_ind,1)
		                  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(11,short_arr_ind,1)
		                  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(12,short_arr_ind,1)
#else
						  short_arr_ind=  (j1-1)*mx + i1
						  Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind,1)
						  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind,1)

						  Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(3,short_arr_ind,1)
						  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(4,short_arr_ind,1)

						  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(5,short_arr_ind,1)
						  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(6,short_arr_ind,1)
						  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(7,short_arr_ind,1)
						  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(8,short_arr_ind,1)
#endif

					  end do
				   end do
			    end do

	 
	 
	 
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
		 
		 if(qp(n).eq.0) return
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
		 
		 call InterpEBfld(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
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

         x1=x + u1*g*dt - 3.5_psn
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
	         xp(n)=r1 + 3.5_psn
	         up(n)=(x1*u1+y1*v1)/r1
             vp(n)=(-y1*u1+x1*v1)/r1	
		 else
			 xp(n)=x1 + 3.5_psn
	         up(n)=u1
			 vp(n)=v1
		 end if  
         yp(n)=yp(n)+dy
		 
		 wp(n)=w0*cinv
		 zp(n)=zp(n) + wp(n)*g*dt
		 	   
		 res=1 
		 		 
	 end subroutine MovePrtlAxis
	 subroutine InterpEBfld(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
		 real(psn) :: x,y,z
		 real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz,dx,dy,dz
		 integer   :: i,j,k,nn,short_arr_ind
#ifdef twoD 
         real(psn), dimension(4) :: wt
#else
         real(psn), dimension(8) :: wt 
#endif  		
			pEx=0.0_psn
			pEy=0.0_psn
			pEz=0.0_psn
			pBx=0.0_psn
			pBy=0.0_psn
			pBz=0.0_psn	

 		    i=x
 		    j=y
 		    dx=x-i
 		    dy=y-j 
#ifndef twoD 
            k=z
            dz=z-k 	
!$OMP SIMD
	      do nn=1,8
		     wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
	      end do 	  

#else
          k=1 
	      do nn=1,4 
		      wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
	      end do     
#endif   

          short_arr_ind= (k-1)*mx*ShortFldArrSizeY  + (j-1)*mx + i 

!$OMP SIMD
#ifdef twoD
	      do nn=1,4
#else
          do nn=1,8
#endif
 			  pEx=pEx+wt(nn)*VecEx(nn,short_arr_ind)
 			  pEy=pEy+wt(nn)*VecEy(nn,short_arr_ind)
 			  pEz=pEz+wt(nn)*VecEz(nn,short_arr_ind)
 			  pBx=pBx+wt(nn)*VecBx(nn,short_arr_ind)
 			  pBy=pBy+wt(nn)*VecBy(nn,short_arr_ind)
 			  pBz=pBz+wt(nn)*VecBz(nn,short_arr_ind)
		  end do 
	 end subroutine InterpEBFld	 
	 
	 
! 	 subroutine UpdateFldCenter
! 		 integer :: n,ind,j,k
! 		 real(psn) :: phi,wt
!
!  		wt=0
!  		do n=0,ny-1
!  			phi= 0.5*dtheta + n*dtheta
!  			wt=wt+abs(sin(phi))
!  		end do
!  		wt=1.0_psn/wt
!
! #ifdef twoD
!          do k=1,1
! #else
!          do k=1,mz
! #endif
!             do j=1,my
! 				Ex(3,j,k)=0.0_psn
! 				By(3,j,k)=0.0_psn
! 				do n=0,ny-1
! 					phi= 0.5*dtheta + n*dtheta
! 					ind=j+n
! 					if(ind.gt.my-3) ind= ind - ny
! 					Ex(3,j,k)=Ex(3,j,k)-Ey(4,ind,k)*sin(phi)*wt
! 					By(3,j,k)=By(3,j,k)+By(4,ind,k)*sin(phi)*wt
! 				end do
! 			end do
!          end do
! 	 end subroutine UpdateFldCenter
	 

!----------------------------------------------------------------------------------------
! Mover for Test Partciels 
!----------------------------------------------------------------------------------------	 
	 	 
	 subroutine MoveTestPrtl
          integer :: n,off
		  real(psn), dimension(VecBlockSize) :: pEx,pEy,pEz,pBx,pBy,pBz,u0,v0,w0,u1,v1,w1,f,g
          integer :: i,j,k
		  real(psn), dimension(VecBlockSize) :: qm
		  real(psn) :: dx,dy,dz
		  integer :: VecBlockSizeThis
#ifdef twoD 
          real(psn), dimension(4) :: FldVec,wt
#else
          real(psn), dimension(8) :: FldVec,wt 
#endif  
	
qm=0.0_psn 
pEx=0.0_psn 
pEy=0.0_psn 
pEz=0.0_psn
pBx=0.0_psn 
pBy=0.0_psn 
pBz=0.0_psn 	 
		 
off=0
  do while(off.lt.used_test_prtl_arr_size) 
	 VecBlockSizeThis=min(VecBlockSize,used_test_prtl_arr_size-off)

	! load q/m for each particles into a vector
	do n=1,VecBlockSizeThis
	    if(flvtp(n+off).eq.0) then ! if condition can possibly be removed by defining flvrqm(0)=0 
		  qm(n)=0
	    else 			
	      qm(n)=flvrqm(flvtp(n+off))*0.5_psn !for optimisation 0.5 is multiplied here itself 
	    end if 
	end do
     !interpolation 		
	do n=1,VecBlockSizeThis 
	  i=xtp(n) 
	  j=ytp(n)
	  dx=xtp(n)-i
	  dy=ytp(n)-j 
#ifndef twoD 
      k=ztp(n)
      dz=ztp(n)-k 		  
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
	  FldVec(1)=FilteredEx(i-1,j,k)+FilteredEx(i,j,k)
	  FldVec(2)=FilteredEx(i,j,k)+FilteredEx(i+1,j,k)
	  FldVec(3)=FilteredEx(i-1,j+1,k)+FilteredEx(i,j+1,k)
	  FldVec(4)=FilteredEx(i,j+1,k)+FilteredEx(i+1,j+1,k)
#ifndef twoD 		  
	  FldVec(5)=FilteredEx(i-1,j,k+1)+FilteredEx(i,j,k+1)
	  FldVec(6)=FilteredEx(i,j,k+1)+FilteredEx(i+1,j,k+1)
	  FldVec(7)=FilteredEx(i-1,j+1,k+1)+FilteredEx(i,j+1,k+1)
	  FldVec(8)=FilteredEx(i,j+1,k+1)+FilteredEx(i+1,j+1,k+1)
#endif 
	  pEx(n)=sum(0.5_psn*wt*FldVec) 

	  FldVec(1)=FilteredEy(i,j-1,k)+FilteredEy(i,j,k)
	  FldVec(2)=FilteredEy(i+1,j-1,k)+FilteredEy(i+1,j,k)
	  FldVec(3)=FilteredEy(i,j,k)+FilteredEy(i,j+1,k)
	  FldVec(4)=FilteredEy(i+1,j,k)+FilteredEy(i+1,j+1,k)
#ifndef twoD 		  
	  FldVec(5)=FilteredEy(i,j-1,k+1)+FilteredEy(i,j,k+1)
	  FldVec(6)=FilteredEy(i+1,j-1,k+1)+FilteredEy(i+1,j,k+1)
	  FldVec(7)=FilteredEy(i,j,k+1)+FilteredEy(i,j+1,k+1)
	  FldVec(8)=FilteredEy(i+1,j,k+1)+FilteredEy(i+1,j+1,k+1)
#endif
	  pEy(n)=sum(0.5_psn*wt*FldVec)

#ifndef twoD 		  
	  FldVec(1)=FilteredEz(i,j,k-1)+FilteredEz(i,j,k)
	  FldVec(2)=FilteredEz(i+1,j,k-1)+FilteredEz(i+1,j,k)
	  FldVec(3)=FilteredEz(i,j+1,k-1)+FilteredEz(i,j+1,k)
	  FldVec(4)=FilteredEz(i+1,j+1,k-1)+FilteredEz(i+1,j+1,k)		  
	  FldVec(5)=FilteredEz(i,j,k)+FilteredEz(i,j,k+1)
	  FldVec(6)=FilteredEz(i+1,j,k)+FilteredEz(i+1,j,k+1)
	  FldVec(7)=FilteredEz(i,j+1,k)+FilteredEz(i,j+1,k+1)
	  FldVec(8)=FilteredEz(i+1,j+1,k)+FilteredEz(i+1,j+1,k+1)
       pEz(n)=sum(0.5_psn*wt*FldVec)
#else
		FldVec(1)=FilteredEz(i,j,k)
		FldVec(2)=FilteredEz(i+1,j,k)
		FldVec(3)=FilteredEz(i,j+1,k)
		FldVec(4)=FilteredEz(i+1,j+1,k)		 
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

!$OMP SIMD 
    do n=1,VecBlockSizeThis    
           !Boris Pusher 
            u0(n)=c*utp(n+off)+pEx(n)
            v0(n)=c*vtp(n+off)+pEy(n)
            w0(n)=c*wtp(n+off)+pEz(n)

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

            utp(n+off)=u0(n)*cinv
            vtp(n+off)=v0(n)*cinv
            wtp(n+off)=w0(n)*cinv

            g(n)=sqc/sqrt(sqc+u0(n)**2+v0(n)**2+w0(n)**2)

            xtp(n+off)=xtp(n+off) + utp(n+off)*g(n)
            ytp(n+off)=ytp(n+off) + vtp(n+off)*g(n) 
            ztp(n+off)=ztp(n+off) + wtp(n+off)*g(n)			   
     end do        
  
! Periodic Boundary conditions 
#ifdef twoD
     do n=1,VecBlockSizeThis
          if(ztp(n+off).gt.zmax) then
            ztp(n+off)=ztp(n+off)-zlen
          else if(ztp(n+off).lt.zmin) then
            ztp(n+off)=zlen+ztp(n+off)
          end if 
	end do 
#endif         

 off=off+VecBlockSizeThis
 end do ! end of the main while loop 
		 
	 end subroutine MoveTestPrtl
     
     
end module movdep