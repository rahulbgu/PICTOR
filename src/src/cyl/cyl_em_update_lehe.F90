module em_update
    use parameters
    use vars
    implicit none

	real(psn), parameter :: beta_xy = 0.0_psn
    real(psn), parameter :: beta_xz = 0.125_psn
    real(psn), parameter :: beta_yx = 0.0_psn
    real(psn), parameter :: beta_yz = 0.125_psn
    real(psn), parameter :: beta_zx = 0.125_psn*(grid_dz/grid_dx)**2
    real(psn), parameter :: beta_zy = 0.125_psn*(grid_dz/(grid_dx*dtheta))**2
    real(psn), parameter :: delta_x = 0 
    real(psn), parameter :: delta_y = 0 
    real(psn), parameter :: delta_z = 0.25_psn*( 1.0_psn - (cinv**2) * (sin(0.5_psn*pi*c))**2 )
    real(psn), parameter :: alpha_x = 1.0_psn - 2.0_psn*beta_xy - 2.0_psn*beta_xz - 3.0_psn*delta_x 
    real(psn), parameter :: alpha_y = 1.0_psn - 2.0_psn*beta_yx - 2.0_psn*beta_yz - 3.0_psn*delta_y
    real(psn), parameter :: alpha_z = 1.0_psn - 2.0_psn*beta_zx - 2.0_psn*beta_zy - 3.0_psn*delta_z
 

    private :: gradX, gradY, gradZ, gradX1

contains 
	
     subroutine UpdateEfield
		 real(psn) :: rmin,r,rp_half,rm_half
		 integer :: i,j,k
         integer :: i1,k1,k2
	
		 i1=3
		 if(procxind.eq.0) i1=4

		 rmin=xborders(procxind) + rshift
#ifndef twoD
          k1=3
          k2=mz-3
#else
          k1=1
          k2=1
#endif

        do k=k1,k2
             do j=3,my-3
                  do i=i1,mx-3
					r=( (i-3.0_psn)+rmin )*grid_dx ! in the code units
					rp_half=r+0.5_psn*grid_dx 
					rm_half=r-0.5_psn*grid_dx

                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))/(rp_half*dtheta) - fldc*grid_inv_dz*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*grid_inv_dz*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*grid_inv_dx*(rp_half*By(i,j,k)-rm_half*By(i-1,j,k))/r - fldc*(Bx(i,j,k)-Bx(i,j-1,k))/(r*dtheta)

					end do
               end do
          end do
		  
		  !update E fld near the axis
		  call UpdateEFldAxis
		  
     end subroutine UpdateEfield 
	 
	 
     subroutine UpdateBfieldHalf
		  real(psn) :: rmin,r,rp_half,rp 
		  integer ::i,j,k
          integer ::i1,k1,k2    
		  real(psn), dimension(mz) :: Bz_ax
	
 		  i1=1
 		  if(procxind.eq.0) i1=4

		  rmin=xborders(procxind)+rshift
#ifndef twoD
          k1=1
          k2=mz-1
#else
          k1=1
          k2=1
#endif
          do k=k1,k2
               do j=1,my-1
                    do i=i1,mx-1
						r=((i-3.0_psn)+rmin)*grid_dx
						rp_half=r+0.5_psn*grid_dx
						rp=r+grid_dx

 						 Bx(i,j,k)=Bx(i,j,k)-fld_halfc*( gradY(Ez,i,j,k) )/(r*dtheta) + fld_halfc*grid_inv_dz*(gradZ(Ey,i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*grid_inv_dz*(gradZ(Ex,i,j,k))+fld_halfc*grid_inv_dx*(gradX(Ez,i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*( grid_inv_dx*gradX1(Ey,i,j,k,rp,r)  - gradY(Ex,i,j,k) )/dtheta )/rp_half
                    end do
               end do
          end do

		  if(inc_axis) then
			  if(procxind.eq.0) then
		          
! 				  !update By at the axis using the standard FDTD equation
! 				  i=3
! 				  do k=k1,k2
! 		               do j=1,my-1
! 						    By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
! 					   end do
! 				  end do
				  
				  !update Bz using the integral equation
				  call UpdateBzAxis(Bz_ax)
				  do k=k1,k2
					  do j=1,my-1
						   Bz(3,j,k)=Bz_ax(k)
					  end do
				  end do

				  !update B fld near the axis
				  call UpdateBFldAxis
			  end if
		  end if
		  

     end subroutine UpdateBfieldHalf
	 

	 
     subroutine AddCurrent
		  real(psn) :: r,rmin,rp_half
          integer:: i,j,k,i1    
		  
		  
		  i1=3
 		  if(procxind.eq.0) i1=4
		  rmin=xborders(procxind)+rshift
#ifndef twoD
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                  do i=i1,mx-3
					     r=((i-3.0_psn)+rmin)*grid_dx
						 rp_half=r+0.5_psn*grid_dx
                         Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)*grid_inv_dz/(rp_half*dtheta)
                         Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)*grid_inv_dx*grid_inv_dz
                         Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)*grid_inv_dx/(r*dtheta)
                    end do
               end do
          end do
		  
		  if(procxind.eq.0) call UpdateEFldAxis
     end subroutine AddCurrent
	 
	 
	 subroutine UpdateEFldAxis
		 integer :: j,k,jp,jm,jpp
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) then
			 
#ifdef twoD
         	do k=1,1 
#else 
         	do k=1,mz 
#endif 	
            	do j=1,my
					
					! using the grid points across and around the axis ddd
					jp=j+ny/4
					jm=j-ny/4
					jpp = j + ny/2
					if(jp.gt.my-3) jp= jp - ny
					if(jpp.gt.my-3) jpp= jpp - ny
					if(jm.lt.3) jm = jm + ny
					
                    Ex(3,j,k) = 0.5_psn*(Ex(4,j,k) - Ex(4,jpp,k))
					!Ex(3,j,k) = 0.5_psn*(Ey(4,jm,k) - Ey(4,jp,k))
					Ex(2,j,k) = -Ex(4,jpp,k)

					Ey(3,j,k) = -Ey(4,jpp,k)
					
					Ez(3,j,k) = Ez(4,jpp,k)
					
					
! 		   			!setting the fld at the axis (as seen by the particles) to zero
! 					Ex(3,j,k)=0.0_psn
! 		   			Ex(2,j,k)=-Ex(4,j,k)
!
! 		   			Ey(3,j,k)=-Ey(4,j,k)
! 		   			Ez(3,j,k)=-Ez(4,j,k)
					
										
				end do
			end do
			

			 
		 end if 
	 end subroutine UpdateEFldAxis
	 
	 subroutine UpdateBFldAxis
		 integer :: j,k,jp,jm,jpp
		 
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) then
#ifdef twoD
         	do k=1,1 
#else 
         	do k=1,mz 
#endif 	
            	do j=1,my
					! using the grid points across and around the axis
					jp=j+ny/4
					jm=j-ny/4
					jpp = j + ny/2
					if(jp.gt.my-3) jp= jp - ny
					if(jpp.gt.my-3) jpp= jpp - ny
					if(jm.lt.3) jm = jm + ny
					
					By(3,j,k) = 0.5_psn*(By(4,j,k)-By(4,jpp,k))
					!By(3,j,k) = 0.5_psn*(Bx(4,jp,k) - Bx(4,jm,k))
					By(2,j,k) = -By(4,jpp,k)
					
					Bx(3,j,k) = -Bx(4,jpp,k)
					
					Bz(2,j,k) = Bz(4,jpp,k)
					
! 					! setting the fld at the axis (as seen by the particles) to zero (a special case of extrapolation)
! 	   			    Bx(3,j,k)=-Bx(4,j,k)
! 	   			    By(3,j,k)=0.0_psn
! 	   			    By(2,j,k)=-By(4,j,k)
! 	   			    !Bz(2,j,k)=2*Bz(3,j,k)-Bz(4,j,k)
! 	   			    Bz(2,j,k)=Bz(3,j,k)
										
				end do
			end do
			
			
			 
		 end if 
		 
	 end subroutine UpdateBFldAxis
	 
	 
	 
	 subroutine UpdateBzAxis(Bz_ax)
		 real(psn), dimension(mz) :: Bz_ax
		 integer :: j,k
		 
		 Bz_ax=0
#ifdef twoD
         do k=1,1 
#else 		 
		 do k=1,mz-1
#endif			 		 
			 do j=3,my-3 !First calcualte del Bz
				 Bz_ax(k)=Bz_ax(k)-fld_halfc*ax_perm_area*Ey(4,j,k) 
			 end do 
		     Bz_ax(k)=Bz(3,3,k)+Bz_ax(k)
		 end do 	
	 end subroutine UpdateBzAxis
	 
	 real(psn) function gradX(F,i,j,k)
	 real(psn), dimension(:,:,:), intent(in) :: F
	 integer, intent(in) :: i,j,k
	 gradX =  alpha_x*( F(i+1,j  ,k  ) - F(i  ,j,  k  )) + &
			  beta_xz*( F(i+1,j  ,k+1) - F(i  ,j  ,k+1)) + &
			  beta_xz*( F(i+1,j  ,k-1) - F(i  ,j  ,k-1)) + &
  end function gradX

  real(psn) function gradX1(F,i,j,k, rp, r)
  real(psn), dimension(:,:,:), intent(in) :: F
  real(psn) :: rp , r
  integer, intent(in) :: i,j,k
  gradX1 =  alpha_x*( rp*F(i+1,j  ,k  ) - r*F(i  ,j,  k  )) + &
		    beta_xz*( rp*F(i+1,j  ,k+1) - r*F(i  ,j  ,k+1)) + &
		    beta_xz*( rp*F(i+1,j  ,k-1) - r*F(i  ,j  ,k-1)) + &
   end function gradX1

  real(psn) function gradY(F,i,j,k)
	 real(psn), dimension(:,:,:), intent(in) :: F
	 integer, intent(in) :: i,j,k
	 gradY =  alpha_y*grid_inv_dy*( F(i  ,j+1,k  ) - F(i  ,j,  k  )) + &
			  beta_yz*grid_inv_dy*( F(i  ,j+1,k+1) - F(i  ,j  ,k+1)) + &
			  beta_yz*grid_inv_dy*( F(i  ,j+1,k-1) - F(i  ,j  ,k-1)) + &
  end function gradY

  real(psn) function gradZ(F,i,j,k)
	 real(psn), dimension(:,:,:), intent(in) :: F
	 integer, intent(in) :: i,j,k
	 gradZ =  alpha_z*( F(i  ,j  ,k+1) - F(i  ,j,  k  )) + &
			  beta_zx*( F(i+1,j  ,k+1) - F(i+1,j,  k  )) + &
			  beta_zx*( F(i-1,j  ,k+1) - F(i-1,j,  k  )) + &
			  beta_zy*( F(i  ,j+1,k+1) - F(i  ,j+1,k  )) + &
			  beta_zy*( F(i  ,j-1,k+1) - F(i  ,j-1,k  )) + &
			  delta_z*( F(i  ,j  ,k+2) - F(i  ,j  ,k-1)) 
  end function gradZ


end module em_update