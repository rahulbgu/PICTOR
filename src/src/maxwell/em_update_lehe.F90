module em_update
    use parameters
    use vars
    use comm_fld
    
    implicit none
    
    real(psn), parameter :: beta_xy = 0.125_psn*(grid_dx/grid_dy)**2
    real(psn), parameter :: beta_xz = 0.125_psn*(grid_dx/grid_dz)**2
    real(psn), parameter :: beta_yx = 0.125_psn
    real(psn), parameter :: beta_yz = 0.0_psn
    real(psn), parameter :: beta_zx = 0.125_psn
    real(psn), parameter :: beta_zy = 0.0_psn
    real(psn), parameter :: delta_x = 0.25_psn*( 1.0_psn - (cinv**2) * (sin(0.5_psn*pi*c))**2 ) 
    real(psn), parameter :: delta_y = 0 
    real(psn), parameter :: delta_z = 0
#ifndef twoD    
    real(psn), parameter :: alpha_x = 1.0_psn - 2.0_psn*beta_xy - 2.0_psn*beta_xz - 3.0_psn*delta_x 
    real(psn), parameter :: alpha_y = 1.0_psn - 2.0_psn*beta_yx - 2.0_psn*beta_yz - 3.0_psn*delta_y
    real(psn), parameter :: alpha_z = 1.0_psn - 2.0_psn*beta_zx - 2.0_psn*beta_zy - 3.0_psn*delta_z
#else
    real(psn), parameter :: alpha_x = 1.0_psn - 2.0_psn*beta_xy  - 3.0_psn*delta_x 
    real(psn), parameter :: alpha_y = 1.0_psn - 2.0_psn*beta_yx  - 3.0_psn*delta_y
    real(psn), parameter :: alpha_z = 1.0_psn 
#endif    

    private :: gradX, gradY, gradZ

contains 
    
    ! subroutine Set_FDTD_Constants
    !     beta_yx = 0.125_psn
    !     beta_zx = 0.125_psn
    !     beta_zy = 0.0_psn
    !     beta_yz = 0.0_psn
    !     beta_xy = 0.125_psn*(grid_dx/grid_dy)**2
    !     beta_xz = 0.125_psn*(grid_dx/grid_dz)**2
    !     delta_y = 0
    !     delta_z = 0
    !     delta_x = 0.25_psn*( 1.0_psn - (cinv**2) * (sin(0.5_psn*pi*c))**2 ) 
    !  end subroutine Set_FDTD_Constants 
	
     subroutine UpdateEfield
		 integer :: i,j,k
        

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j,k)
#endif				
			do k=3,mz-3
#else 
            do k=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j)
#endif
#endif 
             do j=3,my-3
!$OMP SIMD				 
                  do i=3,mx-3
! #ifndef twoD
!                     Ex(i,j,k)=Ex(i,j,k)+fldc*gradY(Bz,i,j-1,k) -fldc*gradZ(By,i,j,k-1)
!                     Ey(i,j,k)=Ey(i,j,k)+fldc*gradZ(Bx,i,j,k-1) -fldc*gradX(Bz,i-1,j,k)
!                     Ez(i,j,k)=Ez(i,j,k)+fldc*gradX(By,i-1,j,k) -fldc*gradY(Bx,i,j-1,k) 
! #else
!                     Ex(i,j,k)=Ex(i,j,k)+fldc*gradY(Bz,i,j-1,k) 
!                     Ey(i,j,k)=Ey(i,j,k)-fldc*gradX(Bz,i-1,j,k)
!                     Ez(i,j,k)=Ez(i,j,k)+fldc*gradX(By,i-1,j,k) -fldc*gradY(Bx,i,j-1,k) 
! #endif

#ifndef twoD
                    Ex(i,j,k)=Ex(i,j,k)+fldc*grid_inv_dy*(Bz(i,j,k)-Bz(i,j-1,k))-fldc*grid_inv_dz*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*grid_inv_dz*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*grid_inv_dx*(By(i,j,k)-By(i-1,j,k))-fldc*grid_inv_dy*(Bx(i,j,k)-Bx(i,j-1,k))
#else
                    Ex(i,j,k)=Ex(i,j,k)+fldc*grid_inv_dy*(Bz(i,j,k)-Bz(i,j-1,k))
                    Ey(i,j,k)=Ey(i,j,k)-fldc*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*grid_inv_dx*(By(i,j,k)-By(i-1,j,k))-fldc*grid_inv_dy*(Bx(i,j,k)-Bx(i,j-1,k))
#endif

                    end do 
               end do
          end do
     end subroutine UpdateEfield 


	 
     subroutine UpdateBfieldHalf
		  integer ::i,j,k 
		     

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j,k)
#endif				
			do k=3,mz-3
#else 
            do k=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j)
#endif
#endif
               do j=3,my-3
!$OMP SIMD				   
                    do i=3,mx-3
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*gradY(Ez,i,j,k)+fld_halfc*gradZ(Ey,i,j,k) 
                         By(i,j,k)=By(i,j,k)-fld_halfc*gradZ(Ex,i,j,k)+fld_halfc*gradX(Ez,i,j,k)
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*gradX(Ey,i,j,k)+fld_halfc*gradY(Ex,i,j,k)
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*gradY(Ez,i,j,k)
                         By(i,j,k)=By(i,j,k)+fld_halfc*gradX(Ez,i,j,k) 
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*gradX(Ey,i,j,k) +fld_halfc*gradY(Ex,i,j,k) 
#endif
                    end do
               end do
          end do 
     end subroutine UpdateBfieldHalf

     real(psn) function gradX(F,i,j,k)
        real(psn), dimension(:,:,:), intent(in) :: F
        integer, intent(in) :: i,j,k
        gradX =  alpha_x*grid_inv_dx*( F(i+1,j  ,k  ) - F(i  ,j,  k  )) + &
                 beta_xy*grid_inv_dx*( F(i+1,j+1,k  ) - F(i  ,j+1,k  )) + &
                 beta_xy*grid_inv_dx*( F(i+1,j-1,k  ) - F(i  ,j-1,k  )) + &
#ifndef twoD                 
                 beta_xz*grid_inv_dx*( F(i+1,j  ,k+1) - F(i  ,j  ,k+1)) + &
                 beta_xz*grid_inv_dx*( F(i+1,j  ,k-1) - F(i  ,j  ,k-1)) + &
#endif
                 delta_x*grid_inv_dx*( F(i+2,j  ,k  ) - F(i-1,j  ,k  )) 
     end function gradX

     real(psn) function gradY(F,i,j,k)
        real(psn), dimension(:,:,:), intent(in) :: F
        integer, intent(in) :: i,j,k
        gradY =  alpha_y*grid_inv_dy*( F(i  ,j+1,k  ) - F(i  ,j,  k  )) + &
                 beta_yx*grid_inv_dy*( F(i+1,j+1,k  ) - F(i+1,j,  k  )) + &
                 beta_yx*grid_inv_dy*( F(i-1,j+1,k  ) - F(i-1,j,  k  )) + &
#ifndef twoD
                 beta_yz*grid_inv_dy*( F(i  ,j+1,k+1) - F(i  ,j  ,k+1)) + &
                 beta_yz*grid_inv_dy*( F(i  ,j+1,k-1) - F(i  ,j  ,k-1)) + &
#endif
                 delta_y*grid_inv_dy*( F(i  ,j+2,k  ) - F(i  ,j-1,k  )) 
     end function gradY

     real(psn) function gradZ(F,i,j,k)
        real(psn), dimension(:,:,:), intent(in) :: F
        integer, intent(in) :: i,j,k
        gradZ =  alpha_z*grid_inv_dz*( F(i  ,j  ,k+1) - F(i  ,j,  k  )) + &
                 beta_zx*grid_inv_dz*( F(i+1,j  ,k+1) - F(i+1,j,  k  )) + &
                 beta_zx*grid_inv_dz*( F(i-1,j  ,k+1) - F(i-1,j,  k  )) + &
                 beta_zy*grid_inv_dz*( F(i  ,j+1,k+1) - F(i  ,j+1,k  )) + &
                 beta_zy*grid_inv_dz*( F(i  ,j-1,k+1) - F(i  ,j-1,k  )) + &
                 delta_z*grid_inv_dz*( F(i  ,j  ,k+2) - F(i  ,j  ,k-1)) 
     end function gradZ

     


	 
     subroutine AddCurrent
          integer:: i,j,k                    

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j,k)
#endif				
		  do k=3,mz-3
#else 
          do k=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j)
#endif
#endif
            do j=3,my-3
!$OMP SIMD				
                  do i=3,mx-3
                        Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)*grid_inv_dy*grid_inv_dz
                        Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)*grid_inv_dx*grid_inv_dz
                        Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)*grid_inv_dy*grid_inv_dx
                    end do
               end do
          end do
     end subroutine AddCurrent
	 
end module em_update