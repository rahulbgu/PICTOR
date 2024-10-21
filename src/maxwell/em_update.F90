module em_update
    use parameters
    use vars
    use comm_fld
    use bc_curved_surf
    implicit none
contains 
	
     subroutine UpdateEfield
		 integer :: i,j,k
        

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j,k)
#endif				
		  do k=2,mz-1
#else 
            do k=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j)
#endif
#endif 
             do j=2,my-1
!$OMP SIMD				 
                  do i=2,mx-1
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
		  
		  if(curved_bc) then 
		      call UpdateBfieldHalf_CFDTD 
			  return 
		  end if       
          
   

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j,k)
#endif				
		  do k=2,mz-1
#else 
            do k=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,j)
#endif
#endif
               do j=2,my-1
!$OMP SIMD				   
                    do i=2,mx-1
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*grid_inv_dy*(Ez(i,j+1,k)-Ez(i,j,k))+fld_halfc*grid_inv_dz*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*grid_inv_dz*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*grid_inv_dx*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*grid_inv_dx*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*grid_inv_dy*(Ex(i,j+1,k)-Ex(i,j,k))
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*grid_inv_dy*(Ez(i,j+1,k)-Ez(i,j,k))
                         By(i,j,k)=By(i,j,k)+fld_halfc*grid_inv_dx*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*grid_inv_dx*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*grid_inv_dy*(Ex(i,j+1,k)-Ex(i,j,k))
#endif
                    end do
               end do
          end do 
     end subroutine UpdateBfieldHalf
	 
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