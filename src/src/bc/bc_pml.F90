module bc_pml
    use parameters
    use vars
    implicit none

    type pml_params
        real(psn) :: kappa_max = 1.0_psn
        real(psn) :: sigma_max = 0.0_psn
        real(psn) :: alpha_max = 1.0_psn
        integer   :: n = 4
        real(psn) :: L = 10.0_psn ! thickness
    end type pml_params
    logical :: no_pml_this_domain = .true.
    type(pml_params), dimension(6) :: pml_params_arr
	real(dbpsn), dimension(6) :: pml_pos	= (/ real(-huge(1),dbpsn), real(huge(1),dbpsn), real(-huge(1),dbpsn), real(huge(1),dbpsn), real(-huge(1),dbpsn), real(huge(1),dbpsn) /)
    integer, dimension(6) :: pml_pos_ind = (/ -huge(1), huge(1), -huge(1), huge(1), -huge(1), huge(1) /)

contains

subroutine InitPMLArr
    allocate(pml_E_psi1%x(mx,my,mz), pml_E_psi1%y(mx,my,mz) , pml_E_psi1%z(mx,my,mz))
    allocate(pml_E_psi2%x(mx,my,mz), pml_E_psi2%y(mx,my,mz) , pml_E_psi2%z(mx,my,mz))
    pml_E_psi1%x = 0; pml_E_psi1%y = 0; pml_E_psi1%z = 0; 
    pml_E_psi2%x = 0; pml_E_psi2%y = 0; pml_E_psi2%z = 0; 

    allocate(pml_B_psi1%x(mx,my,mz), pml_B_psi1%y(mx,my,mz) , pml_B_psi1%z(mx,my,mz))
    allocate(pml_B_psi2%x(mx,my,mz), pml_B_psi2%y(mx,my,mz) , pml_B_psi2%z(mx,my,mz))
    pml_B_psi1%x = 0; pml_B_psi1%y = 0; pml_B_psi1%z = 0; 
    pml_B_psi2%x = 0; pml_B_psi2%y = 0; pml_B_psi2%z = 0; 
end subroutine InitPMLArr

subroutine DeallocatePMLArr
    deallocate(pml_E_psi1%x, pml_E_psi1%y , pml_E_psi1%z)
    deallocate(pml_E_psi2%x, pml_E_psi2%y , pml_E_psi2%z)
    deallocate(pml_B_psi1%x, pml_B_psi1%y , pml_B_psi1%z)
    deallocate(pml_B_psi2%x, pml_B_psi2%y , pml_B_psi2%z)
end subroutine DeallocatePMLArr

subroutine SetPMLParameters(Side,Thickness)
    integer   :: Side
    real(psn) :: Thickness

    pmlBC = .true.

    pml_params_arr(side)%n = 4 !default
    pml_params_arr(side)%L = Thickness
    pml_params_arr(side)%kappa_max = 1.0_psn
    pml_params_arr(side)%sigma_max = 1.0!0.75*(pml_params_arr(side)%n+1)/(150.0*pi*grid_dx)
    pml_params_arr(side)%alpha_max = 9.0e-4

end subroutine SetPMLParameters


subroutine Update_PML_Position
    integer :: n

    !update position of the attenuating boundaries
    do n=1,6
        if(bc_face(n)%type_fld.eq.'pml') then 
            pml_pos(n) = bc_face(n)%pos_fld
            if(mod(n,2).eq.0) pml_pos_ind(n) = ceiling(pml_pos(n)) + 1
            if(mod(n,2).ne.0) pml_pos_ind(n) = floor(pml_pos(n)) - 1
        end if
    end do

    no_pml_this_domain =  1+xborders(procxind)-3 .gt. pml_pos_ind(1) .and. mx+xborders(procxind)-3 .lt. pml_pos_ind(2) .and. & 
                          1+yborders(procyind)-3 .gt. pml_pos_ind(3) .and. my+yborders(procyind)-3 .lt. pml_pos_ind(4) .and. & 
                          1+zborders(proczind)-3 .gt. pml_pos_ind(5) .and. mz+zborders(proczind)-3 .lt. pml_pos_ind(6) 

end subroutine Update_PML_Position
 

logical function pml_outside_domain(mx,my,mz,sx,sy,sz,xborders,yborders,zborders)
    integer       :: mx,my,mz
    integer       :: sx,sy,sz
    integer, dimension(0:sx) :: xborders
    integer, dimension(0:sy) :: yborders
    integer, dimension(0:sz) :: zborders

    pml_outside_domain =  1+xborders(procxind)-3 .gt. pml_pos_ind(1) .and. mx+xborders(procxind)-3 .lt. pml_pos_ind(2) .and. & 
                          1+yborders(procyind)-3 .gt. pml_pos_ind(3) .and. my+yborders(procyind)-3 .lt. pml_pos_ind(4) .and. & 
                          1+zborders(proczind)-3 .gt. pml_pos_ind(5) .and. mz+zborders(proczind)-3 .lt. pml_pos_ind(6) 

end function pml_outside_domain

subroutine UpdateElcFld_PML
    integer :: i,j,k 
    logical :: inside, outside
    real(dbpsn) :: xg,yg,zg
    real(psn)   :: a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z

    ! the entire subdomain is in the interior and does not need to updated any fld
    if(no_pml_this_domain) return 
    
#ifndef twoD
    do k=2,mz-1
#else
    do k=1,1
#endif
        do j=2,my-1
            do i=2,mx-1

                xg = i + xborders(procxind)-3
                yg = j + yborders(procyind)-3 
                zg = k + zborders(proczind)-3
                
                !check if the grid point is in the interior, and not in the PML
                inside =  xg .gt. pml_pos_ind(1) .and. xg .lt. pml_pos_ind(2) .and. &
                          yg .gt. pml_pos_ind(3) .and. yg .lt. pml_pos_ind(4) .and. &
                          zg .gt. pml_pos_ind(5) .and. zg .lt. pml_pos_ind(6)
                
                if(inside) cycle



                !calculate all coeff. 
                call pml_ab(xg,0,a_x,b_x)
                call pml_ab(yg,1,a_y,b_y)
                call pml_ab(zg,2,a_z,b_z)
                c_x = pml_c(xg,0)
                c_y = pml_c(yg,1)
                c_z = pml_c(zg,2)

                !update psi
                pml_E_psi1%x(i,j,k) = a_y*pml_E_psi1%x(i,j,k) + b_y*fldc*grid_inv_dy*(Bz(i,j,k)-Bz(i,j-1,k))
                pml_E_psi1%z(i,j,k) = a_x*pml_E_psi1%z(i,j,k) + b_x*fldc*grid_inv_dx*(By(i,j,k)-By(i-1,j,k))
#ifndef twoD
                pml_E_psi1%y(i,j,k) = a_z*pml_E_psi1%y(i,j,k) + b_z*fldc*grid_inv_dz*(Bx(i,j,k)-Bx(i,j,k-1))                
                pml_E_psi2%x(i,j,k) = a_z*pml_E_psi2%x(i,j,k) + b_z*fldc*grid_inv_dz*(By(i,j,k)-By(i,j,k-1))
#endif
                pml_E_psi2%y(i,j,k) = a_x*pml_E_psi2%y(i,j,k) + b_x*fldc*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k))
                pml_E_psi2%z(i,j,k) = a_y*pml_E_psi2%z(i,j,k) + b_y*fldc*grid_inv_dy*(Bx(i,j,k)-Bx(i,j-1,k))

                
#ifndef twoD
                Ex(i,j,k)=Ex(i,j,k) + c_y*grid_inv_dy*(Bz(i,j,k)-Bz(i,j-1,k)) &
                                    - c_z*grid_inv_dz*(By(i,j,k)-By(i,j,k-1)) &
                                    + pml_E_psi1%x(i,j,k) - pml_E_psi2%x(i,j,k)

                Ey(i,j,k)=Ey(i,j,k) + c_z*grid_inv_dz*(Bx(i,j,k)-Bx(i,j,k-1)) &
                                    - c_x*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k)) &
                                    + pml_E_psi1%y(i,j,k) - pml_E_psi2%y(i,j,k)

#else
                Ex(i,j,k)=Ex(i,j,k) + c_y*grid_inv_dy*(Bz(i,j,k)-Bz(i,j-1,k)) &
                                    + pml_E_psi1%x(i,j,k)
                Ey(i,j,k)=Ey(i,j,k) - c_x*grid_inv_dx*(Bz(i,j,k)-Bz(i-1,j,k)) &
                                    - pml_E_psi2%y(i,j,k)
#endif

                Ez(i,j,k)=Ez(i,j,k) + c_x*grid_inv_dx*(By(i,j,k)-By(i-1,j,k)) &
                                    - c_y*grid_inv_dy*(Bx(i,j,k)-Bx(i,j-1,k)) &
                                    + pml_E_psi1%z(i,j,k) - pml_E_psi2%z(i,j,k)
                
                outside =  xg + 0.5_psn .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg + 0.5_psn .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg           .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg           .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg           .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg           .gt. pml_pos(6) + pml_params_arr(6)%L
                
                if(outside) Ex(i,j,k) = 0.0_psn

                outside =  xg           .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg           .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg + 0.5_psn .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg + 0.5_psn .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg           .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg           .gt. pml_pos(6) + pml_params_arr(6)%L
     
                if(outside) Ey(i,j,k) = 0.0_psn

                outside =  xg           .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg           .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg           .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg           .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg + 0.5_psn .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg + 0.5_psn .gt. pml_pos(6) + pml_params_arr(6)%L
     
                if(outside) Ez(i,j,k) = 0.0_psn

            end do
        end do
    end do
		

end subroutine UpdateElcFld_PML

subroutine UpdateMagFld_PML
    integer :: i,j,k 
    logical :: inside, outside
    real(dbpsn) :: xg,yg,zg
    real(psn)   :: a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z
    

    call Update_PML_Position

    !the entire subdomain is in the interior and does not need to updated any fld
    if(no_pml_this_domain) then
        if (allocated(pml_E_psi1%x)) call DeallocatePMLArr
        return
    else
        if(.not.allocated(pml_E_psi1%x)) call InitPMLArr
    end if 

    
#ifndef twoD
    do k=2,mz-1
#else
    do k=1,1
#endif
        do j=2,my-1
            do i=2,mx-1
                
                !check if the grid point is in the interior, and not in the PML
                inside =  i+xborders(procxind)-3 .gt. pml_pos_ind(1) .and. i+xborders(procxind)-3 .lt. pml_pos_ind(2) .and. & 
                          j+yborders(procyind)-3 .gt. pml_pos_ind(3) .and. j+yborders(procyind)-3 .lt. pml_pos_ind(4) .and. & 
                          k+zborders(proczind)-3 .gt. pml_pos_ind(5) .and. k+zborders(proczind)-3 .lt. pml_pos_ind(6)  
                
                if(inside) cycle

                xg = i + xborders(procxind)-3 + 0.5_dbpsn
                yg = j + yborders(procyind)-3 + 0.5_dbpsn
                zg = k + zborders(proczind)-3 + 0.5_dbpsn
                
                !calculate all coeff. 
                call pml_ab(xg,0,a_x,b_x)
                call pml_ab(yg,1,a_y,b_y)
                call pml_ab(zg,2,a_z,b_z)
                c_x = pml_c(xg,0)
                c_y = pml_c(yg,1)
                c_z = pml_c(zg,2)

!                 !update psi
                pml_B_psi1%x(i,j,k) = a_y*pml_B_psi1%x(i,j,k) + b_y*fldc*grid_inv_dy*(Ez(i,j+1,k)-Ez(i,j,k)) 
                pml_B_psi1%z(i,j,k) = a_x*pml_B_psi1%z(i,j,k) + b_x*fldc*grid_inv_dx*(Ey(i+1,j,k)-Ey(i,j,k))
#ifndef twoD
                pml_B_psi1%y(i,j,k) = a_z*pml_B_psi1%y(i,j,k) + b_z*fldc*grid_inv_dz*(Ex(i,j,k+1)-Ex(i,j,k))                
                pml_B_psi2%x(i,j,k) = a_z*pml_B_psi2%x(i,j,k) + b_z*fldc*grid_inv_dz*(Ey(i,j,k+1)-Ey(i,j,k))
#endif
                pml_B_psi2%y(i,j,k) = a_x*pml_B_psi2%y(i,j,k) + b_x*fldc*grid_inv_dx*(Ez(i+1,j,k)-Ez(i,j,k))
                pml_B_psi2%z(i,j,k) = a_y*pml_B_psi2%z(i,j,k) + b_y*fldc*grid_inv_dy*(Ex(i,j+1,k)-Ex(i,j,k))

                
#ifndef twoD
                Bx(i,j,k)=Bx(i,j,k) - c_y*grid_inv_dy*(Ez(i,j+1,k)-Ez(i,j,k)) &
                                    + c_z*grid_inv_dz*(Ey(i,j,k+1)-Ey(i,j,k)) &
                                    - pml_B_psi1%x(i,j,k) + pml_B_psi2%x(i,j,k)

                By(i,j,k)=By(i,j,k) - c_z*grid_inv_dz*(Ex(i,j,k+1)-Ex(i,j,k)) &
                                    + c_x*grid_inv_dx*(Ez(i+1,j,k)-Ez(i,j,k)) &
                                    - pml_B_psi1%y(i,j,k) + pml_B_psi2%y(i,j,k)

#else
                Bx(i,j,k)=Bx(i,j,k) - c_y*grid_inv_dy*(Ez(i,j+1,k)-Ez(i,j,k)) &
                                    - pml_B_psi1%x(i,j,k) 
                By(i,j,k)=By(i,j,k) + c_x*grid_inv_dx*(Ez(i+1,j,k)-Ez(i,j,k)) &
                                    + pml_B_psi2%y(i,j,k)

#endif
                Bz(i,j,k)=Bz(i,j,k) - c_x*grid_inv_dx*(Ey(i+1,j,k)-Ey(i,j,k)) &
                                    + c_y*grid_inv_dy*(Ex(i,j+1,k)-Ex(i,j,k)) &
                                    - pml_B_psi1%z(i,j,k) + pml_B_psi2%z(i,j,k)
                
                outside =  xg - 0.5_psn .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg - 0.5_psn .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg           .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg           .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg           .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg           .gt. pml_pos(6) + pml_params_arr(6)%L
                
                if(outside) Bx(i,j,k) = 0.0_psn

                outside =  xg           .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg           .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg - 0.5_psn .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg - 0.5_psn .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg           .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg           .gt. pml_pos(6) + pml_params_arr(6)%L
    
                if(outside) By(i,j,k) = 0.0_psn

                outside =  xg           .lt. pml_pos(1) - pml_params_arr(1)%L .or. xg           .gt. pml_pos(2) + pml_params_arr(2)%L .or. &
                           yg           .lt. pml_pos(3) - pml_params_arr(3)%L .or. yg           .gt. pml_pos(4) + pml_params_arr(4)%L .or. &
                           zg - 0.5_psn .lt. pml_pos(5) - pml_params_arr(5)%L .or. zg - 0.5_psn .gt. pml_pos(6) + pml_params_arr(6)%L
    
                if(outside) Bz(i,j,k) = 0.0_psn


            end do
        end do
    end do
		

end subroutine UpdateMagFld_PML




real(psn) function pml_c(x,axis)
    integer, intent(in)     :: axis
    real(dbpsn), intent(in) :: x
    real(psn)               :: kappa, d
    
    pml_c = 0.0_psn

    if(x .gt. pml_pos(2*axis+2)) then
        d = x - pml_pos(2*axis+2)
        kappa = pml_kappa(d,2*axis+2)
        pml_c = fldc*(1.0_psn/kappa - 1.0_psn ) 
    end if 
    if(x.lt.pml_pos(2*axis+1)) then
        d = pml_pos(2*axis+1) - x 
        kappa = pml_kappa(d,2*axis+1)
        pml_c = fldc*(1.0_psn/kappa - 1.0_psn ) 
    end if
    
end function pml_c

subroutine pml_ab(x,axis, a , b)
    integer, intent(in)     :: axis
    real(dbpsn), intent(in) :: x
    real(psn),intent(inout) :: a, b
    real(psn)               :: alpha, sigma, kappa, d
    
    a = 0.0_psn
    b = 0.0_psn

    if(x .gt. pml_pos(2*axis+2)) then
        d = x - pml_pos(2*axis+2) 
        alpha = pml_alpha(d,2*axis+2)
        sigma = pml_sigma(d,2*axis+2)
        kappa = pml_kappa(d,2*axis+2)
        a = exp(-(sigma/kappa)-alpha)
        b = ( sigma / (sigma*kappa + alpha*kappa**2) ) * (a - 1.0_psn)
    end if 
    if(x.lt.pml_pos(2*axis+1)) then
        d = pml_pos(2*axis+1) - x 
        alpha = pml_alpha(d,2*axis+1)
        sigma = pml_sigma(d,2*axis+1)
        kappa = pml_kappa(d,2*axis+1)
        a = exp(-(sigma/kappa)-alpha)
        b = ( sigma / (sigma*kappa + alpha*kappa**2) ) * (a - 1.0_psn)
    end if 
    
end subroutine pml_ab

real(psn) function pml_kappa(d,side)
    real(psn), intent(in) :: d
    integer, intent(in)   :: side
    pml_kappa = min( pml_params_arr(side)%kappa_max , 1.0_psn + ((d/pml_params_arr(side)%L)**pml_params_arr(side)%n)*(pml_params_arr(side)%kappa_max - 1.0_psn)  )
end function pml_kappa 

real(psn) function pml_sigma(d,side)
    real(psn), intent(in) :: d
    integer, intent(in)   :: side
    pml_sigma = min( pml_params_arr(side)%sigma_max , pml_params_arr(side)%sigma_max*(d/pml_params_arr(side)%L)**pml_params_arr(side)%n )
end function pml_sigma

real(psn) function pml_alpha(d,side)
    real(psn), intent(in) :: d
    integer, intent(in)   :: side
    pml_alpha = pml_params_arr(side)%alpha_max
end function pml_alpha




end module bc_pml 