module em_update
    use parameters
    use vars
	use cyl_vars
    implicit none
contains 
	
     subroutine UpdateEfield
		 real(psn) :: rmin,r,rp_half,rm_half
		 integer :: i,j,k
         integer :: i1,k1,k2
	
		 i1=3
		 if(procxind(proc).eq.0) i1=4

		 rmin=rborders(procxind(proc))
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
					r=(i-3.0_psn)+rmin
					rp_half=r+0.5_psn
					rm_half=r-0.5_psn
#ifndef twoD
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))/(rp_half*dtheta) - fldc*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(rp_half*By(i,j,k)-rm_half*By(i-1,j,k))/r - fldc*(Bx(i,j,k)-Bx(i,j-1,k))/(r*dtheta)
#else
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))/(rp_half*dtheta)
                    Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*( (rp_half*By(i,j,k)-rm_half*By(i-1,j,k)) - (Bx(i,j,k)-Bx(i,j-1,k))/dtheta  )/r
#endif
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
	
 		  i1=1
 		  if(procxind(proc).eq.0) i1=4

		  rmin=rborders(procxind(proc))
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
						r=(i-3.0_psn)+rmin
						rp_half=r+0.5_psn
						rp=r+1.0_psn
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))/(r*dtheta) + fld_halfc*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*( (rp*Ey(i+1,j,k)-r*Ey(i,j,k)) - (Ex(i,j+1,k)-Ex(i,j,k))/dtheta )/rp_half
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))/(r*dtheta)
                         By(i,j,k)=By(i,j,k)+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
						 Bz(i,j,k)=Bz(i,j,k)-fld_halfc*( (rp*Ey(i+1,j,k)-r*Ey(i,j,k)) - (Ex(i,j+1,k)-Ex(i,j,k))/dtheta )/rp_half
#endif
                    end do
               end do
          end do

		  if(inc_axis) then
			  if(procxind(proc).eq.0) then
		          
				  !update By at the axis using the standard FDTD equation
				  i=3
				  do k=k1,k2
		               do j=1,my-1
						    By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
					   end do 
				  end do 
				  
				  !update Bz using the integral equation
				  call UpdateBzAxis
				  do k=k1,k2
					  do j=1,my-1
						   Bz(3,j,k)=Bz_ax(k)
					  end do
				  end do
				  
				  !update B fld near the axis
				  call UpdateBFldAxis
			  end if
		  end if
		  
		  ! Temp BC, extrapolate and B a the ext boundary to zero 
		  !call Bfld_exterior

     end subroutine UpdateBfieldHalf
	 
	 subroutine Bfld_exterior
		 integer   :: rmax_local
 		 
		 rmax_local=BC_Rmax_Fld-rborders(procxind(proc))+3
 		 rmax_local=max(rmax_local,1)
		 
 		if(rmax_local.le.mx) then 
 			By(rmax_local,:,:)  = -By(rmax_local-1,:,:)
			By(rmax_local+1,:,:)= -By(rmax_local-2,:,:)
			
 			Bx(rmax_local,:,:)  = 0.0_psn
			Bx(rmax_local+1,:,:)= 0.0_psn
			
 		end if 
	 
	 
	 end subroutine Bfld_exterior
	 
     subroutine AddCurrent
		  real(psn) :: r,rmin,rp_half
          integer:: i,j,k,i1    
		                
 	
		  i1=3
 		  if(procxind(proc).eq.0) i1=4
		  rmin=rborders(procxind(proc))
#ifndef twoD
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                  do i=i1,mx-3
					     r=(i-3.0_psn)+rmin
						 rp_half=r+0.5_psn
                         Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)/(rp_half*dtheta)
                         Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)
                         Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)/(r*dtheta)
                    end do
               end do
          end do
		  
		  if(procxind(proc).eq.0) call UpdateEFldAxis
     end subroutine AddCurrent
	 
	 
	 subroutine UpdateEFldAxis
		 integer :: j,k,jp,jm,jpp
		 if((inc_axis.eqv..true.).and.(procxind(proc).eq.0)) then
			 
#ifdef twoD
         	do k=1,1 
#else 
         	do k=1,mz 
#endif 	
            	do j=1,my
					
! 					! using the grid points across and around the axis, similar to the papers in IEEE
! 					jp=j+ny/4
! 					jm=j-ny/4
! 					jpp = j + ny/2
! 					if(jp.gt.my-3) jp= jp - ny
! 					if(jpp.gt.my-3) jpp= jpp - ny
! 					if(jm.lt.3) jm = jm + ny
!
! 					Ex(3,j,k) = 0.5_psn*(Ey(4,jm,k) - Ey(4,jp,k))
! 					Ex(2,j,k) = -Ex(4,jpp,k)
!
! 					Ey(3,j,k) = -Ey(4,jpp,k)
! 					Ez(3,j,k) = Ez(4,jpp,k)
					
		   			! setting the fld at the axis (as seen by the particles) to zero 
					Ex(3,j,k)=0.0_psn
		   			Ex(2,j,k)=-Ex(4,j,k)

		   			Ey(3,j,k)=-Ey(4,j,k)
		   			Ez(3,j,k)=-Ez(4,j,k)
					
										
				end do
			end do
			

			 
		 end if 
	 end subroutine UpdateEFldAxis
	 
	 subroutine UpdateBFldAxis

		 integer :: j,k,jpp
		 if((inc_axis.eqv..true.).and.(procxind(proc).eq.0)) then
#ifdef twoD
         	do k=1,1 
#else 
         	do k=1,mz 
#endif 	
            	do j=1,my
! 					! using the grid points across and around the axis, similar to the papers in IEEE
! 					jpp = j + ny/2
! 					if(jpp.gt.my-3) jpp= jpp - ny
!
! 					Bx(3,j,k) = -Bx(4,jpp,k)
! 					By(2,j,k) = -By(4,jpp,k)
! 					Bz(2,j,k) = Bz(4,jpp,k)
					
					! setting the fld at the axis (as seen by the particles) to zero (a special case of extrapolation)
	   			    Bx(3,j,k)=-Bx(4,j,k)
	   			    By(3,j,k)=0.0_psn
	   			    By(2,j,k)=-By(4,j,k)
	   			    !Bz(2,j,k)=2*Bz(3,j,k)-Bz(4,j,k)
	   			    Bz(2,j,k)=Bz(3,j,k)
										
				end do
			end do
			
			
			 
		 end if 
		 
	 end subroutine UpdateBFldAxis
	 
	 
	 
	 
	 subroutine UpdateBzAxis
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
	 
	  
	 

end module em_update