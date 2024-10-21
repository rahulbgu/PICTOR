module cyl_filter
    use parameters
    use vars
	use cyl_vars
	use comm_fldprtl
	use cyl_comm_fldprtl
	use cyl_bc
	implicit none 
contains 
    
	subroutine smoothen_current_cyl
         integer :: ni
         do ni=1,curr_filter
              call MovingAverageFilter_cyl(Jx,0.5_psn)
              call MovingAverageFilter_cyl(Jy,0.0_psn)
              call MovingAverageFilter_cyl(Jz,0.0_psn)
         end do
    end subroutine smoothen_current_cyl
	
     subroutine MovingAverageFilter_cyl(J0,shift)
          real(psn),dimension(mx,my,mz) :: J0
          real(psn),dimension(mx,my,mz) ::FldTemp
          real(psn) :: shift
		  integer :: i,j,k,i1
		  real(psn) :: r,rmin, wt0, wt1
	
		  i1=3
 		  if(procxind.eq.0) i1=4
		  rmin=rborders(procxind)
		
		  FldTemp=0.0_psn
          call ExchangeZXEdgeCurrent1(J0)                    
#ifndef twoD
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                   do i=i1,mx-3
					     r=(i-3.0_psn)+rmin+shift
						 wt0=1.0_psn - 0.5_psn/(1.0_psn+r)
						 wt1= 0.25_psn/(1.0_psn+r) 
                         FldTemp(i,j,k)=wt1*J0(i,j-1,k)+wt0*J0(i,j,k)+wt1*J0(i,j+1,k)
                    end do
               end do
          end do
          J0=FldTemp      
     end subroutine MovingAverageFilter_cyl
	 
end module cyl_filter