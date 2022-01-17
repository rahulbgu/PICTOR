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
	 		  if(inc_axis) then 
				  call AxisCurrentBC
			  else 
				  call FoldInCurrentRmin
			  end if 
	 		  call CurrentBC_Outer
         end do
    end subroutine smoothen_current_cyl
	
     subroutine MovingAverageFilter_cyl(J0,dx)
          real(psn),dimension(mx,my,mz) :: J0
          real(psn),dimension(mx,my,mz) ::FldTemp
          real(psn) :: fm,fp,dx
		  integer :: i,j,k
	

		
		FldTemp=0.0_psn
		call ExchangeYZEdgeCurrent1(J0)
		if(inc_axis) call ExchangeYZEdgeCurrent1_Axis(J0)  		
#ifndef twoD
         do k=3,mz-3
#else
        do k=1,1
#endif
              do j=3,my-3
                  do i=3,mx-3
					     !fp=(i-3.0_psn+rborders(procxind(proc))+dx)/(i-2.0_psn+rborders(procxind(proc))+dx)
						 !fm=(i-3.0_psn+rborders(procxind(proc))+dx)/(i-4.0_psn+rborders(procxind(proc))+dx)
						 !if((procxind(proc).eq.0).and.i.eq.4) fm=0
                         fp=1.0_psn
						 fm=1.0_psn
						 FldTemp(i,j,k)=wtm1*J0(i-1,j,k)*fm +wt0*J0(i,j,k)+ wtp1*J0(i+1,j,k)*fp
                  end do
               end do
          end do
          J0=FldTemp
		  
          call ExchangeZXEdgeCurrent1(J0)                    
#ifndef twoD
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                   do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*J0(i,j-1,k)+wt0*J0(i,j,k)+wtp1*J0(i,j+1,k)
                    end do
               end do
          end do
          J0=FldTemp
#ifndef twoD          

          call ExchangeXYEdgeCurrent1(J0)
          do k=3,mz-3
               do j=3,my-3
                    do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*J0(i,j,k-1)+wt0*J0(i,j,k)+wtp1*J0(i,j,k+1)
                    end do
               end do
          end do
          J0=FldTemp
#endif        
     end subroutine MovingAverageFilter_cyl
	 
	 subroutine FoldInCurrentRmax
		 integer :: ind
		 real(psn) :: f
		 ind=BC_Rmax_Prtl-rborders(procxind(proc))+xmin
		 !f=(BC_Rmax_Prtl+0.5_psn)/(BC_Rmax_Prtl-0.5_psn) 
 		 f=(BC_Rmax_Prtl-0.5_psn)/(BC_Rmax_Prtl+0.5_psn) 
		 if(ind.ge.2.and.ind.le.mx-2) then
 		 	Jy(ind,:,:)=Jy(ind,:,:)+Jy(ind+1,:,:)*f 
 		 	Jz(ind,:,:)=Jz(ind,:,:)+Jz(ind+1,:,:)*f
 			Jx(ind,:,:)=0.0_psn
 			Jy(ind+1,:,:)=0.0_psn
 			Jz(ind+1,:,:)=0.0_psn
 		 end if	 
	 end subroutine FoldInCurrentRmax
	 
	 subroutine FoldInCurrentRmin
		integer :: ind
		ind=BC_Rmin_Prtl-rborders(procxind(proc))+xmin
 		if(ind.ge.3.and.ind.le.mx-2) then
 			!To ensure that the current is deposited on right place for reflected particles
 			Jx(ind,:,:)=Jx(ind,:,:)-Jx(ind-1,:,:)
 		 	Jy(ind,:,:)=Jy(ind,:,:)+Jy(ind-1,:,:)
 		 	Jz(ind,:,:)=Jz(ind,:,:)+Jz(ind-1,:,:)
 			Jx(ind-1,:,:)=0.0_psn
 			Jy(ind-1,:,:)=0.0_psn
 			Jz(ind-1,:,:)=0.0_psn
 		end if
	 end subroutine FoldInCurrentRmin 
	 
end module cyl_filter