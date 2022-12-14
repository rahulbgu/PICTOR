module current
     use parameters
     use vars
     use comm_fld
 
     implicit none
	 
	 logical :: ext_current_present = .false. ! no external current by default 

	 procedure(vector_global), pointer :: J_ext => null()
contains 
     
     subroutine ResetCurrent
          Jx=0.0_psn
          Jy=0.0_psn
          Jz=0.0_psn
     end subroutine ResetCurrent
	 
!      subroutine UpdateCurrentsAllEdges
!           call ExchangeYZEdgeCurrent
!           call AddImportedCurrentYZ
!           call ExchangeZXEdgeCurrent
!           call AddImportedCurrentZX
! #ifndef twoD
!           call ExchangeXYEdgeCurrent
!           call AddImportedCurrentXY
! #endif
!      end subroutine UpdateCurrentsAllEdges

     
     subroutine smoothen_current
          integer :: ni

		  do ni=1,curr_filter
			  
			   call MovingAverageFilter(Jx,Jy,Jz)
			   
			   call CurrentFilter_BC
		  
		  end do

     end subroutine smoothen_current
!--------------------------------------------------------------------------------------------------     
! Additional Subroutines: used for various advanced purposes, not essential for basic PIC  
!--------------------------------------------------------------------------------------------------     
    subroutine SetFilteredEfield
          integer :: ni
          FilteredEx=Ex
          FilteredEy=Ey
          FilteredEz=Ez
          if(nMoverEMfilter.gt.0) then 
               do ni=1,nMoverEMfilter
				    
					call MovingAverageFilter(FilteredEx,FilteredEy,FilteredEz)
					
               end do 
			   
			   call SendRecvFlds(FilteredEx,FilteredEy,FilteredEz)
    
         end if
     end subroutine SetFilteredEfield
     
!--------------------------------------------------------------------------------------------------     
! Filters: used to smoothen grid quantitites
!--------------------------------------------------------------------------------------------------

     subroutine MovingAverageFilter(Fldx,Fldy,Fldz)
		  real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
		  integer :: k1,k2
		  k1 = 2; k2 = mz-2;
#ifdef twoD
		  k1 =1; k2=1;  
#endif		  		  
		  call SendRecvFldsYZ(Fldx,Fldy,Fldz)		  
		  call MovingAverageX(Fldx,Fldy,Fldz, 2,mx-2,2,my-2,k1,k2)
		  
		  call SendRecvFldsZX(Fldx,Fldy,Fldz)
		  call MovingAverageY(Fldx,Fldy,Fldz, 2,mx-2,2,my-2,k1,k2)
		  
#ifndef twoD
          call SendRecvFldsXY(Fldx,Fldy,Fldz)		  
		  call MovingAverageZ(Fldx,Fldy,Fldz, 2,mx-2,2,my-2,2,mz-2)
#endif		  
       
     end subroutine MovingAverageFilter
	 
	 subroutine MovingAverageX(Fldx,Fldy,Fldz, i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
		 integer :: i1,i2,j1,j2,k1,k2
		 
		 call MovingAverageX1(Fldx,i1,i2,j1,j2,k1,k2)
		 call MovingAverageX1(Fldy,i1,i2,j1,j2,k1,k2)
		 call MovingAverageX1(Fldz,i1,i2,j1,j2,k1,k2)
		 
	 end subroutine MovingAverageX 
	 
	 subroutine MovingAverageX1(Fld,i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fld
		 real(psn), dimension(mx,my,mz) :: Fld_filtered
		 integer :: i1,i2,j1,j2,k1,k2
		 integer :: i,j,k
		 
#ifndef twoD
         do k=k1,k2
#else
         do k=1,1
#endif
         	do j=j1,j2
            	do i=i1,i2
                      Fld_filtered(i,j,k)=wtm1*Fld(i-1,j,k)+wt0*Fld(i,j,k)+wtp1*Fld(i+1,j,k)
                end do
            end do
         end do
		 
		 Fld(i1:i2,j1:j2,k1:k2) = Fld_filtered(i1:i2,j1:j2,k1:k2)
		  
	 end subroutine MovingAverageX1
	 
	 subroutine MovingAverageY(Fldx,Fldy,Fldz, i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
		 integer :: i1,i2,j1,j2,k1,k2
		 
		 call MovingAverageY1(Fldx,i1,i2,j1,j2,k1,k2)
		 call MovingAverageY1(Fldy,i1,i2,j1,j2,k1,k2)
		 call MovingAverageY1(Fldz,i1,i2,j1,j2,k1,k2)
		 
	 end subroutine MovingAverageY 
	 
	 subroutine MovingAverageY1(Fld,i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fld
		 real(psn), dimension(mx,my,mz) :: Fld_filtered
		 integer :: i1,i2,j1,j2,k1,k2
		 integer :: i,j,k
		 
#ifndef twoD
         do k=k1,k2
#else
         do k=1,1
#endif
      	 	do j=j1,j2
            	do i=i1,i2
                   Fld_filtered(i,j,k)=wtm1*Fld(i,j-1,k)+wt0*Fld(i,j,k)+wtp1*Fld(i,j+1,k)
                end do
            end do
         end do
		 
		 Fld(i1:i2,j1:j2,k1:k2) = Fld_filtered(i1:i2,j1:j2,k1:k2)
	 end subroutine MovingAverageY1
	 
	 subroutine MovingAverageZ(Fldx,Fldy,Fldz, i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
		 integer :: i1,i2,j1,j2,k1,k2
		 
		 call MovingAverageZ1(Fldx,i1,i2,j1,j2,k1,k2)
		 call MovingAverageZ1(Fldy,i1,i2,j1,j2,k1,k2)
		 call MovingAverageZ1(Fldz,i1,i2,j1,j2,k1,k2)
		 
	 end subroutine MovingAverageZ 
	 
	 subroutine MovingAverageZ1(Fld,i1,i2,j1,j2,k1,k2)
		 real(psn), dimension(mx,my,mz) :: Fld
		 real(psn), dimension(mx,my,mz) :: Fld_filtered
		 integer :: i1,i2,j1,j2,k1,k2
		 integer :: i,j,k
		 
         do k=k1,k2
              do j=j1,j2
                   do i=i1,i2
                        Fld_filtered(i,j,k)=wtm1*Fld(i,j,k-1)+wt0*Fld(i,j,k)+wtp1*Fld(i,j,k+1)
                   end do
              end do
         end do

		 Fld(i1:i2,j1:j2,k1:k2) = Fld_filtered(i1:i2,j1:j2,k1:k2)
	 end subroutine MovingAverageZ1
	 
!      subroutine MovingAverageFilter(J0)
!           real(psn),dimension(mx,my,mz) :: J0
!           real(psn),dimension(mx,my,mz) ::FldTemp
!           integer :: i,j,k
!
! 		FldTemp=0.0_psn
! 		call ExchangeYZEdgeCurrent1(J0)
! #ifndef twoD
!         do k=2,mz-2
! #else
!         do k=1,1
! #endif
!               do j=2,my-2
!                   do i=2,mx-2
!                          FldTemp(i,j,k)=wtm1*J0(i-1,j,k)+wt0*J0(i,j,k)+wtp1*J0(i+1,j,k)
!                     end do
!                end do
!           end do
!           J0=FldTemp
!
!           call ExchangeZXEdgeCurrent1(J0)
! #ifndef twoD
!           do k=2,mz-2
! #else
!           do k=1,1
! #endif
!             do j=2,my-2
!                    do i=2,mx-2
!                          FldTemp(i,j,k)=wtm1*J0(i,j-1,k)+wt0*J0(i,j,k)+wtp1*J0(i,j+1,k)
!                     end do
!                end do
!           end do
!
!           J0=FldTemp
! #ifndef twoD
!           call ExchangeXYEdgeCurrent1(J0)
!           do k=2,mz-2
!                do j=2,my-2
!                     do i=2,mx-2
!                          FldTemp(i,j,k)=wtm1*J0(i,j,k-1)+wt0*J0(i,j,k)+wtp1*J0(i,j,k+1)
!                     end do
!                end do
!           end do
!           J0=FldTemp
! #endif
!      end subroutine MovingAverageFilter
     
	 
	 subroutine CurrentFilter_BC
		 if(bc_face(2)%pos_prtl.gt.0) call FoldInCurrentRight(bc_face(2)%pos_prtl)
		 if(bc_face(1)%pos_prtl.gt.0) call FoldInCurrentLeft(bc_face(1)%pos_prtl)
	 end subroutine CurrentFilter_BC
	 
	 subroutine FoldInCurrentRight(pos)
		 real(dbpsn) :: pos
		 integer :: ind_local
		 
		 ! Assuming current ->0 near the injector, let the current filter into a buffer domain between the inecjor and Fld BC 
		 if(bc_face(2)%type_prtl.eq.'iflw') return 
		 
		 ind_local = floor(pos-xborders(procxind)+3)
		 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
			 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local+1,:,:)
			 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local+1,:,:)
			 Jy(ind_local+1,:,:) = 0 
			 Jz(ind_local+1,:,:) = 0
		 end if 

		 
		 ind_local = floor(bc_face(2)%pos_prtl-0.5_psn -xborders(procxind)+3)
		 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
			 Jx(ind_local,:,:) = Jx(ind_local,:,:) + Jx(ind_local+1,:,:)
			 Jx(ind_local+1,:,:) = 0
		 end if 
	 end subroutine FoldInCurrentRight
	 
	 subroutine FoldInCurrentLeft(pos)
		 real(dbpsn) :: pos
		 integer :: ind_local
		 
		 ind_local = ceiling( pos -xborders(procxind)+3 )
		 if(ind_local.ge.2.and.ind_local.le.mx) then 
			 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local-1,:,:)
			 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local-1,:,:)
			 Jy(ind_local-1,:,:) = 0 
			 Jz(ind_local-1,:,:) = 0
		 end if 

		 ind_local = ceiling(pos-0.5_psn -xborders(procxind)+3)
		 if(ind_local.ge.2.and.ind_local.le.mx) then 
			 Jx(ind_local,:,:) = Jx(ind_local,:,:) + Jx(ind_local-1,:,:)
			 Jx(ind_local-1,:,:) = 0
		 end if 
	 end subroutine FoldInCurrentLeft
	 
	 !--------------------------------------------------------------------------------------------------     
	 ! External Current
	 !--------------------------------------------------------------------------------------------------	
	 
     subroutine AddExternalCurrent
 		integer :: i,j,k
 		real(dbpsn) :: x,y,z
 		real(psn) :: j_x,j_y,j_z

        if(ext_current_present.eqv..false.) return
		
#ifdef twoD
         do k=1,1
#else 	
 	     do k=1,mz
#endif 
         do j=1,my 
 		    do i=1,mx
 				x= i -3.0_psn + xborders(procxind)
 				y= j -3.0_psn + yborders(procyind)
#ifdef twoD    
                z=0.0_psn
#else				 				
 				z= k - 3.0_psn + zborders(proczind)
#endif 				
 			    call J_ext(x + 0.5_dbpsn,y,z,j_x,j_y,j_z)
 				Jx(i,j,k)= Jx(i,j,k) + j_x
 			    call J_ext(x,y + 0.5_dbpsn,z,j_x,j_y,j_z)
 				Jy(i,j,k)= Jy(i,j,k) + j_y
#ifdef twoD				
 				call J_ext(x,y,0.0_dbpsn,j_x,j_y,j_z)
#else				
 				call J_ext(x,y,z + 0.5_dbpsn,j_x,j_y,j_z)
#endif				
			
 				Jz(i,j,k)= Jz(i,j,k) + j_z
 		    end do 
         end do  
         end do 
 	end subroutine AddExternalCurrent
	 

end module current