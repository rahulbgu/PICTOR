module fields
     use parameters
     use vars
     use comm_fldprtl
#ifdef cyl	 
	 use cyl_filter
#endif	 
     implicit none
	 
	 logical :: ext_current_present = .false. ! no external current by default 

	 procedure(vector_global), pointer :: J_ext => null()
contains 
     
     subroutine ResetCurrent
          Jx=0.0_psn
          Jy=0.0_psn
          Jz=0.0_psn
     end subroutine ResetCurrent
	 
     subroutine UpdateCurrentsAllEdges
          call ExchangeYZEdgeCurrent
          call AddImportedCurrentYZ
          call ExchangeZXEdgeCurrent
          call AddImportedCurrentZX          
#ifndef twoD
          call ExchangeXYEdgeCurrent
          call AddImportedCurrentXY
#endif           
     end subroutine UpdateCurrentsAllEdges
     subroutine AddImportedCurrentYZ
          integer:: i,j,k,i1

       do k=1,mz
            do j=1,my
                do i=3,5  
                    Jx(i,j,k)=Jx(i,j,k)+buff_lJx(i-2,j,k)
                    Jy(i,j,k)=Jy(i,j,k)+buff_lJy(i-2,j,k)
                    Jz(i,j,k)=Jz(i,j,k)+buff_lJz(i-2,j,k)
                 end do
            end do
        end do 
      !right edges 
     do k=1,mz
           do j=1,my
            i1=0
           do i=mx-4,mx-2  
            i1=i1+1
               Jx(i,j,k)=Jx(i,j,k)+buff_rJx(i1,j,k)
               Jy(i,j,k)=Jy(i,j,k)+buff_rJy(i1,j,k)
               Jz(i,j,k)=Jz(i,j,k)+buff_rJz(i1,j,k)
            end do
       end do
    end do 
     end subroutine AddImportedCurrentYZ
     subroutine AddImportedCurrentZX
          integer:: i,j,k,j1
       do k=1,mz
            do j=3,5
                do i=1,mx  
                    Jx(i,j,k)=Jx(i,j,k)+buff_bJx(i,j-2,k)
                    Jy(i,j,k)=Jy(i,j,k)+buff_bJy(i,j-2,k)
                    Jz(i,j,k)=Jz(i,j,k)+buff_bJz(i,j-2,k)
                 end do
            end do
        end do 
      !top edges 
     do k=1,mz
         j1=0           
           do j=my-4,my-2
              j1=j1+1
           do i=1,mx  
               Jx(i,j,k)=Jx(i,j,k)+buff_tJx(i,j1,k)
               Jy(i,j,k)=Jy(i,j,k)+buff_tJy(i,j1,k)
               Jz(i,j,k)=Jz(i,j,k)+buff_tJz(i,j1,k)
            end do
       end do
    end do 
     end subroutine AddImportedCurrentZX
     subroutine AddImportedCurrentXY
        Jx(1:mx,1:my,3:5)=Jx(1:mx,1:my,3:5)+buff_dJx(1:mx,1:my,1:3)
        Jy(1:mx,1:my,3:5)=Jy(1:mx,1:my,3:5)+buff_dJy(1:mx,1:my,1:3)
        Jz(1:mx,1:my,3:5)=Jz(1:mx,1:my,3:5)+buff_dJz(1:mx,1:my,1:3)
        
        Jx(1:mx,1:my,mz-4:mz-2)=Jx(1:mx,1:my,mz-4:mz-2)+buff_uJx(1:mx,1:my,1:3)
        Jy(1:mx,1:my,mz-4:mz-2)=Jy(1:mx,1:my,mz-4:mz-2)+buff_uJy(1:mx,1:my,1:3)
        Jz(1:mx,1:my,mz-4:mz-2)=Jz(1:mx,1:my,mz-4:mz-2)+buff_uJz(1:mx,1:my,1:3)
     end subroutine AddImportedCurrentXY
     
     subroutine smoothen_current
          integer :: ni
#ifdef cyl		  
		  call smoothen_current_cyl
#else 
		  do ni=1,curr_filter
			   call MovingAverageFilter(Jx)
			   call MovingAverageFilter(Jy)
			   call MovingAverageFilter(Jz)
			   
			   call CurrentFilter_BC
		  end do
#endif 
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
                    call MovingAverageFilter(FilteredEx)
                    call MovingAverageFilter(FilteredEy)
                    call MovingAverageFilter(FilteredEz)
               end do 
               call ExchangeYZEdgeField(FilteredEx,FilteredEy,FilteredEz)
               call ExchangeZXEdgeField(FilteredEx,FilteredEy,FilteredEz)
#ifndef twoD      
               call ExchangeXYEdgeField(FilteredEx,FilteredEy,FilteredEz)
#endif     
         end if
     end subroutine SetFilteredEfield
     
!--------------------------------------------------------------------------------------------------     
! Filters: used to smoothen grid quantitites
!--------------------------------------------------------------------------------------------------
     subroutine MovingAverageFilter(J0)
          real(psn),dimension(mx,my,mz) :: J0
          real(psn),dimension(mx,my,mz) ::FldTemp
          integer :: i,j,k
		
		FldTemp=0.0_psn
		call ExchangeYZEdgeCurrent1(J0) 		
#ifndef twoD
         do k=3,mz-3
#else
        do k=1,1
#endif
              do j=3,my-3
                  do i=3,mx-3
                         FldTemp(i,j,k)=wtm1*J0(i-1,j,k)+wt0*J0(i,j,k)+wtp1*J0(i+1,j,k)
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
     end subroutine MovingAverageFilter
     

     subroutine SyncCurrentEdges
          !used in savedata: this subroutine makes sures that current value at the outer edges, which are not 
          !updated on this proc, are updated from other proc before the data is ready to be saved 
          call ExchangeYZEdgeCurrent1(Jx)
          call ExchangeYZEdgeCurrent1(Jy)
          call ExchangeYZEdgeCurrent1(Jz)
           
          call ExchangeZXEdgeCurrent1(Jx)                    
          call ExchangeZXEdgeCurrent1(Jy)                    
          call ExchangeZXEdgeCurrent1(Jz)     
#ifndef twoD           
          call ExchangeXYEdgeCurrent1(Jx)                    
          call ExchangeXYEdgeCurrent1(Jy)                    
          call ExchangeXYEdgeCurrent1(Jz)    
#endif           
     end subroutine SyncCurrentEdges
	 
	 subroutine CurrentFilter_BC
		 if(BC_Xmax_Prtl.gt.0) call FoldInCurrentRight
		 if(BC_Xmin_Prtl.gt.0) call FoldInCurrentLeft
	 end subroutine CurrentFilter_BC
	 
	 subroutine FoldInCurrentRight
		 integer :: ind_local
		 
		 
		 if(BC_Xmax_Prtl_Type.eq.'iflw') return ! Assuming current ->0 near the injector, let the current filter into a buffer domain between the inecjor and Fld BC 
		 
		 ind_local = floor(BC_Xmax_Prtl-xborders(procxind(proc))+3)
		 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
			 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local+1,:,:)
			 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local+1,:,:)
			 Jy(ind_local+1,:,:) = 0 
			 Jz(ind_local+1,:,:) = 0
		 end if 

		 
		 ind_local = floor(BC_Xmax_Prtl-0.5_psn -xborders(procxind(proc))+3)
		 if(ind_local.ge.1.and.ind_local.le.mx-1) then 
			 Jx(ind_local,:,:) = Jx(ind_local,:,:) + Jx(ind_local+1,:,:)
			 Jx(ind_local+1,:,:) = 0
		 end if 
	 end subroutine FoldInCurrentRight
	 
	 subroutine FoldInCurrentLeft
		 integer :: ind_local
		 
		 ind_local = ceiling( BC_Xmin_Prtl  -xborders(procxind(proc))+3 )
		 if(ind_local.ge.2.and.ind_local.le.mx) then 
			 Jy(ind_local,:,:) = Jy(ind_local,:,:) + Jy(ind_local-1,:,:)
			 Jz(ind_local,:,:) = Jz(ind_local,:,:) + Jz(ind_local-1,:,:)
			 Jy(ind_local-1,:,:) = 0 
			 Jz(ind_local-1,:,:) = 0
		 end if 

		 ind_local = ceiling(BC_Xmin_Prtl-0.5_psn -xborders(procxind(proc))+3)
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
 				x= i -3.0_psn + xborders(procxind(proc))
 				y= j -3.0_psn + yborders(procyind(proc))
#ifdef twoD    
                z=0.0_psn
#else				 				
 				z= k - 3.0_psn + zborders(proczind(proc))
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
	 

end module fields