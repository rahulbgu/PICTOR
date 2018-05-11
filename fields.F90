module fields
     use parameters
     use vars
     use comm_fldprtl
     implicit none
contains 
     subroutine UpdateEfield
		 integer :: i,j,k
          integer :: i1,j1,k1,i2,j2,k2
          
          i1=3
          i2=mx-3
          j1=3
          j2=my-3

#ifndef twoD
          k1=3
          k2=mz-3
#else          
          k1=1
          k2=1
#endif

        do k=k1,k2
             do j=j1,j2
                  do i=i1,i2
#ifndef twoD
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))-fldc*(By(i,j,k)-By(i,j,k-1))
                    Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
#else
                    Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))
                    Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
                    Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
#endif
                    end do 
               end do
          end do
     end subroutine UpdateEfield 
     subroutine UpdateBfieldHalf
		  integer ::i,j,k
          integer :: i1,j1,k1,i2,j2,k2          
          i1=1
          i2=mx-1
          j1=1
          j2=my-1
          
#ifndef twoD
          k1=1
          k2=mz-1
#else
        k1=1
        k2=1
#endif    

          do k=k1,k2
               do j=j1,j2
                    do i=i1,i2
#ifndef twoD
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))+fld_halfc*(Ey(i,j,k+1)-Ey(i,j,k))
                         By(i,j,k)=By(i,j,k)-fld_halfc*(Ex(i,j,k+1)-Ex(i,j,k))+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
#else
                         Bx(i,j,k)=Bx(i,j,k)-fld_halfc*(Ez(i,j+1,k)-Ez(i,j,k))
                         By(i,j,k)=By(i,j,k)+fld_halfc*(Ez(i+1,j,k)-Ez(i,j,k))
                         Bz(i,j,k)=Bz(i,j,k)-fld_halfc*(Ey(i+1,j,k)-Ey(i,j,k))+fld_halfc*(Ex(i,j+1,k)-Ex(i,j,k))
#endif
                    end do
               end do
          end do 
     end subroutine UpdateBfieldHalf
     
          
     
!---------------------------  END OF FLD UPDATE SUBROUTINES --------------------------------------------

!=======================================================================================================
!
! Radiation Boundary condition: Surface subroutine of TRISTAN can be used to implement radiaiton boundary condition 
!=======================================================================================================

! subroutine Surface(bx,by,bz,ex,ey,ez,ix,iy,iz,mx,my,mz,ind)
! 	real(psn), dimension(1):: bx,by,bz,ex,ey,ez
! 	integer :: ind,ix,iy,iz,mx,my,mz
! 	real(psn) :: rs,s,os
!     rs=2.0_psn*c/(1.0_psn+c)
!     s=0.4142136_psn
!     os=0.5_psn*(1.0_psn-s)*rs
!
!
!
!      do m=ind+iz*(mz-1),ind+iz*(mz-1)+iy*(my-2),iy
!         do n=m+ix,m+ix*(mx-2),ix
!            bxp(n)=bxp(n)+rs*(bx(n-iz)-bxp(n)+s*(bzp(n)-bz(n-ix)))-os*(ez(n+iy)-ezp(n))-(os-c)*(ez(n+iy-iz)-ez(n-iz))-c*(eyp(n)-ey(n-iz))
! 	    end do
! 	end do
!
! 	do m=ind+iz*(mz-1),ind+iz*(mz-1)+ix*(mx-2),ix
!         do n=m+iy,m+iy*(my-2),iy
!            byp(n)=byp(n)+rs*(by(n-iz)-byp(n)+s*(bzp(n)-bz(n-iy)))+os*(ez(n+ix)-ezp(n))+(os-c)*(ez(n+ix-iz)-ez(n-iz))+c*(exp(n)-ex(n-iz))
! 	    end do
! 	end do
!
! end subroutine Surface
     

     subroutine AddCurrent
          integer:: i,j,k                    

#ifndef twoD 
          do k=3,mz-3
#else
          do k=1,1
#endif
            do j=3,my-3
                  do i=3,mx-3

                         Ex(i,j,k)=Ex(i,j,k)-Jx(i,j,k)
                         Ey(i,j,k)=Ey(i,j,k)-Jy(i,j,k)
                         Ez(i,j,k)=Ez(i,j,k)-Jz(i,j,k)
                    end do
               end do
          end do
     end subroutine AddCurrent
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
          !call FoldInCurrent
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
        Jx(:,:,3:5)=Jx(:,:,3:5)+buff_dJx(:,:,1:3)
        Jy(:,:,3:5)=Jy(:,:,3:5)+buff_dJy(:,:,1:3)
        Jz(:,:,3:5)=Jz(:,:,3:5)+buff_dJz(:,:,1:3)
        
        Jx(:,:,mz-4:mz-2)=Jx(:,:,mz-4:mz-2)+buff_uJx(:,:,1:3)
        Jy(:,:,mz-4:mz-2)=Jy(:,:,mz-4:mz-2)+buff_uJy(:,:,1:3)
        Jz(:,:,mz-4:mz-2)=Jz(:,:,mz-4:mz-2)+buff_uJz(:,:,1:3)
     end subroutine AddImportedCurrentXY
     
     subroutine smoothen_current
          integer :: ni
          do ni=1,curr_filter
               call MovingAverageFilter(Jx)
               call MovingAverageFilter(Jy)
               call MovingAverageFilter(Jz)
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
	 
	 
!---------------------------------------------------------------------------------------------------------	 
! !  The following filtering scheme is used in Shock Simulations
! !---------------------------------------------------------------------------------------------------------
! subroutine smoothen_current_subdomain(nf,x1,x2)
! 	 integer, intent(IN) :: nf,x1,x2 !x1 and x2 are in local cordinates
! 	 !integer :: imax,imin
! 	 integer :: n
!      do n=1,nf
! ! 		  if((x1.ge.3).and.(x1.le.mx-2)) then
! ! 			  Jx(x1-1,:,:)=-Jx(x1+1,:,:)
! ! 			  Jy(x1-1,:,:)=Jy(x1,:,:)
! ! 			  Jz(x1-1,:,:)=Jz(x1,:,:)
! ! 		  end if
! !  		  if((x2.ge.2).and.(x2.le.mx-2)) then
! ! 			  Jx(x2+1,:,:)=-Jx(x2,:,:)
! ! 			  Jy(x2+1,:,:)=Jy(x2,:,:)
! ! 			  Jz(x2+1,:,:)=Jz(x2,:,:)
! ! 		  end if
!           call MovingAverageFilterSubDomain(Jx,x1,x2)
!           call MovingAverageFilterSubDomain(Jy,x1,x2)
!           call MovingAverageFilterSubDomain(Jz,x1,x2)
! ! 		  if((x1.ge.3).and.(x1.le.mx-2)) then
! ! 			  Jx(x1-1,:,:)=0.0_psn
! ! 			  Jy(x1-1,:,:)=0.0_psn
! ! 			  Jz(x1-1,:,:)=0.0_psn
! ! 		  end if
! !  		  if((x2.ge.2).and.(x2.le.mx-2)) then
! ! 			  Jx(x2+1,:,:)=0.0_psn
! ! 			  Jy(x2+1,:,:)=0.0_psn
! ! 			  Jz(x2+1,:,:)=0.0_psn
! ! 		  end if
!      end do
! end subroutine smoothen_current_subdomain
!
!
! !The following scheme does not work! There some chagne accumulation at the current edges
! 	 subroutine MovingAverageFilterSubDomain(J0,x1,x2)
!           real(psn),dimension(mx,my,mz), intent(INOUT):: J0
! 		  integer, intent(IN) :: x1,x2
! 		  integer :: imin,imax
!           real(psn),dimension(mx,my,mz) ::FldTemp
!           integer :: i,j,k
!
! 	      imin=max(3,x1-1)
! 	      imax=min(mx-3,x2-1)
!
! 		FldTemp=0.0_psn
! 		call ExchangeYZEdgeCurrent1(J0)
! #ifndef twoD
!         do k=3,mz-3
! #else
!         do k=1,1
! #endif
!                do j=3,my-3
!                   do i=imin,imax
! 						 FldTemp(i,j,k)=wtm*J0(i-1,j,k)+wt0*J0(i,j,k)+wtp*J0(i+1,j,k)
!                     end do
!                end do
!           end do
! 		  J0=FldTemp
!           call ExchangeZXEdgeCurrent1(J0)
! 		  if((x1.ge.2).and.(x1.le.mx-2)) then
! 			    J0(x1,:,:)=J0(x1,:,:)+J0(x1-1,:,:)
! 		        J0(x1-1,:,:)=0.0_psn
! 		  end if
! 		  if((x2.ge.2).and.(x2.le.mx-2)) then
! 		       J0(x2,:,:)=J0(x2,:,:)+J0(x2+1,:,:)
! 	           J0(x2+1,:,:)=0.0_psn
! 		  end if
!           call ExchangeZXEdgeCurrent1(J0)
!
! #ifndef twoD
!           do k=3,mz-3
! #else
!           do k=1,1
! #endif
!             do j=3,my-3
!                    do i=imin,imax
!                          FldTemp(i,j,k)=wtm*J0(i,j-1,k)+wt0*J0(i,j,k)+wtp*J0(i,j+1,k)
!                     end do
!                end do
!           end do
!           J0=FldTemp
! #ifndef twoD
!
!           call ExchangeXYEdgeCurrent1(J0)
!           do k=3,mz-3
!                do j=3,my-3
!                     do i=imin,imax
!                          FldTemp(i,j,k)=wtm*J0(i,j,k-1)+wt0*J0(i,j,k)+wtp*J0(i,j,k+1)
!                     end do
!                end do
!           end do
!           J0=FldTemp
! #endif
!      end subroutine MovingAverageFilterSubDomain
 
	 
     
!      subroutine CopyCurrentLayerXY(zs,zd,Jthis) !Not needed
!           real(psn),dimension(mx,my,mz) :: Jthis
!           integer:: zs,zd,i,j
!           do j=1,my
!                do i=1,mx
!                     Jthis(i,j,zd)=Jthis(i,j,zs)
!                end do
!           end do
!      end subroutine CopyCurrentLayerXY
     
     !      subroutine FoldInCurrent
     ! !           call AddCurrentLayerZX(my,5)
     ! !           call AddCurrentLayerZX(my-1,4)
     ! !           call AddCurrentLayerZX(my-2,3)
     ! !           call AddCurrentLayerZX(1,my-4)
     ! !           call AddCurrentLayerZX(2,my-3)
     ! #ifndef twoD
     !           call AddCurrentLayerXY(mz,5)
     !           call AddCurrentLayerXY(mz-1,4)
     !           call AddCurrentLayerXY(mz-2,3)
     !           call AddCurrentLayerXY(1,mz-4)
     !           call AddCurrentLayerXY(2,mz-3)
     ! #endif
     !
     !           !comment out if more than one proc
     ! !           call AddCurrentLayerYZ(mx,5)
     ! !           call AddCurrentLayerYZ(mx-1,4)
     ! !           call AddCurrentLayerYZ(mx-2,3)
     ! !           call AddCurrentLayerYZ(1,mx-4)
     ! !           call AddCurrentLayerYZ(2,mx-3)
     !      end subroutine FoldInCurrent

     !      subroutine UpdateEBfldsXYEdge ! Not needed
     ! #ifndef twoD
     !           call CopyEBFldsEdgeXY(3,mz-2)
     !           call CopyEBFldsEdgeXY(4,mz-1)
     !           call CopyEBFldsEdgeXY(5,mz)
     !           call CopyEBFldsEdgeXY(mz-3,2)
     !           call CopyEBFldsEdgeXY(mz-4,1)
     ! #endif
     !      end subroutine UpdateEBFldsXYEdge
     !      subroutine CopyEBFldsEdgeXY(zs,zd) ! Not needed
     !           integer:: zs,zd,i,j
     !           do j=1,my
     !               do i=1,mx
     !                     Ex(i,j,zd)=Ex(i,j,zs)
     !                     Ey(i,j,zd)=Ey(i,j,zs)
     !                     Ez(i,j,zd)=Ez(i,j,zs)
     !                     Bx(i,j,zd)=Bx(i,j,zs)
     !                     By(i,j,zd)=By(i,j,zs)
     !                     Bz(i,j,zd)=Bz(i,j,zs)
     !                end do
     !           end do
     !      end subroutine CopyEBFldsEdgeXY

     !      subroutine CopyLayerYZ(xs,xd)
     !           integer:: xs,xd,j,k
     !           do k=1,mz
     !               do j=1,my
     !                     Ex(xd,j,k)=Ex(xs,j,k)
     !                     Ey(xd,j,k)=Ey(xs,j,k)
     !                     Ez(xd,j,k)=Ez(xs,j,k)
     !                     Bx(xd,j,k)=Bx(xs,j,k)
     !                     By(xd,j,k)=By(xs,j,k)
     !                     Bz(xd,j,k)=Bz(xs,j,k)
     !                end do
     !           end do
     !      end subroutine CopyLayerYZ
     !      subroutine CopyFldEdgeZX(ys,yd)
     !           integer:: ys,yd,i,k
     !          do k=1,mz
     !              do i=1,mx
     !
     !                     Ex(i,yd,k)=Ex(i,ys,k)
     !                     Ey(i,yd,k)=Ey(i,ys,k)
     !                     Ez(i,yd,k)=Ez(i,ys,k)
     !                     Bx(i,yd,k)=Bx(i,ys,k)
     !                     By(i,yd,k)=By(i,ys,k)
     !                     Bz(i,yd,k)=Bz(i,ys,k)
     !                end do
     !           end do
     !      end subroutine CopyFldEdgeZX
     !      subroutine AddCurrentLayerZX(ys,yd)
     !           integer:: ys,yd,i,k
     !           do k=1,mz
     !               do i=1,mx
     !                     Jx(i,yd,k)=Jx(i,yd,k)+Jx(i,ys,k)
     !                     Jy(i,yd,k)=Jy(i,yd,k)+Jy(i,ys,k)
     !                     Jz(i,yd,k)=Jz(i,yd,k)+Jz(i,ys,k)
     !                end do
     !           end do
     !      end subroutine AddCurrentLayerZX
     !      subroutine AddCurrentLayerYZ(xs,xd)
     !           integer:: xs,xd,j,k
     !           do k=1,mz
     !                do j=1,my
     !                     Jx(xd,j,k)=Jx(xd,j,k)+Jx(xs,j,k)
     !                     Jy(xd,j,k)=Jy(xd,j,k)+Jy(xs,j,k)
     !                     Jz(xd,j,k)=Jz(xd,j,k)+Jz(xs,j,k)
     !                end do
     !           end do
     !      end subroutine AddCurrentLayerYZ
     !      subroutine AddCurrentLayerXY(zs,zd) !Not Needed
     !           integer:: zs,zd,i,j
     !           do j=1,my
     !                do i=1,mx
     !                     Jx(i,j,zd)=Jx(i,j,zd)+Jx(i,j,zs)
     !                     Jy(i,j,zd)=Jy(i,j,zd)+Jy(i,j,zs)
     !                     Jz(i,j,zd)=Jz(i,j,zd)+Jz(i,j,zs)
     !                end do
     !           end do
     !      end subroutine AddCurrentLayerXY
     
     
!------------------------------------------------------------------------------------------------------
! The following subroutines are useful to apply filter to Electromagentic field compoents 
! a) only the transverse compoents of electric field is smoothened to satisfy div E = \rho and div B=0
!------------------------------------------------------------------------------------------------------     
     
!      subroutine FilterEMfld
!           integer :: ni
!           do ni=1,nEMfilter
!                call FilterEMfldX(Ey)
!                call FilterEMfldX(Ez)
!                call FilterEMfldY(Ex)
!                call FilterEMfldY(Ez)
! #ifndef twoD
!                call FilterEMfldZ(Ex)
!                call FilterEMfldZ(Ey)
! #endif
!                call FilterEMfldX(By)
!                call FilterEMfldX(Bz)
!                call FilterEMfldY(Bx)
!                call FilterEMfldY(Bz)
! #ifndef twoD
!                call FilterEMfldZ(Bx)
!                call FilterEMfldZ(By)
! #endif
!            !now update all the outer edges
!                call ExchangeYZEdgeEMFld1(Ey)
!                call ExchangeYZEdgeEMFld1(Ez)
!                call ExchangeZXEdgeEMfld1(Ex)
!                call ExchangeZXEdgeEMfld1(Ez)
! #ifndef twoD
!             call CopyFldLayerXY(3,mz-2,Ex)
!                call CopyFldLayerXY(mz-3,2,Ex)
!             call CopyFldLayerXY(3,mz-2,Ey)
!                call CopyFldLayerXY(mz-3,2,Ey)
! #endif
!                call ExchangeYZEdgeEMFld1(By)
!                call ExchangeYZEdgeEMFld1(Bz)
!                call ExchangeZXEdgeEMfld1(Bx)
!                call ExchangeZXEdgeEMfld1(Bz)
! #ifndef twoD
!             call CopyFldLayerXY(3,mz-2,Bx)
!                call CopyFldLayerXY(mz-3,2,Bx)
!             call CopyFldLayerXY(3,mz-2,By)
!                call CopyFldLayerXY(mz-3,2,By)
! #endif
!           end do
!      end subroutine FilterEMfld
!      subroutine FilterEMfldX(Fld)
!           real(psn),dimension(mx,my,mz) :: Fld
!           real(psn),dimension(mx,my,mz) :: Fldtemp
!           integer :: i,j
! #ifndef twoD
!           do k=3,mz-3
! #else
!         do k=1,1
! #endif
!                do j=3,my-3
!                     do i=3,mx-3
!                          Fldtemp(i,j,k)=wtm*Fld(i-1,j,k)+wt0*Fld(i,j,k)+wtp*Fld(i+1,j,k)
!                     end do
!                end do
!           end do
!      end subroutine FilterEMfldX
!      subroutine FilterEMfldY(Fld)
!           real(psn),dimension(mx,my,mz) :: Fld
!           real(psn),dimension(mx,my,mz) :: Fldtemp
!           integer :: i,j
! #ifndef twoD
!           do k=3,mz-3
! #else
!         do k=1,1
! #endif
!                do j=3,my-3
!                     do i=3,mx-3
!                          Fldtemp(i,j,k)=wtm*Fld(i,j-1,k)+wt0*Fld(i,j,k)+wtp*Fld(i,j+1,k)
!                     end do
!                end do
!           end do
!      end subroutine FilterEMfldY
!      subroutine FilterEMfldZ
!           real(psn),dimension(mx,my,mz) :: Fld
!           real(psn),dimension(mx,my,mz) :: Fldtemp
!           integer :: i,j
!
!           do k=3,mz-3
!                do j=3,my-3
!                     do i=3,mx-3
!                          Fldtemp(i,j,k)=wtm*Fld(i,j,k-1)+wt0*Fld(i,j,k)+wtp*Fld(i,j,k+1)
!                     end do
!                end do
!           end do
!      end subroutine FilterEMfldZ
!
!      subroutine CopyFldLayerXY(zs,zd,Fld)
!           real(psn),dimension(mx,my,mz) :: Fld
!           integer:: zs,zd,i,j
!           do j=1,my
!                do i=1,mx
!                     Fld(i,j,zd)=Fld(i,j,zs)
!                end do
!           end do
!      end subroutine CopyFldLayerXY


!========================================================================================================     
! The following subroutines are useful for the Hybrid1 case 
!========================================================================================================     
! #ifdef Hybrid1
!
! subroutine UpdateEfieldHybrid1
!           integer :: i1,j1,k1,i2,j2,k2
!          real(psn) :: thisBx,thisBy,thisBz
!           real(psn) :: curlBx,curlBy,curlBz
!
!           i1=3
!           i2=mx-3
!           j1=3
!           j2=my-3
!
! #ifndef twoD
!           k1=3
!           k2=mz-3
! #else
!           k1=1
!           k2=1
! #endif
!         do k=k1,k2
!              do j=j1,j2
!                   do i=i1,i2
! #ifndef twoD
!                     thisBx=(  (Bx(i,j,k)+Bx(i,j,k-1)+Bx(i,j-1,k)+Bx(i,j-1,k-1)) +
!                                    (Bx(i+1,j,k)+Bx(i+1,j,k-1)+Bx(i+1,j-1,k)+Bx(i+1,j-1,k-1)  )*0.125
!                          thisBy=(By(i,j,k)+By(i,j,k-1))*0.5
!                          thisBz=(Bz(i,j,k)+Bz(i,j-1,k))*0.5
!                          curlBx=(Bz(i,j,k)-Bz(i,j-1,k)) - (By(i,j,k)-By(i,j,k-1))
!                          curlBy=
!
!                          Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))-fldc*(By(i,j,k)-By(i,j,k-1))
!                          Ey(i,j,k)=Ey(i,j,k)+fldc*(Bx(i,j,k)-Bx(i,j,k-1))-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
!                          Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
! #else
!                     Ex(i,j,k)=Ex(i,j,k)+fldc*(Bz(i,j,k)-Bz(i,j-1,k))
!                     Ey(i,j,k)=Ey(i,j,k)-fldc*(Bz(i,j,k)-Bz(i-1,j,k))
!                     Ez(i,j,k)=Ez(i,j,k)+fldc*(By(i,j,k)-By(i-1,j,k))-fldc*(Bx(i,j,k)-Bx(i,j-1,k))
! #endif
!                     end do
!                end do
!           end do
! end subroutine UpdateEfieldHybrid1
!
! #endif
!
     
     
!=========================================================================================================
! The following subroutines are useful for the Hybrid2 case 
! some of them can be combined with normal subroutines written for PIC case, check for code length optimisations,but maybe 
!it helps the code run faster 
!=========================================================================================================
! #ifdef Hybrid2
!
! ! subroutine UpdateEfieldHybrid2
! !      integer :: i1,j1,k1,i2,j2,k2
! !      real(psn)::mean_DEN
! !
! !
! !           i1=3
! !           i2=mx-3
! !           j1=3
! !           j2=my-3
! !
! ! #ifndef twoD
! !           k1=3
! !           k2=mz-3
! ! #else
! !           k1=1
! !           k2=1
! ! #endif
! !         do k=k1,k2
! !              do j=j1,j2
! !                   do i=i1,i2
! !                          !if(Deny_hbd2(i,j,k).eq.0) print*,'In field1 Zero Deny',Deny_hbd2(i,j,k),'at',proc,'in',i,j,k
! !
! ! #ifndef twoD
! !                         Ex(i,j,k)= (Jy(i,j,k)+Jy(i+1,j,k))*(Bz(i,j,k)+Bz_ext0) + (Jy(i,j-1,k)+Jy(i+1,j-1,k))*(Bz(i,j-1,k)+Bz_ext0)&
! !                                    -(Jz(i,j,k)+Jz(i+1,j,k))*(By(i,j,k)+By_ext0) - (Jz(i,j,k-1)+Jz(i+1,j,k-1))*(By(i,j,k-1)+By_ext0)
! !                          Ey(i,j,k)= (Jz(i,j,k)+Jz(i,j+1,k))*(Bx(i,j,k)+Bx_ext0) + (Jz(i,j,k-1)+Jz(i,j+1,k-1))*(Bx(i,j,k-1)+Bx_ext0)&
! !                                    -(Jx(i,j,k)+Jx(i,j+1,k))*(Bz(i,j,k)+Bz_ext0) - (Jx(i-1,j,k)+Jx(i-1,j+1,k))*(Bz(i-1,j,k)+Bz_ext0)
! !                         Ez(i,j,k)= (Jx(i,j,k)+Jx(i,j,k+1))*(By(i,j,k)+By_ext0) + (Jx(i-1,j,k)+Jx(i-1,j,k+1))*(By(i-1,j,k)+By_ext0)&
! !                                    -(Jy(i,j,k)+Jy(i,j,k+1))*(Bx(i,j,k)+Bx_ext0) - (Jy(i,j-1,k)+Jy(i,j-1,k+1))*(Bx(i,j-1,k)+Bx_ext0)
! !
! !                          Ex(i,j,k)= -Ex(i,j,k)&
! !                                     +( Bx(i,j,k)    -Bx(i,j,k-1)  +Bz(i-1,j,k)   +Bx(i+1,j,k)  -Bx(i+1,j,k-1)  -Bz(i+1,j,k)    )*(Bz(i,j,k)+Bz_ext0)&
! !                                       +( Bx(i,j-1,k)  -Bx(i,j-1,k-1)+Bz(i-1,j-1,k) +Bx(i+1,j-1,k)-Bx(i+1,j-1,k-1)-Bz(i+1,j-1,k)  )*(Bz(i,j-1,k)+Bz_ext0)&
! !                                       -(-By(i-1,j,k)  -Bx(i,j,k)    +Bx(i,j-1,k)   +By(i+1,j,k)  -Bx(i+1,j,k)    +Bx(i+1,j-1,k)  )*(By(i,j,k)+By_ext0)&
! !                                       -(-By(i-1,j,k-1)-Bx(i,j,k-1)  +Bx(i,j-1,k-1) +By(i+1,j,k-1)-Bx(i+1,j,k-1)  +Bx(i+1,j-1,k-1))*(By(i,j,k-1)+By_ext0)&
! !
! !                          Ey(i,j,k)= -Ey(i,j,k)&
! !                                     +( By(i,j,k)    -By(i-1,j,k)  +Bx(i,j-1,k)   +By(i,j+1,k)  -By(i-1,j+1,k)  -Bx(i,j+1,k)    )*(Bx(i,j,k)+Bx_ext0)&
! !                                       +( By(i,j,k-1)  -By(i-1,j,k-1)+Bx(i,j-1,k-1) +By(i,j+1,k-1)-By(i-1,j+1,k-1)-Bx(i,j+1,k-1)  )*(Bx(i,j,k-1)+Bx_ext0)&
! !                                       -(-Bz(i,j-1,k)  -By(i,j,k)    +By(i,j,k-1)   +Bz(i,j+1,k)  -By(i,j+1,k)    +By(i,j+1,k-1)  )*(Bz(i,j,k)+Bz_ext0)&
! !                                       -(-Bz(i-1,j-1,k)-By(i-1,j,k)  +By(i-1,j,k-1) +Bz(i-1,j+1,k)-By(i-1,j+1,k)  +By(i-1,j+1,k-1))*(Bz(i-1,j,k)+Bz_ext0)
! !
! !                          Ez(i,j,k)= -Ez(i,j,k)&
! !                                     +( Bz(i,j,k)    -Bz(i,j-1,k)  +By(i,j,k-1)   +Bz(i,j,k+1)  -Bz(i,j-1,k+1)  -By(i,j,k+1)    )*(By(i,j,k)+By_ext0)&
! !                                       +( Bz(i-1,j,k)  -Bz(i-1,j-1,k)+By(i-1,j,k-1) +Bz(i-1,j,k+1)-Bz(i-1,j-1,k+1)-By(i-1,j,k+1)  )*(By(i-1,j,k)+By_ext0)&
! !                                       -(-Bx(i,j,k-1)  -Bz(i,j,k)    +Bz(i-1,j,k)   +Bx(i,j,k+1)  -Bz(i,j,k+1)    +Bz(i-1,j,k+1)  )*(Bx(i,j,k)+Bx_ext0)&
! !                                       -(-Bx(i,j-1,k-1)-Bz(i,j-1,k)  +Bz(i-1,j-1,k) +Bx(i,j-1,k+1)-Bz(i,j-1,k+1)  +Bz(i-1,j-1,k+1))*(Bx(i,j-1,k)+Bx_ext0)
! !
! ! !                          Ex= 0.25*Ex/Denx_hbd2
! ! !                          Ey= 0.25*Ey/Deny_hbd2
! ! !                          Ez= 0.25*Ez/Denz_hbd2
! !
! ! #else
! !
! !                          Ex(i,j,k)= (Jy(i,j,k)+Jy(i+1,j,k))*(Bz(i,j,k)+Bz_ext0) + (Jy(i,j-1,k)+Jy(i+1,j-1,k))*(Bz(i,j-1,k)+Bz_ext0)&
! !                                    -2*(Jz(i,j,k)+Jz(i+1,j,k))*(By(i,j,k)+By_ext0)
! !                          Ey(i,j,k)= 2*(Jz(i,j,k)+Jz(i,j+1,k))*(Bx(i,j,k)+Bx_ext0)&
! !                                    -(Jx(i,j,k)+Jx(i,j+1,k))*(Bz(i,j,k)+Bz_ext0) - (Jx(i-1,j,k)+Jx(i-1,j+1,k))*(Bz(i-1,j,k)+Bz_ext0)
! !                          Ez(i,j,k)= (Jx(i,j,k))*(By(i,j,k)+By_ext0) + (Jx(i-1,j,k))*(By(i-1,j,k)+By_ext0)&
! !                                    -(Jy(i,j,k))*(Bx(i,j,k)+By_ext0) - (Jy(i,j-1,k))*(Bx(i,j-1,k)+Bx_ext0)
! !
! !                          Ex(i,j,k)= -Ex(i,j,k)&
! !                                       +( Bz(i-1,j,k)  -Bz(i+1,j,k)    )*(Bz(i,j,k)+Bz_ext0)&
! !                                         +( Bz(i-1,j-1,k)-Bz(i+1,j-1,k)  )*(Bz(i,j-1,k)+Bz_ext0)&
! !                                       -2*(-By(i-1,j,k)  -Bx(i,j,k)    +Bx(i,j-1,k)   +By(i+1,j,k)  -Bx(i+1,j,k)    +Bx(i+1,j-1,k)  )*(By(i,j,k)+By_ext0)
! !
! !
! !                          Ey(i,j,k)= -Ey(i,j,k)&
! !                                     +2*( By(i,j,k)    -By(i-1,j,k)  +Bx(i,j-1,k)   +By(i,j+1,k)  -By(i-1,j+1,k)  -Bx(i,j+1,k)    )*(Bx(i,j,k)+Bx_ext0)&
! !                                         -(-Bz(i,j-1,k)  +Bz(i,j+1,k)   )*(Bz(i,j,k)+Bz_ext0)&
! !                                         -(-Bz(i-1,j-1,k)+Bz(i-1,j+1,k) )*(Bz(i-1,j,k)+Bz_ext0)
! !
! !                          Ez(i,j,k)= -Ez(i,j,k)&
! !                                     +( Bz(i,j,k)    -Bz(i,j-1,k)  )*(By(i,j,k)+By_ext0)&
! !                                       +( Bz(i-1,j,k)  -Bz(i-1,j-1,k))*(By(i-1,j,k)+By_ext0)&
! !                                       -(-Bz(i,j,k)    +Bz(i-1,j,k)  )*(Bx(i,j,k)+Bx_ext0)&
! !                                       -(-Bz(i,j-1,k)  +Bz(i-1,j-1,k))*(Bx(i,j-1,k)+Bx_ext0)
! !
! ! !                          Ex= 0.25*Ex/Denx_hbd2
! ! !                          Ey= 0.25*Ey/Deny_hbd2
! ! !                          Ez= 0.5 *Ez/Deny_hbd2
! ! !find zero density
! ! if(Denx_hbd2(i,j,k).eq.0) print*,'Zero Denx',Denx_hbd2(i,j,k),'at',proc,'in',i,j,k
! ! if(Deny_hbd2(i,j,k).eq.0) print*,'In field Zero Deny',Deny_hbd2(i,j,k),'at',proc,'in',i,j,k
! ! if(Denz_hbd2(i,j,k).eq.0) print*,'Zero Denz',Denz_hbd2(i,j,k),'at',proc,'in',i,j,k
! !
! !
! ! #endif
! ! mean_DEN=1
! ! Ex(i,j,k)= 0.25*Ex(i,j,k)/mean_DEN!/(Denx_hbd2(i,j,k)*c)
! ! Ey(i,j,k)= 0.25*Ey(i,j,k)/mean_DEN!/(Deny_hbd2(i,j,k)*c)
! ! Ez(i,j,k)= 0.5 *Ez(i,j,k)/mean_DEN!/(Deny_hbd2(i,j,k)*c)
! ! !if(isnan(Ex(i,j,k))) Ex(i,j,k)=0
! ! !if(isnan(Ey(i,j,k))) Ey(i,j,k)=0
! ! !if(isnan(Ez(i,j,k))) Ez(i,j,k)=0
! !
! !
! !
! !                     end do
! !                end do
! !           end do
! !
! !
! !           !Ex= 0.25*Ex/Denx_hbd2
! !           !Ey= 0.25*Ey/Deny_hbd2
! !           !Ez= 0.5 *Ez/Deny_hbd2
! !
! ! end subroutine UpdateEfieldHybrid2
! !
!
!
! subroutine UpdateEfieldHybrid2
!      real(psn) :: thisBx,thisBy,thisBz
!      real(psn) :: curlBx,curlBy,curlBz
!      integer   :: i1,j1,k1,i2,j2,k2
!
!
!      i1=3
!      i2=mx-3
!      j1=3
!      j2=my-3
! #ifndef twoD
!      k1=3
!      k2=mz-3
! #else
!      k1=1
!      k2=1
! #endif
!     do k=k1,k2
!         do j=j1,j2
!              do i=i1,i2
!                !magnetic field at (i,j,k)
!                thisBx=( Bx(i,j,k)+Bx(i,j,k-1)+Bx(i,j-1,k)+Bx(i,j-1,k-1) )*0.25 +Bx_ext0
!                thisBy=( By(i,j,k)+By(i,j,k-1)+By(i-1,j,k)+By(i-1,j,k-1) )*0.25 +By_ext0
!                thisBz=( Bz(i,j,k)+Bz(i,j-1,k)+Bz(i-1,j,k)+Bz(i-1,j-1,k) )*0.25 +Bz_ext0
!
!                !curl of magnetic field at (i,j,k)
!                curlBx=( (Bz(i,j,k)  -Bz(i,j-1,k))   - (By(i,j,k)  -By(i,j,k-1))   &
!                        +(Bz(i-1,j,k)-Bz(i-1,j-1,k)) - (By(i-1,j,k)-By(i-1,j,k-1)) )*0.5
!                curlBy=( (Bx(i,j,k)  -Bx(i,j,k-1))   - (Bz(i,j,k)  -Bz(i-1,j,k))   &
!                        +(Bx(i,j-1,k)-Bx(i,j-1,k-1)) - (Bz(i,j-1,k)-Bz(i-1,j-1,k)) )*0.5
!                curlBz=( (By(i,j,k)  -By(i-1,j,k))   - (Bx(i,j,k)  -Bx(i,j-1,k))   &
!                        +(By(i,j,k-1)-By(i-1,j,k-1)) - (Bx(i,j,k-1)-Bx(i,j-1,k-1)) )*0.5
!
!                Ex(i,j,k)= (-Jy(i,j,k)+curlBy)*thisBz - (-Jz(i,j,k)+curlBz)*thisBy
!                Ey(i,j,k)= (-Jz(i,j,k)+curlBz)*thisBx - (-Jx(i,j,k)+curlBx)*thisBz
!                Ez(i,j,k)= (-Jx(i,j,k)+curlBx)*thisBy - (-Jy(i,j,k)+curlBy)*thisBx
!               end do
!           end do
!      end do
!
!
!
! end subroutine UpdateEfieldHybrid2
!
!
!
!
!
!
!
!
!
!
! subroutine ResetDensityHybrid2
!      Denx_hbd2=0
!      Deny_hbd2=0
!      Denz_hbd2=0
! end subroutine ResetDensityHybrid2
!
! subroutine smoothen_densityHybrid2
!      integer :: ni
!      do ni=1,curr_filter
!           !uncomment when more than one proc
!           call avg_curr(Denx_hbd2)
!           call avg_curr(Deny_hbd2)
!           call avg_curr(Denz_hbd2)
!      end do
! end subroutine smoothen_densityHybrid2
!
! subroutine UpdateDensityAllEdgesHybrid2
!      call ExchangeYZEdgeDensity_Hybrid2
!       call AddImportedDensityYZ_Hybrid2
!      call ExchangeZXEdgeDensity_Hybrid2
!       call AddImportedDensityZX_Hybrid2
!      call FoldInDensity_Hybrid2
! end subroutine UpdateDensityAllEdgesHybrid2
!
! subroutine AddImportedDensityYZ_Hybrid2
!      integer:: i,j,k,i1
!    do k=1,mz
!        do j=1,my
!            do i=3,5
!                Denx_hbd2(i,j,k)=Denx_hbd2(i,j,k)+buff_lDenx_hbd2(i-2,j,k)
!                Deny_hbd2(i,j,k)=Deny_hbd2(i,j,k)+buff_lDeny_hbd2(i-2,j,k)
!                Denz_hbd2(i,j,k)=Denz_hbd2(i,j,k)+buff_lDenz_hbd2(i-2,j,k)
!             end do
!        end do
!    end do
!  !right edges
!  do k=1,mz
!     do j=1,my
!        i1=0
!       do i=mx-4,mx-2
!        i1=i1+1
!           Denx_hbd2(i,j,k)=Denx_hbd2(i,j,k)+buff_rDenx_hbd2(i1,j,k)
!           Deny_hbd2(i,j,k)=Deny_hbd2(i,j,k)+buff_rDeny_hbd2(i1,j,k)
!           Denz_hbd2(i,j,k)=Denz_hbd2(i,j,k)+buff_rDenz_hbd2(i1,j,k)
!        end do
!    end do
! end do
! end subroutine AddImportedDensityYZ_Hybrid2
! subroutine AddImportedDensityZX_Hybrid2
!      integer:: i,j,k,j1
!    do k=1,mz
!        do j=3,5
!            do i=1,mx
!                Denx_hbd2(i,j,k)=Denx_hbd2(i,j,k)+buff_bDenx_hbd2(i,j-2,k)
!                Deny_hbd2(i,j,k)=Deny_hbd2(i,j,k)+buff_bDeny_hbd2(i,j-2,k)
!                Denz_hbd2(i,j,k)=Denz_hbd2(i,j,k)+buff_bDenz_hbd2(i,j-2,k)
!             end do
!        end do
!    end do
!  !top edges
!  do k=1,mz
!     j1=0
!     do j=my-4,my-2
!        j1=j1+1
!       do i=1,mx
!           Denx_hbd2(i,j,k)=Denx_hbd2(i,j,k)+buff_tDenx_hbd2(i,j1,k)
!           Deny_hbd2(i,j,k)=Deny_hbd2(i,j,k)+buff_tDeny_hbd2(i,j1,k)
!           Denz_hbd2(i,j,k)=Denz_hbd2(i,j,k)+buff_tDenz_hbd2(i,j1,k)
!        end do
!    end do
! end do
! end subroutine AddImportedDensityZX_Hybrid2
!
! subroutine FoldInDensity_Hybrid2
! #ifndef twoD
!           call AddDensityLayerXY_Hybrid2(mz,5)
!           call AddDensityLayerXY_Hybrid2(mz-1,4)
!           call AddDensityLayerXY_Hybrid2(mz-2,3)
!           call AddDensityLayerXY_Hybrid2(1,mz-4)
!           call AddDensityLayerXY_Hybrid2(2,mz-3)
! #endif
! end subroutine FoldInDensity_Hybrid2
! subroutine AddDensityLayerXY_Hybrid2(zs,zd)
!      integer:: zs,zd,i,j
!      do j=1,my
!           do i=1,mx
!                Denx_hbd2(i,j,zd)=Denx_hbd2(i,j,zd)+Denx_hbd2(i,j,zs)
!                Deny_hbd2(i,j,zd)=Deny_hbd2(i,j,zd)+Deny_hbd2(i,j,zs)
!                Denz_hbd2(i,j,zd)=Denz_hbd2(i,j,zd)+Denz_hbd2(i,j,zs)
!           end do
!      end do
! end subroutine AddDensityLayerXY_Hybrid2
!
! subroutine CheckZeroDensity_Hybrid2
!      integer :: i,j,k
!      do k=1,1
!           do j=3,my-3
!                do i=3,mx-3
!                     !if(isnan(Denx_hbd2(i,j,k))) print*,'Denx is', Denx_hbd2(i,j,k), 'at',proc
!                     !if(isnan(Deny_hbd2(i,j,k))) print*,'Deny is', Deny_hbd2(i,j,k), 'at',proc
!                     !if(isnan(Denz_hbd2(i,j,k))) print*,'Denz is', Denz_hbd2(i,j,k), 'at',proc
!                     if(Deny_hbd2(i,j,k).eq.0) print*,'Deny is', Deny_hbd2(i,j,k), 'at',proc,'in',i,j,k
!
!                end do
!           end do
!      end do
! end subroutine CheckZeroDensity_Hybrid2
!
! !subroutines to update current on the outer edges since they are needed to compute electric field
! subroutine UpdateCurrentOuterEdges_Hybrid2
!     call ExchangeYZEdgeCurrent1(Jx)
!     call ExchangeYZEdgeCurrent1(Jy)
!     call ExchangeYZEdgeCurrent1(Jz)
!
!      call ExchangeZXEdgeCurrent1(Jx)
!      call ExchangeZXEdgeCurrent1(Jy)
!      call ExchangeZXEdgeCurrent1(Jz)
! #ifndef twoD
!     call CopyCurrentLayerXY(3,mz-2,Jx)
!     call CopyCurrentLayerXY(mz-3,2,Jx)
!     call CopyCurrentLayerXY(3,mz-2,Jy)
!     call CopyCurrentLayerXY(mz-3,2,Jy)
!     call CopyCurrentLayerXY(3,mz-2,Jz)
!     call CopyCurrentLayerXY(mz-3,2,Jz)
! #endif
!
! end subroutine UpdateCurrentOuterEdges_Hybrid2
!
!
!
! #endif

!------------------------------------------------------------------------------------------------------
!The following subroutines are useful for the Hybrid3 method 
!------------------------------------------------------------------------------------------------------     


! #ifdef Hybrid3
!
! subroutine UpdateEfieldHybrid3
!      real(psn) :: thisBx,thisBy,thisBz
!      real(psn) :: curlBx,curlBy,curlBz
!      integer   :: i1,j1,k1,i2,j2,k2
!
!
!      i1=3
!      i2=mx-3
!      j1=3
!      j2=my-3
! #ifndef twoD
!      k1=3
!      k2=mz-3
! #else
!      k1=1
!      k2=1
! #endif
!     do k=k1,k2
!         do j=j1,j2
!              do i=i1,i2
! #ifndef twoD
!                !magnetic field at (i,j,k)
!                thisBx=( Bx(i,j,k)+Bx(i,j,k-1)+Bx(i,j-1,k)+Bx(i,j-1,k-1) )*0.25 +Bx_ext0
!                thisBy=( By(i,j,k)+By(i,j,k-1)+By(i-1,j,k)+By(i-1,j,k-1) )*0.25 +By_ext0
!                thisBz=( Bz(i,j,k)+Bz(i,j-1,k)+Bz(i-1,j,k)+Bz(i-1,j-1,k) )*0.25 +Bz_ext0
!
!                !curl of magnetic field at (i,j,k)
!                curlBx=( (Bz(i,j,k)  -Bz(i,j-1,k))   - (By(i,j,k)  -By(i,j,k-1))   &
!                        +(Bz(i-1,j,k)-Bz(i-1,j-1,k)) - (By(i-1,j,k)-By(i-1,j,k-1)) )*0.5
!                curlBy=( (Bx(i,j,k)  -Bx(i,j,k-1))   - (Bz(i,j,k)  -Bz(i-1,j,k))   &
!                        +(Bx(i,j-1,k)-Bx(i,j-1,k-1)) - (Bz(i,j-1,k)-Bz(i-1,j-1,k)) )*0.5
!                curlBz=( (By(i,j,k)  -By(i-1,j,k))   - (Bx(i,j,k)  -Bx(i,j-1,k))   &
!                        +(By(i,j,k-1)-By(i-1,j,k-1)) - (Bx(i,j,k-1)-Bx(i,j-1,k-1)) )*0.5
!
!                Ex(i,j,k)= (-Jy(i,j,k)+curlBy)*thisBz - (-Jz(i,j,k)+curlBz)*thisBy
!                Ey(i,j,k)= (-Jz(i,j,k)+curlBz)*thisBx - (-Jx(i,j,k)+curlBx)*thisBz
!                Ez(i,j,k)= (-Jx(i,j,k)+curlBx)*thisBy - (-Jy(i,j,k)+curlBy)*thisBx
! #else
!                !magnetic field at (i,j,k)
!                thisBx=( Bx(i,j,k)+Bx(i,j-1,k) )*0.5 + Bx_ext0
!                thisBy=( By(i,j,k)+By(i-1,j,k) )*0.5 + By_ext0
!                thisBz=( Bz(i,j,k)+Bz(i,j-1,k)+Bz(i-1,j,k)+Bz(i-1,j-1,k) )*0.25 +Bz_ext0
!
!                !curl of magnetic field at (i,j,k)
!                curlBx=( (Bz(i,j,k)  -Bz(i,j-1,k))   &
!                        +(Bz(i-1,j,k)-Bz(i-1,j-1,k)) )*0.5
!                curlBy=( - (Bz(i,j,k)  -Bz(i-1,j,k))   &
!                         - (Bz(i,j-1,k)-Bz(i-1,j-1,k)) )*0.5
!                curlBz=( (By(i,j,k)  -By(i-1,j,k))   - (Bx(i,j,k)  -Bx(i,j-1,k)) )
!
!                gEx_hbd3(i,j,k)=( (-Jy(i,j,k)+c*curlBy)*thisBz - (-Jz(i,j,k)+c*curlBz)*thisBy  )/(c*Den_hbd3(i,j,k))
!                gEy_hbd3(i,j,k)=( (-Jz(i,j,k)+c*curlBz)*thisBx - (-Jx(i,j,k)+c*curlBx)*thisBz  )/(c*Den_hbd3(i,j,k))
!                gEz_hbd3(i,j,k)=( (-Jx(i,j,k)+c*curlBx)*thisBy - (-Jy(i,j,k)+c*curlBy)*thisBx  )/(c*Den_hbd3(i,j,k))
!
!                if(thisBy.ne.0) print*,'by is not zero',thisBy
!                !Den_hbd3(i,j,k)=thisBy
!
! #endif
!               end do
!           end do
!      end do
!
! !update the outer edges of the electric field gEx which used later for computing Ex, edges on other proc are communicated
! #ifndef twoD
!         gEx_hbd3(:,:,mz-2)=gEx_hbd3(:,:,3)
!         gEy_hbd3(:,:,mz-2)=gEy_hbd3(:,:,3)
!         gEz_hbd3(:,:,mz-2)=gEz_hbd3(:,:,3)
! #endif
!
! end subroutine UpdateEfieldHybrid3
!
! subroutine UpdateDensityAllEdgesHybrid3
!      call ExchangeYZEdgeDensity_Hybrid3
!     call AddImportedDensityYZ_Hybrid3
!
!      call ExchangeZXEdgeDensity_Hybrid3
!      call AddImportedDensityZX_Hybrid3
!
! #ifndef twoD
!     Den_hbd3(3:5,:,:)=Den_hbd3(3:5)+Den_hbd3(mx-4:mx-2,:,:)
!      Den_hbd3(mx-4:mx,:,:)=Den_hbd3(3:5,:,:)
! #endif
! end subroutine UpdateDensityAllEdgesHybrid3
!
!
! subroutine UpdateEdgeEfieldHybrid3
!      integer   :: i1,j1,k1,i2,j2,k2
!      i1=3
!      i2=mx-3
!      j1=3
!      j2=my-3
! #ifndef twoD
!      k1=3
!      k2=mz-3
! #else
!      k1=1
!      k2=1
! #endif
! !Now compute the electric field at the Grid points
!      do k=k1,k2
!          do j=j1,j2
!               do i=i1,i2
!                Ex(i,j,k)=( gEx_hbd3(i,j,k) + gEx_hbd3(i+1,j,k) )*0.5
!                Ey(i,j,k)=( gEy_hbd3(i,j,k) + gEy_hbd3(i,j+1,k) )*0.5
! #ifndef twoD
!                Ez(i,j,k)=( gEz_hbd3(i,j,k) + gEz_hbd3(i,j,k+1) )*0.5
! #else
!             Ez(i,j,k)=  gEz_hbd3(i,j,k)
! #endif
!                end do
!           end do
!      end do
!
! end subroutine UpdateEdgeEfieldHybrid3
!
! subroutine AddImportedDensityYZ_Hybrid3
!      integer:: i,j,k,i1
!    do k=1,mz
!        do j=1,my
!            do i=3,5
!                Den_hbd3(i,j,k)=Den_hbd3(i,j,k)+buff_lDen_hbd3(i-2,j,k)
!             end do
!        end do
!    end do
!  !right edges
!  do k=1,mz
!     do j=1,my
!        i1=0
!       do i=mx-4,mx-2
!        i1=i1+1
!           Den_hbd3(i,j,k)=Den_hbd3(i,j,k)+buff_rDen_hbd3(i1,j,k)
!        end do
!    end do
! end do
! end subroutine AddImportedDensityYZ_Hybrid3
! subroutine AddImportedDensityZX_Hybrid3
!      integer:: i,j,k,j1
!    do k=1,mz
!        do j=3,5
!            do i=1,mx
!                Den_hbd3(i,j,k)=Den_hbd3(i,j,k)+buff_bDen_hbd3(i,j-2,k)
!             end do
!        end do
!    end do
!  !top edges
!  do k=1,mz
!     j1=0
!     do j=my-4,my-2
!        j1=j1+1
!       do i=1,mx
!           Den_hbd3(i,j,k)=Den_hbd3(i,j,k)+buff_tDen_hbd3(i,j1,k)
!        end do
!    end do
! end do
! end subroutine AddImportedDensityZX_Hybrid3
! subroutine smoothen_densityHybrid3
!      integer :: ni
!      do ni=1,curr_filter
!           call avg_curr(Den_hbd3)
!      end do
! end subroutine smoothen_densityHybrid3
!
!
!
! #endif



end module fields