! This module is designed to facilitate transfer of particles between MPI proc.
module particles 
     use parameters 
     use vars
     use memory
     use comm_fldprtl 
	 use loadprtlout
     implicit none
contains
           
     subroutine AppendParticles
          integer :: i 

          !np=np+(linp_count+rinp_count+binp_count+tinp_count+uinp_count+dinp_count) !p_count= count(no.) of p (particles) that arrived from another domain 
		  !if(np.gt.prtl_arr_size) call ReshapePrtlArr(int(1.1*np+100))

          do i=1,rinp_count
			   rinp(i)%x=rinp(i)%x+xlen
			   call InsertNewPrtl(rinp(i)%x,rinp(i)%y,rinp(i)%z,rinp(i)%u,rinp(i)%v,rinp(i)%w,rinp(i)%q,rinp(i)%tag,rinp(i)%flv,rinp(i)%var1)	
          end do 
          do i=1,linp_count
               linp(i)%x=linp(i)%x-aint(linp(i)%x)+xmin
			   call InsertNewPrtl(linp(i)%x,linp(i)%y,linp(i)%z,linp(i)%u,linp(i)%v,linp(i)%w,linp(i)%q,linp(i)%tag,linp(i)%flv,linp(i)%var1)	
          end do
          do i=1,tinp_count
               tinp(i)%y=tinp(i)%y+ylen
			   call InsertNewPrtl(tinp(i)%x,tinp(i)%y,tinp(i)%z,tinp(i)%u,tinp(i)%v,tinp(i)%w,tinp(i)%q,tinp(i)%tag,tinp(i)%flv,tinp(i)%var1)	
          end do
          do i=1,binp_count
               binp(i)%y=binp(i)%y-aint(binp(i)%y)+ymin
			   call InsertNewPrtl(binp(i)%x,binp(i)%y,binp(i)%z,binp(i)%u,binp(i)%v,binp(i)%w,binp(i)%q,binp(i)%tag,binp(i)%flv,binp(i)%var1)	
          end do
#ifndef twoD
          do i=1,uinp_count
               uinp(i)%z=uinp(i)%z+zlen
			   call InsertNewPrtl(uinp(i)%x,uinp(i)%y,uinp(i)%z,uinp(i)%u,uinp(i)%v,uinp(i)%w,uinp(i)%q,uinp(i)%tag,uinp(i)%flv,uinp(i)%var1)	
          end do
          do i=1,dinp_count
               dinp(i)%z=dinp(i)%z-aint(dinp(i)%z)+zmin
			   call InsertNewPrtl(dinp(i)%x,dinp(i)%y,dinp(i)%z,dinp(i)%u,dinp(i)%v,dinp(i)%w,dinp(i)%q,dinp(i)%tag,dinp(i)%flv,dinp(i)%var1)	
		  end do       
#endif          
     end subroutine AppendParticles
     subroutine AppendTestParticles
          integer :: i 
          do i=rinp_count+1,rinp_count+rintp_count
			   rinp(i)%x=rinp(i)%x+xlen
			   call InsertNewTestPrtl(rinp(i)%x,rinp(i)%y,rinp(i)%z,rinp(i)%u,rinp(i)%v,rinp(i)%w,rinp(i)%q,rinp(i)%tag,rinp(i)%flv,rinp(i)%var1)	
          end do 
          do i=linp_count+1,linp_count+lintp_count
               linp(i)%x=linp(i)%x-aint(linp(i)%x)+xmin
			   call InsertNewTestPrtl(linp(i)%x,linp(i)%y,linp(i)%z,linp(i)%u,linp(i)%v,linp(i)%w,linp(i)%q,linp(i)%tag,linp(i)%flv,linp(i)%var1)	
          end do
          do i=tinp_count+1,tinp_count+tintp_count
               tinp(i)%y=tinp(i)%y+ylen
			   call InsertNewTestPrtl(tinp(i)%x,tinp(i)%y,tinp(i)%z,tinp(i)%u,tinp(i)%v,tinp(i)%w,tinp(i)%q,tinp(i)%tag,tinp(i)%flv,tinp(i)%var1)	
          end do
          do i=binp_count+1,binp_count+bintp_count
               binp(i)%y=binp(i)%y-aint(binp(i)%y)+ymin
			   call InsertNewTestPrtl(binp(i)%x,binp(i)%y,binp(i)%z,binp(i)%u,binp(i)%v,binp(i)%w,binp(i)%q,binp(i)%tag,binp(i)%flv,binp(i)%var1)	
          end do
#ifndef twoD
          do i=uinp_count+1,uinp_count+uintp_count
               uinp(i)%z=uinp(i)%z+zlen
			   call InsertNewTestPrtl(uinp(i)%x,uinp(i)%y,uinp(i)%z,uinp(i)%u,uinp(i)%v,uinp(i)%w,uinp(i)%q,uinp(i)%tag,uinp(i)%flv,uinp(i)%var1)	
          end do
          do i=dinp_count+1,dinp_count+dintp_count
               dinp(i)%z=dinp(i)%z-aint(dinp(i)%z)+zmin
			   call InsertNewTestPrtl(dinp(i)%x,dinp(i)%y,dinp(i)%z,dinp(i)%u,dinp(i)%v,dinp(i)%w,dinp(i)%q,dinp(i)%tag,dinp(i)%flv,dinp(i)%var1)	
		  end do       
#endif          
     end subroutine AppendTestParticles
     
     subroutine ExchangePrtl
		  call StartTimer(9)
          lcross=0
          rcross=0
          tcross=0
          bcross=0
          ucross=0
          dcross=0
          call LoadPrtlOutliers
          lpcross=lcross
          rpcross=rcross
          tpcross=tcross
          bpcross=bcross
          upcross=ucross
          dpcross=dcross
		  np=np-(lpcross+rpcross+tpcross+bpcross+upcross+dpcross)
          call LoadTestPrtlOutliers
          ltpcross=lcross-lpcross
          rtpcross=rcross-rpcross
          ttpcross=tcross-tpcross
          btpcross=bcross-bpcross
          utpcross=ucross-upcross
          dtpcross=dcross-dpcross
		  ntp=ntp-(ltpcross+rtpcross+ttpcross+btpcross+utpcross+dtpcross)
          
          call StopTimer(9)     
          call SendRecvPrtlSize !test particles are transferred in the same array 
          call UpdateTransferInSize
          call SendRecvPrtl
     end subroutine ExchangePrtl
	 



     subroutine UpdateTransferInSize
            if(linp_size.lt.(linp_count+lintp_count)) call ReshapeTransferInArr(linp,linp_size,int((linp_count+lintp_count)*1.1+100))
            if(rinp_size.lt.(rinp_count+rintp_count)) call ReshapeTransferInArr(rinp,rinp_size,int((rinp_count+rintp_count)*1.1+100))
            if(tinp_size.lt.(tinp_count+tintp_count)) call ReshapeTransferInArr(tinp,tinp_size,int((tinp_count+tintp_count)*1.1+100))
            if(binp_size.lt.(binp_count+bintp_count)) call ReshapeTransferInArr(binp,binp_size,int((binp_count+bintp_count)*1.1+100))
#ifndef twoD
            if(uinp_size.lt.(uinp_count+uintp_count)) call ReshapeTransferInArr(uinp,uinp_size,int((uinp_count+uintp_count)*1.1+100))
            if(dinp_size.lt.(dinp_count+dintp_count)) call ReshapeTransferInArr(dinp,dinp_size,int((dinp_count+dintp_count)*1.1+100))
#endif             
    end subroutine UpdateTransferInSize 
	
	![Not In use :: outdated]
     subroutine RegulatePrtlMemory
          integer :: hist_max
          
          np_history(hist_ind)=np
          hist_ind=hist_ind+1
          if(hist_ind.gt.size_hist_bin_size) hist_ind=1
          
          !reduce the size of particle array to optimise the performance
          if(modulo(t,size_hist_bin_size).eq.0) then
             hist_max=maxval(np_history)
             if(hist_max.lt.int(prtl_arr_size*0.9)) call ReshapePrtlArr(int(hist_max*1.1))
         end if                        
     end subroutine RegulatePrtlMemory
	 
 
!-----------
! The following subroutines were desgined to work with single tthread 
! New subroutines supposrts openmp 
! REshaping of transfer arr size is removed  
 
!
!  #ifdef twoD
!       subroutine LoadPrtlOutliers
!  		 implicit none
!            integer:: i
!
!            do i=1,used_prtl_arr_size
!                 if(qp(i).eq.0) cycle
!                 if(xp(i).lt.xmin) then
!                      if((yp(i).le.(ymax+xmin-xp(i))).and.(yp(i).ge.(ymin-xmin+xp(i)))) then
!                           lcross=lcross+1
!                           if(lcross.gt.loutp_size) call ReshapeTransferOutArr(loutp,loutp_size,int(loutp_size*1.1+10))
!                           call LoadPrtl(loutp,loutp_size,lcross,i)
!                           call DeletePrtl(i)
!                      end if
!                 else if(xp(i).gt.xmax) then
!                      if((yp(i).le.(ymax+xp(i)-xmax)).and.(yp(i).ge.(ymin+xmax-xp(i)))) then
!                           rcross=rcross+1
!                           if(rcross.gt.routp_size) call ReshapeTransferOutArr(routp,routp_size,int(routp_size*1.1+10))
!  						 call LoadPrtl(routp,routp_size,rcross,i)
!                           call DeletePrtl(i)
!                  end if
!                 end if
!                 if(yp(i).lt.ymin) then
!                      if((xp(i).gt.(xmin-ymin+yp(i))).and.(xp(i).lt.(xmax+ymin-yp(i)))) then
!                           bcross=bcross+1
!                           if(bcross.gt.boutp_size) call ReshapeTransferOutArr(boutp,boutp_size,int(boutp_size*1.1+10))
!  						 call LoadPrtl(boutp,boutp_size,bcross,i)
!                           call DeletePrtl(i)
!                      end if
!                 else if(yp(i).gt.ymax) then
!                      if((xp(i).gt.(xmin-yp(i)+ymax)).and.(xp(i).lt.(xmax+yp(i)-ymax))) then
!                           tcross=tcross+1
!                           if(tcross.gt.toutp_size) call ReshapeTransferOutArr(toutp,toutp_size,int(toutp_size*1.1+10))
!  						 call LoadPrtl(toutp,toutp_size,tcross,i)
!                           call DeletePrtl(i)
!                      end if
!                 end if
!           end do
!
!       end subroutine LoadPrtlOutliers
!
!       subroutine LoadTestPrtlOutliers
!  		 implicit none
!            integer:: i
!            do i=1,used_test_prtl_arr_size
!                 if(qtp(i).eq.0) cycle
!                 if(xtp(i).lt.xmin) then
!                      if((ytp(i).le.(ymax+xmin-xtp(i))).and.(ytp(i).ge.(ymin-xmin+xtp(i)))) then
!                           lcross=lcross+1
!                           if(lcross.gt.loutp_size) call ReshapeTransferOutArr(loutp,loutp_size,int(loutp_size*1.1+10))
!                           call LoadTestPrtl(loutp,loutp_size,lcross,i)
!                           call DeleteTestPrtl(i)
!                      end if
!                 else if(xtp(i).gt.xmax) then
!                      if((ytp(i).le.(ymax+xtp(i)-xmax)).and.(ytp(i).ge.(ymin+xmax-xtp(i)))) then
!                           rcross=rcross+1
!                           if(rcross.gt.routp_size) call ReshapeTransferOutArr(routp,routp_size,int(routp_size*1.1+10))
!  						 call LoadTestPrtl(routp,routp_size,rcross,i)
!                           call DeleteTestPrtl(i)
!                  end if
!                 end if
!                 if(ytp(i).lt.ymin) then
!                      if((xtp(i).gt.(xmin-ymin+ytp(i))).and.(xtp(i).lt.(xmax+ymin-ytp(i)))) then
!                           bcross=bcross+1
!                           if(bcross.gt.boutp_size) call ReshapeTransferOutArr(boutp,boutp_size,int(boutp_size*1.1+10))
!  						 call LoadTestPrtl(boutp,boutp_size,bcross,i)
!                           call DeleteTestPrtl(i)
!                      end if
!                 else if(ytp(i).gt.ymax) then
!                      if((xtp(i).gt.(xmin-ytp(i)+ymax)).and.(xtp(i).lt.(xmax+ytp(i)-ymax))) then
!                           tcross=tcross+1
!                           if(tcross.gt.toutp_size) call ReshapeTransferOutArr(toutp,toutp_size,int(toutp_size*1.1+10))
!  						 call LoadTestPrtl(toutp,toutp_size,tcross,i)
!                           call DeleteTestPrtl(i)
!                      end if
!                 end if
!           end do
!       end subroutine LoadTestPrtlOutliers
!
!
!
!
!  #else
!       subroutine LoadPrtlOutliers
!  		 implicit none
!            integer:: i
!
!            !print*,'at',proc,'pfree in',pfree_in,'pfree out',pfree_out
!            do i=1,used_prtl_arr_size
!                 if(qp(i).eq.0) cycle
!                 if(xp(i).lt.xmin) then
!                      if((yp(i).le.(ymax+xmin-xp(i))).and.(yp(i).ge.(ymin-xmin+xp(i)))) then
!                           if(zp(i).le.(zmax+xmin-xp(i)).and.(zp(i).ge.(zmin-xmin+xp(i)))) then
!                               lcross=lcross+1
!                               if(lcross.gt.loutp_size) call ReshapeTransferOutArr(loutp,loutp_size,int(loutp_size*1.1+10))
!                               call LoadPrtl(loutp,loutp_size,lcross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 else if(xp(i).gt.xmax) then
!                      if((yp(i).le.(ymax+xp(i)-xmax)).and.(yp(i).ge.(ymin+xmax-xp(i)))) then
!                           if((zp(i).le.(zmax+xp(i)-xmax)).and.(zp(i).ge.(zmin+xmax-xp(i)))) then
!                               rcross=rcross+1
!                               if(rcross.gt.routp_size) call ReshapeTransferOutArr(routp,routp_size,int(routp_size*1.1+10))
!  							 call LoadPrtl(routp,routp_size,rcross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 end if
!                 if(yp(i).lt.ymin) then
!                      if((xp(i).gt.(xmin-ymin+yp(i))).and.(xp(i).lt.(xmax+ymin-yp(i)))) then
!                           if((zp(i).ge.(zmin-ymin+yp(i))).and.(zp(i).le.(zmax+ymin-yp(i)))) then
!                               bcross=bcross+1
!                               if(bcross.gt.boutp_size) call ReshapeTransferOutArr(boutp,boutp_size,int(boutp_size*1.1+10))
!  							 call LoadPrtl(boutp,boutp_size,bcross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 else if(yp(i).gt.ymax) then
!                      if((xp(i).gt.(xmin-yp(i)+ymax)).and.(xp(i).lt.(xmax+yp(i)-ymax))) then
!                           if((zp(i).ge.(zmin-yp(i)+ymax)).and.(zp(i).le.(zmax+yp(i)-ymax))) then
!                               tcross=tcross+1
!                               if(tcross.gt.toutp_size) call ReshapeTransferOutArr(toutp,toutp_size,int(toutp_size*1.1+10))
!  							 call LoadPrtl(toutp,toutp_size,tcross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 end if
!
!
!                 if(zp(i).lt.zmin) then
!                      if((xp(i).gt.(xmin-zmin+zp(i))).and.(xp(i).lt.(xmax+zmin-zp(i)))) then
!                           if((yp(i).gt.(ymin-zmin+zp(i))).and.(yp(i).lt.(ymax+zmin-zp(i)))) then
!                               dcross=dcross+1
!                               if(dcross.gt.doutp_size) call ReshapeTransferOutArr(doutp,doutp_size,int(doutp_size*1.1+10))
!  							 call LoadPrtl(doutp,doutp_size,dcross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 else if(zp(i).gt.zmax) then
!                      if((xp(i).gt.(xmin-zp(i)+zmax)).and.(xp(i).lt.(xmax+zp(i)-zmax))) then
!                           if((yp(i).gt.(ymin-zp(i)+zmax)).and.(yp(i).lt.(ymax+zp(i)-zmax))) then
!                               ucross=ucross+1
!                               if(ucross.gt.uoutp_size) call ReshapeTransferOutArr(uoutp,uoutp_size,int(uoutp_size*1.1+10))
!  							 call LoadPrtl(uoutp,uoutp_size,ucross,i)
!                               call DeletePrtl(i)
!                          end if
!                      end if
!                 end if
!           end do
!       end subroutine LoadPrtlOutliers
!
!       subroutine LoadTestPrtlOutliers
!  		 implicit none
!            integer:: i
!
!            do i=1,used_test_prtl_arr_size
!                 if(qtp(i).eq.0) cycle
!                 if(xtp(i).lt.xmin) then
!                      if((ytp(i).le.(ymax+xmin-xtp(i))).and.(ytp(i).ge.(ymin-xmin+xtp(i)))) then
!                           if(ztp(i).le.(zmax+xmin-xtp(i)).and.(ztp(i).ge.(zmin-xmin+xtp(i)))) then
!                               lcross=lcross+1
!                               if(lcross.gt.loutp_size) call ReshapeTransferOutArr(loutp,loutp_size,int(loutp_size*1.1+10))
!                               call LoadTestPrtl(loutp,loutp_size,lcross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 else if(xtp(i).gt.xmax) then
!                      if((ytp(i).le.(ymax+xtp(i)-xmax)).and.(ytp(i).ge.(ymin+xmax-xtp(i)))) then
!                           if((ztp(i).le.(zmax+xtp(i)-xmax)).and.(ztp(i).ge.(zmin+xmax-xtp(i)))) then
!                               rcross=rcross+1
!                               if(rcross.gt.routp_size) call ReshapeTransferOutArr(routp,routp_size,int(routp_size*1.1+10))
!  							 call LoadTestPrtl(routp,routp_size,rcross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 end if
!                 if(ytp(i).lt.ymin) then
!                      if((xtp(i).gt.(xmin-ymin+ytp(i))).and.(xtp(i).lt.(xmax+ymin-ytp(i)))) then
!                           if((ztp(i).ge.(zmin-ymin+ytp(i))).and.(ztp(i).le.(zmax+ymin-ytp(i)))) then
!                               bcross=bcross+1
!                               if(bcross.gt.boutp_size) call ReshapeTransferOutArr(boutp,boutp_size,int(boutp_size*1.1+10))
!  							 call LoadTestPrtl(boutp,boutp_size,bcross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 else if(ytp(i).gt.ymax) then
!                      if((xtp(i).gt.(xmin-ytp(i)+ymax)).and.(xtp(i).lt.(xmax+ytp(i)-ymax))) then
!                           if((ztp(i).ge.(zmin-ytp(i)+ymax)).and.(ztp(i).le.(zmax+ytp(i)-ymax))) then
!                               tcross=tcross+1
!                               if(tcross.gt.toutp_size) call ReshapeTransferOutArr(toutp,toutp_size,int(toutp_size*1.1+10))
!  							 call LoadTestPrtl(toutp,toutp_size,tcross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 end if
!
!
!                 if(ztp(i).lt.zmin) then
!                      if((xtp(i).gt.(xmin-zmin+ztp(i))).and.(xtp(i).lt.(xmax+zmin-ztp(i)))) then
!                           if((ytp(i).gt.(ymin-zmin+ztp(i))).and.(ytp(i).lt.(ymax+zmin-ztp(i)))) then
!                               dcross=dcross+1
!                               if(dcross.gt.doutp_size) call ReshapeTransferOutArr(doutp,doutp_size,int(doutp_size*1.1+10))
!  							 call LoadTestPrtl(doutp,doutp_size,dcross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 else if(ztp(i).gt.zmax) then
!                      if((xtp(i).gt.(xmin-ztp(i)+zmax)).and.(xtp(i).lt.(xmax+ztp(i)-zmax))) then
!                           if((ytp(i).gt.(ymin-ztp(i)+zmax)).and.(ytp(i).lt.(ymax+ztp(i)-zmax))) then
!                               ucross=ucross+1
!                               if(ucross.gt.uoutp_size) call ReshapeTransferOutArr(uoutp,uoutp_size,int(uoutp_size*1.1+10))
!  							 call LoadTestPrtl(uoutp,uoutp_size,ucross,i)
!                               call DeleteTestPrtl(i)
!                          end if
!                      end if
!                 end if
!           end do
!       end subroutine LoadTestPrtlOutliers
!  #endif
	 
	 
	 
end module particles