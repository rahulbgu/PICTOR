! This module is designed to facilitate transfer of particles between MPI proc.
module particles 
     use parameters 
     use vars
     use memory
     use comm_fldprtl
#ifdef cyl
     use cyl_comm_fldprtl
#endif	  
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
          call LoadTestPrtlOutliers !Test Particles are also loaded in the same array as particles
          ltpcross=lcross-lpcross
          rtpcross=rcross-rpcross
          ttpcross=tcross-tpcross
          btpcross=bcross-bpcross
          utpcross=ucross-upcross
          dtpcross=dcross-dpcross
		  ntp=ntp-(ltpcross+rtpcross+ttpcross+btpcross+utpcross+dtpcross)
          
          call StopTimer(9)     
          call SendRecvPrtlSize !test particles are transferred in the same array 
#ifdef cyl
          if(inc_axis) call ExchangePrtlAxis
#endif		  
          call UpdateTransferInSize
          call SendRecvPrtl
		  call UpdateTransferOutSize
     end subroutine ExchangePrtl
	 
	
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
	  
	 
end module particles