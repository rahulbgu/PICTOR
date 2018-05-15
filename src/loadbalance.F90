module loadbalance
     use parameters
     use vars
     use memory
	 use movdep
     use comm_loadbalance
     use fields, only : ExchangeYZEdgeField, ExchangeZXEdgeField, ExchangeXYEdgeField !use of theses subroutine can be avoided  
     implicit none 
     integer :: mx_new,my_new,mz_new !new size of the fld arrays, used in load redistribution 
#ifdef twoD	 
     real, dimension(0:nSubDomainsX*nSubDomainsY-1) :: TotalTime
#else 
     real, dimension(0:nSubDomainsX*nSubDomainsY*nSubDomainsZ-1) :: TotalTime
#endif	 
	 real, dimension(0:nSubDomainsX-1) :: TotalTimeX 
      
contains 
!------------------------------------------------------------------------------
!The following subroutines are helpful in the case of inhomogenous plasma,such as shock, when physical domain on proc is allowed to change 
!------------------------------------------------------------------------------

!Important note::: "The following subroutine is very outdated"
! It was written for 2D domain decomposition and now needs to be updated for 3D domain decomposition and new schemes



!  subroutine LoadBalanceHomogeneous
!       integer :: i,j,l,m,i0,j0,i1,i2
!       real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1) :: exectime_grid
!       integer,dimension(0:nSubDomainsX-1,0:nSubDomainsY-1) :: proc_swap
!       real, dimension(0:nSubDomainsX-1) :: mean_TotalTime_column
!       real :: time_temp
!       integer, dimension(nproc) :: slow_proc,replacement_proc
!       integer :: loadbalance_profiling_ind
!       logical :: xborders_changed
!       save loadbalance_profiling_ind
!       data loadbalance_profiling_ind /1/
!       xborders_changed=.false.
!       !profile the performance history of all processors before making any decison about balancing the laod
!       if(modulo(t,loadbalance_profiling_period).eq.0) then
!            if(loadbalance_profiling_ind.gt.100) loadbalance_profiling_ind=1
!            hist_TotalTime(loadbalance_profiling_ind)=real(exec_time(4))/xlen
!            loadbalance_profiling_ind=loadbalance_profiling_ind+1
!            mean_TotalTime=average(hist_TotalTime,100)
!       end if
!
!        if(modulo(t,domain_change_period).eq.0) then
!            exectime_grid=0
!            exectime_grid(procxind(proc),procyind(proc))=mean_TotalTime
!            call BcastExecTimeAllHomogeneous(exectime_grid)
!
!          !first calculate the average execution time of all proc in a column (i=const.)
!            mean_TotalTime_column=0
!            do i=0,nSubDomainsX-1
!                do j=0,nSubDomainsY-1
!                     mean_TotalTime_column(i)=mean_TotalTime_column(i)+exectime_grid(i,j)
!                end do
!                    mean_TotalTime_column(i)=mean_TotalTime_column(i)/nSubDomainsY
!            end do
!            !print*, 'column average', mean_TotalTime_column
!
!
!            ! prepare a list of underperforming proc
!            l=0
!            do i=0,nSubDomainsX-1
!                 do j=0,nSubDomainsY-1
!                      if(exectime_grid(i,j).gt.mean_TotalTime_column(i)*1.1) then
!                           l=l+1
!                           slow_proc(l)=proc_grid(i,j)
!                      end if
!                 end do
!            end do
!            if((l.gt.0).and.(proc.eq.0)) print *,'Proc found to slow down the performance are:',slow_proc(1:l)
! !            if(proc.eq.0) then
! !            do i=0,nSubDomainsX-1
! !                 do j=0,nSubDomainsY-1
! !                  print*,'AT PROC',proc_grid(i,j), 'time',exectime_grid(i,j)
! !                 end do
! !            end do
! !           end if
!
!            !find a well performing proc that would replace underperforming proc
!          proc_swap=0
!            do m=1,l
!                 replacement_proc(m)=-1
!                 do i=0,procxind(slow_proc(m))-1
!                      do j=0,nSubDomainsY-1
!                          if(proc_swap(i,j).eq.0) then
!                                if(exectime_grid(i,j)*1.05.lt.exectime_grid(procxind(slow_proc(m)),procyind(slow_proc(m)))) then
!                                     proc_swap(i,j)=1
!                                     proc_swap(procxind(slow_proc(m)),procyind(slow_proc(m)))=1
!                                     replacement_proc(m)=proc_grid(i,j)
!                                end if
!                           end if
!                           if(replacement_proc(m).ge.0) exit
!                      end do
!                      if(replacement_proc(m).ge.0) exit
!                 end do
!            end do
!
!            do m=1,l
!                 if(replacement_proc(m).ge.0) then
!                      proc_grid(procxind(replacement_proc(m)),procyind(replacement_proc(m)))=slow_proc(m)
!                      proc_grid(procxind(slow_proc(m)),procyind(slow_proc(m)))=replacement_proc(m)
!                      time_temp=exectime_grid(procxind(replacement_proc(m)),procyind(replacement_proc(m)))
!                      exectime_grid(procxind(replacement_proc(m)),procyind(replacement_proc(m)))=exectime_grid(procxind(slow_proc(m)),procyind(slow_proc(m)))
!                  exectime_grid(procxind(slow_proc(m)),procyind(slow_proc(m)))=time_temp
!
!                      i1=procxind(replacement_proc(m))
!                      i2=procyind(replacement_proc(m))
!                      procxind(replacement_proc(m))=procxind(slow_proc(m))
!                      procyind(replacement_proc(m))=procyind(slow_proc(m))
!                      procxind(slow_proc(m))=i1
!                      procyind(slow_proc(m))=i2
!                      call SwapDomain(slow_proc(m),replacement_proc(m))
!                      if(proc.eq.0) print*,'Swapped slow proc', slow_proc(m),'With',replacement_proc(m)
!                end if
!            end do
!            !if(l.gt.0) print*,'replacement proc is',m,replacement_proc
!
!
!
!            !put the exetime of the slowest proc in mean_TotalTime_column
!            mean_TotalTime_column=0
!            do i=0,nSubDomainsX-1
!                 do j=0,nSubDomainsY-1
!                      if(exectime_grid(i,j).gt.mean_TotalTime_column(i)) mean_TotalTime_column(i)=exectime_grid(i,j)
!                 end do
!            end do
!
!            !Adjust the x-boundaries to redistribute load
!            time_temp=0
!            do i=0,nSubDomainsX-1
!                 time_temp=time_temp+1/mean_TotalTime_column(i)
!            end do
!            time_temp=(xborders(nSubDomainsX)-xborders(0))/time_temp
!
!            xborders_new(0)=xborders(0)
!            do i=0,nSubDomainsX-2
!                 !xborders_new(i+1)=xborders_new(i)+max(1,int(time_temp*(xborders(i+1)-xborders(i))/mean_TotalTime_column(i)))
!                 xborders_new(i+1)=xborders_new(i)+max(4,fsave_ratio*int(time_temp/(mean_TotalTime_column(i)*fsave_ratio)))
!                 xborders_new(i+1)=min(xborders_new(i+1),xborders(nSubDomainsX)-max(fsave_ratio,4)*(nSubDomainsX-i-1)) !to set maximum limit on domain enlargement
!                 !the above scheme is due to the way save data is written, must be modified
!
!                 if(abs(xborders_new(i+1)-xborders(i+1)).gt.0) xborders_changed=.true.
!            end do
!            xborders_new(nSubDomainsX)=xborders(nSubDomainsX)
!
!            if(xborders_changed) then
!                 if(proc.eq.0) then
!                     print*,'Changing the Xborders'
!                      print*,'New Borders',xborders_new
!                      print*,'Old Borders',xborders
!                      print*, 'proc_grid is ',proc_grid
!                 end if
!                 call SetNewXBorders
!            end if
!
!        end if
!
!  end subroutine LoadBalanceHomogeneous

!-------------------------------------------------------------------------------------
! The following load-balancing scheme adjusts y-borders only
!-------------------------------------------------------------------------------------
 subroutine LoadBalanceY
      integer :: i,j,k,dYdomain 
#ifdef twoD
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: exectime_grid
#else 
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: exectime_grid
#endif 	  

      real, dimension(0:nSubDomainsY-1) :: Max_TotalTime_column
      real :: time_temp
      logical :: yborders_changed
 
      yborders_changed=.false.

      if(modulo(t,domain_change_period).eq.0) then
           exectime_grid=0
           dYdomain=yborders(procyind(proc)+1)-yborders(procyind(proc))
		   exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/dYdomain
           call BcastExecTimeAllHomogeneous(exectime_grid)
         !first calculate the average execution time of all proc in a column (i=const.)
           Max_TotalTime_column=0
           do i=0,nSubDomainsX-1
               do j=0,nSubDomainsY-1
#ifdef twoD
                   do k=0,0
#else				   
				   do k=0,nSubDomainsZ-1
#endif
                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(j)) Max_TotalTime_column(j)=exectime_grid(i,j,k)
				   end do 
               end do
           end do
          
           !Adjust the x-boundaries to redistribute load, and do not include last proc along x in load balancing 
           time_temp=0
           do j=0,nSubDomainsY-1
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(j)
           end do
           time_temp=(yborders(nSubDomainsY)-yborders(0))/time_temp
		   
           yborders_new(0)=yborders(0)
           do j=0,nSubDomainsY-2
                   !xborders_new(i+1)=xborders_new(i)+max(1,int(time_temp*(xborders(i+1)-xborders(i))/mean_TotalTime_column(i)))
                   yborders_new(j+1)=yborders_new(j)+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(j)*fsave_ratio)))
                   yborders_new(j+1)=min(yborders_new(j+1),yborders(nSubDomainsY)-max(fsave_ratio,4)*(nSubDomainsY-j-1)) !to set maximum limit on domain enlargement		
           end do
		   !xborders_new(nSubDomainsX-1)=xborders(nSubDomainsX-1)
           yborders_new(nSubDomainsY)=yborders(nSubDomainsY)
		   
		   
	       do i=0,nSubDomainsY
			    if(abs(yborders_new(i)-yborders(i)).gt.0) yborders_changed=.true.
		   end do  
		  
           if(yborders_changed) then
                if(proc.eq.0) then 
					print*,'Adjusting the Yborders to balance the Load ...'
				    !print*,'xborders',yborders
				    !print*,'xborders_new',yborders_new
					!print*,'sec_time',Max_TotalTime_column
				end if
                call SetNewYBorders
	            call ReorderPrtlArr
				call ReorderTestPrtlArr
           end if
		   

       end if

 end subroutine LoadBalanceY
 
 subroutine LoadBalanceYRecn
      integer :: i,j,k,dYdomain 
#ifdef twoD
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: exectime_grid
#else 
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: exectime_grid
#endif 	  

      real, dimension(0:nSubDomainsY-1) :: Max_TotalTime_column
      real :: time_temp
      logical :: yborders_changed
 
      yborders_changed=.false.
	  if(nSubDomainsY.le.3) return ! too few proc. for load balancing
	  

      if(modulo(t,domain_change_period).eq.0) then
           exectime_grid=0
           dYdomain=yborders(procyind(proc)+1)-yborders(procyind(proc))
		   exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/dYdomain
           call BcastExecTimeAllHomogeneous(exectime_grid)
         !first calculate the average execution time of all proc in a column (i=const.)
           Max_TotalTime_column=0
           do i=0,nSubDomainsX-1
               do j=0,nSubDomainsY-1
#ifdef twoD
                   do k=0,0
#else				   
				   do k=0,nSubDomainsZ-1
#endif
                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(j)) Max_TotalTime_column(j)=exectime_grid(i,j,k)
				   end do 
               end do
           end do
          
           !Adjust y boundaries in the lower half
           time_temp=0
           do j=0,nSubDomainsY/2-1
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(j)
           end do
           time_temp=(yborders(nSubDomainsY/2)-yborders(0))/time_temp

		   
           yborders_new(0)=yborders(0)
           do j=0,nSubDomainsY/2-2
                   yborders_new(j+1)=yborders_new(j)+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(j)*fsave_ratio)))				   
                   yborders_new(j+1)=min(yborders_new(j+1),yborders(nSubDomainsY)-max(fsave_ratio,4)*(nSubDomainsY-j-1)) !to set maximum limit on domain enlargement		
           end do
           yborders_new(nSubDomainsY/2)=yborders(nSubDomainsY/2)
		   !now adjust boundaries above the mid plane
           time_temp=0
           do j=nSubDomainsY/2,nSubDomainsY-1
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(j)
           end do
           time_temp=(yborders(nSubDomainsY)-yborders(nSubDomainsY/2))/time_temp
		   
           do j=nSubDomainsY/2,nSubDomainsY-2
                   yborders_new(j+1)=yborders_new(j)+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(j)*fsave_ratio)))
                   yborders_new(j+1)=min(yborders_new(j+1),yborders(nSubDomainsY)-max(fsave_ratio,4)*(nSubDomainsY-j-1)) !to set maximum limit on domain enlargement		
           end do
           yborders_new(nSubDomainsY)=yborders(nSubDomainsY)
		   
		   
		   
	       do i=0,nSubDomainsY
			    if(abs(yborders_new(i)-yborders(i)).gt.0) yborders_changed=.true.
		   end do  
		  
           if(yborders_changed) then
                if(proc.eq.0) then 
					print*,'Adjusting the Yborders to balance the Load ...'
				    !print*,'xborders',yborders
				    !print*,'xborders_new',yborders_new
				end if
                call SetNewYBorders
	            call ReorderPrtlArr
				call ReorderTestPrtlArr
           end if
		   

       end if

 end subroutine LoadBalanceYRecn

 subroutine LoadBalanceShock(xinj,xwall)
      integer :: i,j,k,xwall,dXdomain
	  real(dbpsn) :: xinj 
#ifdef twoD
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: exectime_grid
#else 
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: exectime_grid
#endif 	  

      real, dimension(0:nSubDomainsX-1) :: Max_TotalTime_column
      real :: time_temp
	  real(psn) :: mean_dx 
      integer :: loadbalance_profiling_ind
	  integer :: xinj_proc
      logical :: xborders_changed
      save loadbalance_profiling_ind
      data loadbalance_profiling_ind /1/
      xborders_changed=.false.
      !profile the performance history of all processors before making any decison about balancing the laod
!       if(modulo(t,loadbalance_profiling_period).eq.0) then
!            if(loadbalance_profiling_ind.gt.10) loadbalance_profiling_ind=1
!            hist_TotalTime(loadbalance_profiling_ind)=real(exec_time(31))
!            loadbalance_profiling_ind=loadbalance_profiling_ind+1
!            mean_TotalTime=average(hist_TotalTime,10)
!       end if

       if(modulo(t,domain_change_period).eq.0) then
           exectime_grid=0
           dXdomain=xborders(procxind(proc)+1)-max(xborders(procxind(proc)),xwall)
		   dXdomain=max(1,dXdomain)
		   exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/dXdomain
           call BcastExecTimeAllHomogeneous(exectime_grid)

         !first calculate the average execution time of all proc in a column (i=const.)
           Max_TotalTime_column=0
           do i=0,nSubDomainsX-1
               do j=0,nSubDomainsY-1
#ifdef twoD
                   do k=0,0
#else				   
				   do k=0,nSubDomainsZ-1
#endif
                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(i)) Max_TotalTime_column(i)=exectime_grid(i,j,k)
				   end do 
               end do
           end do
          
		   !Find the proc. index where injector lies 
		   do i=0,nSubDomainsX-1
			   if((xinj.ge.xborders(i)).and.(xinj.lt.xborders(i+1))) then 
				   xinj_proc=i 
			   end if 
		   end do 
           !Adjust the x-boundaries to redistribute load, and do not include last proc along x in load balancing 
           time_temp=0
           do i=0,nSubDomainsX-1
			    if(i.eq.xinj_proc) exit 
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(i)
           end do
           !time_temp=(xborders(nSubDomainsX-1)-xborders(0))/time_temp
           !time_temp=(xborders(xinj_proc)-xwall)/time_temp
           time_temp=(xinj-xwall)/time_temp
		   
		   

           xborders_new(0)=xborders(0)
           do i=0,nSubDomainsX-1
			   if(i+1.lt.xinj_proc) then 			 
                   !xborders_new(i+1)=xborders_new(i)+max(1,int(time_temp*(xborders(i+1)-xborders(i))/mean_TotalTime_column(i)))
                   xborders_new(i+1)=max(xwall,xborders_new(i))+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(i)*fsave_ratio)))
                   xborders_new(i+1)=min(xborders_new(i+1),xborders(xinj_proc)-max(fsave_ratio,4)*(xinj_proc-i-1)) !to set maximum limit on domain enlargement		
			   else
				   xborders_new(i+1)=xborders(i+1)
			   end if 	
           end do
		   !xborders_new(nSubDomainsX-1)=xborders(nSubDomainsX-1)
           !xborders_new(nSubDomainsX)=xborders(nSubDomainsX)
		   
		   !mean_dx=real(xborders(xinj_proc)-xwall)/real(3*xinj_proc) !Half of average domain size on each proc
		   !print*,'mean_dx',mean_dx,xborders(xinj_proc),fsave_ratio*floor(mean_dx/fsave_ratio)
		   !if(xinj.gt.xborders(xinj_proc)+mean_dx) then 
			   !if(xinj_proc.gt.0) 
			   xborders_new(xinj_proc)=xborders(xinj_proc)+fsave_ratio*floor(real(xinj-xborders(xinj_proc))/fsave_ratio)
			   !end if 
		   
		   
		   
! 		   if(proc.eq.0) then
! 		print*,'xborders',xborders
! 		print*,'xborders_new',xborders_new
! 		if(proc.eq.0) print*,'Time',Max_TotalTime_column,exec_time(32)
! 	end if
	
	       do i=0,nSubDomainsX
			    if(abs(xborders_new(i)-xborders(i)).gt.0) xborders_changed=.true.
		   end do  
		  
           if(xborders_changed) then
                if(proc.eq.0) then 
					print*,'Adjusting the Xborders to balance the Load ...'
				    !print*,'xborders',xborders
				    !print*,'xborders_new',xborders_new
				end if
                call SetNewXBorders
	            call ReorderPrtlArr
				call ReorderTestPrtlArr
           end if

       end if

 end subroutine LoadBalanceShock
 
 
 

 subroutine LoadBalanceShockTwoStream(xinj)
	  implicit none 
	  integer :: i,j,k
	  real(dbpsn) :: xinj 
#ifdef twoD
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: exectime_grid
#else 
      real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: exectime_grid
#endif 	  

      real, dimension(0:nSubDomainsX-1) :: Max_TotalTime_column
      real :: time_temp
	  real(psn) :: mean_dx 
      integer :: loadbalance_profiling_ind
	  integer :: xinj_proc_left,xinj_proc_right
      logical :: xborders_changed
      save loadbalance_profiling_ind
      data loadbalance_profiling_ind /1/
      xborders_changed=.false.
      !profile the performance history of all processors before making any decison about balancing the laod
!       if(modulo(t,loadbalance_profiling_period).eq.0) then
!            if(loadbalance_profiling_ind.gt.10) loadbalance_profiling_ind=1
!            hist_TotalTime(loadbalance_profiling_ind)=real(exec_time(31))
!            loadbalance_profiling_ind=loadbalance_profiling_ind+1
!            mean_TotalTime=average(hist_TotalTime,10)
!       end if

       if(modulo(t,domain_change_period).eq.0) then
           exectime_grid=0
		   exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/(xborders(procxind(proc)+1)-xborders(procxind(proc)))
           call BcastExecTimeAllHomogeneous(exectime_grid)

         !first calculate the average execution time of all proc in a column (i=const.)
           Max_TotalTime_column=0
           do i=0,nSubDomainsX-1
               do j=0,nSubDomainsY-1
#ifdef twoD
                   do k=0,0
#else				   
				   do k=0,nSubDomainsZ-1
#endif
                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(i)) Max_TotalTime_column(i)=exectime_grid(i,j,k)
				   end do 
               end do
           end do
          
		   !Find the proc. index where injector lies 
		   do i=0,nSubDomainsX-1
			   if((xinj.ge.xborders(i)).and.(xinj.lt.xborders(i+1))) then 
				   xinj_proc_right=i 
			   end if 
			   if((-xinj.gt.xborders(i)).and.(-xinj.le.xborders(i+1))) then 
				   xinj_proc_left=i 
			   end if 
		   end do 
		   
           !----Adjust the x-boundaries at Secons  Half of the box redistribute load, and do not include last proc along x in load balancing 
           time_temp=0
           do i=nSubDomainsX/2,nSubDomainsX-1
			    if(i.eq.xinj_proc_right) exit 
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(i)
           end do
           !time_temp=(xborders(nSubDomainsX-1)-xborders(0))/time_temp
           time_temp=(xborders(xinj_proc_right)-xborders(nSubDomainsX/2))/time_temp
		   
           xborders_new(nSubDomainsX/2)=xborders(nSubDomainsX/2)
           do i=nSubDomainsX/2,nSubDomainsX-1
			   if(i+1.lt.xinj_proc_right) then 			 
                   xborders_new(i+1)=xborders_new(i)+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(i)*fsave_ratio)))
                   xborders_new(i+1)=min(xborders_new(i+1),xborders(xinj_proc_right)-max(fsave_ratio,4)*(xinj_proc_right-i-1)) !to set maximum limit on domain enlargement		
			   else
				   xborders_new(i+1)=xborders(i+1)
			   end if 	
           end do

           mean_dx=real(2*(xinj_proc_right-nSubDomainsX/2))
		   if(mean_dx.ne.0) then !this is to ensure that mean_dx does not become NaN
			   mean_dx=real(xborders(xinj_proc_right)-xborders(nSubDomainsX/2))/mean_dx !Half of average domain size on each proc
			   if(xinj.gt.xborders(xinj_proc_right)+mean_dx) then 
				   if(xinj_proc_right.gt.nSubDomainsX/2) xborders_new(xinj_proc_right)=xborders(xinj_proc_right)+fsave_ratio*floor(real(xinj-xborders(xinj_proc_right))/fsave_ratio)
			   end if 
	       end if 
		   !-----Now adjust the domains for Left half -----   
           time_temp=0
           do i=nSubDomainsX/2-1,0,-1
			    if(i.eq.xinj_proc_left) exit 
                time_temp=time_temp+1.0_psn/Max_TotalTime_column(i)
           end do
           !time_temp=(xborders(nSubDomainsX-1)-xborders(0))/time_temp
           time_temp=(xborders(nSubDomainsX/2)-xborders(xinj_proc_left+1))/time_temp
		   

           do i=nSubDomainsX/2-1,0,-1
			   if(i-1.gt.xinj_proc_left) then 			 
                   xborders_new(i)=xborders_new(i+1)-max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(i)*fsave_ratio)))
                   xborders_new(i)=max(xborders_new(i),xborders(xinj_proc_left+1)+max(fsave_ratio,4)*(i-xinj_proc_left-1)) !to set maximum limit on domain enlargement		
			   else
				   xborders_new(i)=xborders(i)
			   end if 	
           end do
           mean_dx=real(2*(nSubDomainsX/2-xinj_proc_left-1)) 
		   if(mean_dx.ne.0) then ! this is to ensure thhat actual mean_dx does not become NaN
			   mean_dx=real(xborders(nSubDomainsX/2)-xborders(xinj_proc_left+1))/mean_dx !Half of average domain size on each proc
			   if(-xinj.lt.xborders(xinj_proc_left+1)-mean_dx) then 
				   if(xinj_proc_left.lt.nSubDomainsX/2-1) xborders_new(xinj_proc_left+1)=xborders(xinj_proc_left+1)-fsave_ratio*floor(real(xborders(xinj_proc_left+1)+xinj)/fsave_ratio)
			   end if 
	       end if 
		   
		  		   
   		   !if(proc.eq.0) print*,'Time',Max_TotalTime_column,exec_time(32)
	
	
	       do i=0,nSubDomainsX
			    if(abs(xborders_new(i)-xborders(i)).gt.0) xborders_changed=.true.
		   end do  
		  
           if(xborders_changed) then
                if(proc.eq.0) then 
					print*,'Adjusting the Xborders to balance the Load ...'
				    !print*,'xborders',xborders
				    !print*,'xborders_new',xborders_new
				end if
                call SetNewXBorders
	            call ReorderPrtlArr
	            call ReorderTestPrtlArr
           end if

       end if

 end subroutine LoadBalanceShockTwoStream
 
 
 
 
 real function average(arr,len)
      integer :: len,i
      real, dimension(len) :: arr
      real :: sum
      sum=0.0_psn
      do i=1,len
           sum=sum+arr(i)
      end do
      average=sum/len     
 end function average

subroutine DefineNewBorders_default
    !call BcastExecTime
     xborders_new=xborders !nothing is changed 
end subroutine DefineNewBorders_default




 subroutine SetNewXBorders
     integer :: i,procxind_this
     integer :: lmost_recv,rmost_recv,lmost_send,rmost_send
     integer, dimension(:), allocatable :: segment_borders_send,segment_borders_recv 

     procxind_this=procxind(proc)

     do i=0,nSubDomainsX-1 !determine left and rightmost proc. needed to communicate with 
          if((xborders_new(procxind_this).ge.xborders(i)).and.(xborders_new(procxind_this).lt.xborders(i+1))) lmost_recv=i
          if((xborders_new(procxind_this+1).le.xborders(i+1)).and.(xborders_new(procxind_this+1).gt.xborders(i))) rmost_recv=i
          if((xborders(procxind_this).ge.xborders_new(i)).and.(xborders(procxind_this).lt.xborders_new(i+1))) lmost_send=i
          if((xborders(procxind_this+1).le.xborders_new(i+1)).and.(xborders(procxind_this+1).gt.xborders_new(i))) rmost_send=i                    
     end do
	 
	 !determine x-index of grid points that need to be transferred 
     allocate(segment_borders_send(lmost_send:rmost_send+1)) 
     allocate(segment_borders_recv(lmost_recv:rmost_recv+1))
     call DetermineLocalXBorders(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)
          
     mx_new=xborders_new(procxind_this+1)-xborders_new(procxind_this)+5!new local domain size 
     my_new=my !this subroutine is designed to change domain only in x-direction
     mz_new=mz

     deallocate(F0)
     allocate(F0(mx_new,my_new,mz_new))

     !Send-Recieve EM field 
	 call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Ex,1)
     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Ey,2)
     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Ez,3)
     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Bx,4)
     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,By,5)
     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Bz,6)
	 if(nMoverEMfilter.gt.0) then 
		 call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,filteredEx,7)
	     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,filteredEy,8)
	     call SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,filteredEz,9)
	 end if
     !Now Send-Recieve particle data
     call SendRecvPrtlSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)
     call SendRecvTestPrtlSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)

     deallocate(segment_borders_send,segment_borders_recv)

     !update all other relevant auxiliary arrays and varaiables local to this domain 
     mx=mx_new
     my=my_new
     mz=mz_new     
     xmax=mx-2
     xlen=xmax-3

     deallocate(Jx,Jy,Jz)
     allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz))
     Jx=0.0_psn!reset the current 
     Jy=0.0_psn
     Jz=0.0_psn 
     deallocate(buff_tJx,buff_tJy,buff_tJz,buff_bJx,buff_bJy,buff_bJz)
     allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz)) 
#ifndef twoD	 
     deallocate(buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz)
	 allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))    
#endif

     call ReshapeShortMoverFldArr
	 	 
	 call ExchangeYZEdgeField(Ex,Ey,Ez)
     call ExchangeYZEdgeField(Bx,By,Bz) 
     call ExchangeZXEdgeField(Ex,Ey,Ez)
     call ExchangeZXEdgeField(Bx,By,Bz)
#ifndef twoD      
     call ExchangeXYEdgeField(Ex,Ey,Ez)
     call ExchangeXYEdgeField(Bx,By,Bz)
#endif 
if(nMoverEMfilter.gt.0) then 
	 call ExchangeYZEdgeField(filteredEx,filteredEy,filteredEz)
     call ExchangeZXEdgeField(filteredEx,filteredEy,filteredEz)
#ifndef twoD      
     call ExchangeXYEdgeField(filteredEx,filteredEy,filteredEz)
#endif 
end if


     xborders=xborders_new !now update the x-borders   
end subroutine SetNewXborders

subroutine DetermineLocalXBorders(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)
     integer :: i,procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv     
     integer, dimension(lmost_send:rmost_send+1)::segment_borders_send
     integer, dimension(lmost_recv:rmost_recv+1)::segment_borders_recv
     ! local boundaries of all segments in old data 
     do i=lmost_send,rmost_send
          segment_borders_send(i)=max(xborders(procxind_this),xborders_new(i))-xborders(procxind_this)+3
          segment_borders_send(i+1)=min(xborders(procxind_this+1),xborders_new(i+1))-xborders(procxind_this)+3 !can be optimised 
     end do
     !local borders of all segments in new data
     do i=lmost_recv,rmost_recv
          segment_borders_recv(i)=max(xborders_new(procxind_this),xborders(i))-xborders_new(procxind_this)+3
          segment_borders_recv(i+1)=min(xborders_new(procxind_this+1),xborders(i+1))-xborders_new(procxind_this)+3  
     end do 
end subroutine DetermineLocalXBorders

subroutine SendRecvFldSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv,Fld,tag)
     integer :: procxind_this,i,j,k,tag
     integer :: lmost_send,rmost_send,lmost_recv,rmost_recv
     integer, dimension(lmost_send:rmost_send+1)::segment_borders_send,borders_send
     integer, dimension(lmost_recv:rmost_recv+1)::segment_borders_recv
     real(psn),dimension(:,:,:), allocatable :: Fld

     j=procyind(proc)
	 k=proczind(proc)
      

    if((procxind_this.ge.lmost_send).and.(procxind_this.le.rmost_send))  then !copy the data that already exist on this proc
          F0(segment_borders_recv(procxind_this):segment_borders_recv(procxind_this+1)-1,:,:)=Fld(segment_borders_send(procxind_this):segment_borders_send(procxind_this+1)-1,:,:)
    end if

     do i=min(lmost_send,lmost_recv),max(rmost_send,rmost_recv)
          if(i.eq.procxind_this) cycle
          if((i.ge.lmost_send).and.(i.le.rmost_send)) then 
               call SendFldToNewProcX(proc_grid(i,j,k),segment_borders_send(i),segment_borders_send(i+1)-1,Fld,tag)
          else if((i.ge.lmost_recv).and.(i.le.rmost_recv)) then 
               call RecvFldFromOldProcX(proc_grid(i,j,k),segment_borders_recv(i),segment_borders_recv(i+1)-1,F0,tag,mx_new,my_new,mz_new) !Recieve data from other proc
          end if
     end do
     
	 deallocate(Fld)
     call move_alloc(F0,Fld)
     allocate(F0(mx_new,my_new,mz_new))     
end subroutine SendRecvFldSegmentsX     
!---------------------------------------------------------------------------------------------------------------------
! Send and recv all prtls from all relavant domains 
! Note that insertion of new particles into the main array is done by serially looking for emptly slots 
! It is different from normal insertion of new unsorted particles at the end of particle array 
! particle sorting is required after the load balance
!--------------------------------------------------------------------------------------------------------------------

subroutine SendRecvPrtlSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)
	 implicit none 
	 integer :: procxind_this,i,j,k,n,psend_size,pcount,incoming_np,precv_size
     integer :: lmost_send,rmost_send,lmost_recv,rmost_recv
     type(particle), dimension(:), allocatable :: psend,precv
     integer, dimension(lmost_send:rmost_send+1)::segment_borders_send,extent_send
     integer, dimension(lmost_recv:rmost_recv+1)::segment_borders_recv
     integer, dimension(lmost_send:rmost_send):: segment_np_send,prtl_start_ind
     integer, dimension(lmost_recv:rmost_recv):: segment_np_recv

     
     !scan through the particle data to get the size of the data to be sent to each proc
     extent_send=segment_borders_send ! an extended border is needed to include corner and edge particles
     extent_send(rmost_send+1)=extent_send(rmost_send+1)+1 
     extent_send(lmost_send)=extent_send(lmost_send)-1
     segment_np_send=0
     do n=1,used_prtl_arr_size
          if(qp(n).eq.0) cycle
          i=lmost_send
          do while(i.le.rmost_send) !this can be optimised, if the need be (if large number segments break away from one proc)
               if((xp(n).ge.extent_send(i)).and.(xp(n).lt.extent_send(i+1))) then
                    segment_np_send(i)=segment_np_send(i)+1
                    exit
               end if 
               i=i+1
          end do 
     end do 

     !find the offset value for each proc
     pcount=1
     prtl_start_ind=1
     psend_size=0
     do i=lmost_send,rmost_send
          if(i.eq.procxind_this) cycle !becuase particle data is not sent to self
          prtl_start_ind(i)=pcount          
          pcount=pcount+segment_np_send(i)
          psend_size=psend_size+segment_np_send(i)          
     end do
     
     allocate(psend(psend_size)) !all outliers would be loaded in psend
     np=np-psend_size !update number of particles in this proc

     !send the size of particle data to all relevant proc
     j=procyind(proc)
	 k=proczind(proc)

     segment_np_recv=0
     incoming_np=0

     do i=min(lmost_send,lmost_recv),max(rmost_send,rmost_recv)
          if(i.eq.procxind_this) cycle 
          if((i.ge.lmost_send).and.(i.le.rmost_send)) then 
               call SendPrtlSizeToNewProc(proc_grid(i,j,k),segment_np_send(i))
          else if((i.ge.lmost_recv).and.(i.le.rmost_recv)) then 
              call RecvPrtlSizeFromOldProc(proc_grid(i,j,k),segment_np_recv(i))
               incoming_np=incoming_np+segment_np_recv(i)     
          end if
     end do
     
     !now load outliers in all segments in psend array, ready to be sent in parts to it's new proc 
     segment_np_send=0 !reinitialise to recount particle, but this time load them in ptemp
     do n=1,used_prtl_arr_size
          if(qp(n).eq.0) cycle
         i=lmost_send
          do while(i.le.rmost_send) 
               if((xp(n).ge.extent_send(i)).and.(xp(n).lt.extent_send(i+1))) then                    
                    if(i.ne.procxind_this) then 
                       xp(n)=xp(n)-segment_borders_send(i)   
					   call LoadPrtl(psend,psend_size,prtl_start_ind(i)+segment_np_send(i),n)  
                       segment_np_send(i)=segment_np_send(i)+1
                       call DeletePrtl(n) !delete particle if the particle is leaving this proc                       
                    else !if particle remains on the same proc, make sure it's x-cordinate is updated
                       xp(n)=xp(n)+xborders(procxind_this)-xborders_new(procxind_this)                     
                    end if
                    exit
               end if 
               i=i+1
          end do 
     end do 
     
     
     
     !Now send the particles to their new proc
     np=np+incoming_np !update number of particles 
     if(np.gt.prtl_arr_size) call ReshapePrtlArr(int(1.05*np)) !make sure that there is enough place to accomodate all new particles

     precv_size=maxval(segment_np_recv)
     allocate(precv(precv_size))

     do i=min(lmost_send,lmost_recv),max(rmost_send,rmost_recv)
          if(i.eq.procxind_this) cycle 
          if((i.ge.lmost_send).and.(i.le.rmost_send)) then 
              call SendPrtlToNewProc(proc_grid(i,j,k),prtl_start_ind(i),prtl_start_ind(i)+segment_np_send(i)-1,segment_np_send(i),psend,psend_size)     
          else if((i.ge.lmost_recv).and.(i.le.rmost_recv)) then 
               call RecvPrtlFromOldProc(proc_grid(i,j,k),segment_np_recv(i),segment_np_recv(i),precv,precv_size,segment_borders_recv(i),0)
          end if
     end do
     
     deallocate(psend) !free the memory that was used to send-recv particles      
     deallocate(precv)
     used_prtl_arr_size=prtl_arr_size
end subroutine SendRecvPrtlSegmentsX


subroutine SendRecvTestPrtlSegmentsX(procxind_this,lmost_send,rmost_send,lmost_recv,rmost_recv,segment_borders_send,segment_borders_recv)
	 implicit none 
	 integer :: procxind_this,i,j,k,n,psend_size,pcount,incoming_np,precv_size
     integer :: lmost_send,rmost_send,lmost_recv,rmost_recv
     type(particle), dimension(:), allocatable :: psend,precv
     integer, dimension(lmost_send:rmost_send+1)::segment_borders_send,extent_send
     integer, dimension(lmost_recv:rmost_recv+1)::segment_borders_recv
     integer, dimension(lmost_send:rmost_send):: segment_np_send,prtl_start_ind
     integer, dimension(lmost_recv:rmost_recv):: segment_np_recv

     
     !scan through the particle data to get the size of the data to be sent to each proc
     extent_send=segment_borders_send ! an extended border is needed to include corner and edge particles
     extent_send(rmost_send+1)=extent_send(rmost_send+1)+1 
     extent_send(lmost_send)=extent_send(lmost_send)-1
     segment_np_send=0
     do n=1,used_test_prtl_arr_size
          if(qtp(n).eq.0) cycle
          i=lmost_send
          do while(i.le.rmost_send) !this can be optimised, if the need be (if large number segments break away from one proc)
               if((xtp(n).ge.extent_send(i)).and.(xtp(n).lt.extent_send(i+1))) then
                    segment_np_send(i)=segment_np_send(i)+1
                    exit
               end if 
               i=i+1
          end do 
     end do 

     !find the offset value for each proc
     pcount=1
     prtl_start_ind=1
     psend_size=0
     do i=lmost_send,rmost_send
          if(i.eq.procxind_this) cycle !becuase particle data is not sent to self
          prtl_start_ind(i)=pcount          
          pcount=pcount+segment_np_send(i)
          psend_size=psend_size+segment_np_send(i)          
     end do
     
     allocate(psend(psend_size)) !all outliers would be loaded in psend
     ntp=ntp-psend_size !update number of particles in this proc

     !send the size of particle data to all relevant proc
     j=procyind(proc)
	 k=proczind(proc)

     segment_np_recv=0
     incoming_np=0

     do i=min(lmost_send,lmost_recv),max(rmost_send,rmost_recv)
          if(i.eq.procxind_this) cycle 
          if((i.ge.lmost_send).and.(i.le.rmost_send)) then 
               call SendPrtlSizeToNewProc(proc_grid(i,j,k),segment_np_send(i))
          else if((i.ge.lmost_recv).and.(i.le.rmost_recv)) then 
              call RecvPrtlSizeFromOldProc(proc_grid(i,j,k),segment_np_recv(i))
               incoming_np=incoming_np+segment_np_recv(i)     
          end if
     end do

     
     !now load outliers in all segments in psend array, ready to be sent in parts to it's new proc 
     segment_np_send=0 !reinitialise to recount particle, but this time load them in ptemp
     do n=1,used_test_prtl_arr_size
          if(qtp(n).eq.0) cycle
         i=lmost_send
          do while(i.le.rmost_send) 
               if((xtp(n).ge.extent_send(i)).and.(xtp(n).lt.extent_send(i+1))) then                    
                    if(i.ne.procxind_this) then 
                       xtp(n)=xtp(n)-segment_borders_send(i)   
					   call LoadTestPrtl(psend,psend_size,prtl_start_ind(i)+segment_np_send(i),n)  
                       segment_np_send(i)=segment_np_send(i)+1
                       call DeleteTestPrtl(n) !delete particle if the particle is leaving this proc                       
                    else !if particle remains on the same proc, make sure it's x-cordinate is updated
                       xtp(n)=xtp(n)+xborders(procxind_this)-xborders_new(procxind_this)                     
                    end if
                    exit
               end if 
               i=i+1
          end do 
     end do 
     
     
     
     !Now send the particles to their new proc
     ntp=ntp+incoming_np !update number of particles 
     if(ntp.gt.test_prtl_arr_size) call ReshapeTestPrtlArr(int(1.05*ntp)) !make sure that there is enough place to accomodate all new particles

     precv_size=maxval(segment_np_recv)
     allocate(precv(precv_size))

     do i=min(lmost_send,lmost_recv),max(rmost_send,rmost_recv)
          if(i.eq.procxind_this) cycle 
          if((i.ge.lmost_send).and.(i.le.rmost_send)) then 
              call SendPrtlToNewProc(proc_grid(i,j,k),prtl_start_ind(i),prtl_start_ind(i)+segment_np_send(i)-1,segment_np_send(i),psend,psend_size)     
          else if((i.ge.lmost_recv).and.(i.le.rmost_recv)) then 
               call RecvTestPrtlFromOldProc(proc_grid(i,j,k),segment_np_recv(i),segment_np_recv(i),precv,precv_size,segment_borders_recv(i),0)
          end if
     end do
     
     deallocate(psend) !free the memory that was used to send-recv particles      
     deallocate(precv)
	 used_test_prtl_arr_size=test_prtl_arr_size	 
end subroutine SendRecvTestPrtlSegmentsX


!----------------------------------------------------------------------------------------------------------------------
! The following subroutines are used to shift yborders for load balancing
!----------------------------------------------------------------------------------------------------------------------

 subroutine SetNewYBorders
     integer :: i,procyind_this
     integer :: bmost_recv,tmost_recv,bmost_send,tmost_send
     integer, dimension(:), allocatable :: segment_borders_send,segment_borders_recv 

     procyind_this=procyind(proc)

     do i=0,nSubDomainsY-1 !determine left and rightmost proc. needed to communicate with 
          if((yborders_new(procyind_this).ge.yborders(i)).and.(yborders_new(procyind_this).lt.yborders(i+1))) bmost_recv=i
          if((yborders_new(procyind_this+1).le.yborders(i+1)).and.(yborders_new(procyind_this+1).gt.yborders(i))) tmost_recv=i
          if((yborders(procyind_this).ge.yborders_new(i)).and.(yborders(procyind_this).lt.yborders_new(i+1))) bmost_send=i
          if((yborders(procyind_this+1).le.yborders_new(i+1)).and.(yborders(procyind_this+1).gt.yborders_new(i))) tmost_send=i                    
     end do
	 
	 !determine y-index of grid points that need to be transferred 
     allocate(segment_borders_send(bmost_send:tmost_send+1)) 
     allocate(segment_borders_recv(bmost_recv:tmost_recv+1))
     call DetermineLocalYBorders(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)
          
     mx_new=mx
	 my_new=yborders_new(procyind_this+1)-yborders_new(procyind_this)+5!new local domain size 
     mz_new=mz
     deallocate(F0)
     allocate(F0(mx_new,my_new,mz_new))

     !Send-Recieve EM field 
	 call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Ex,1)
     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Ey,2)
     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Ez,3)
     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Bx,4)
     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,By,5)
     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Bz,6)
	 if(nMoverEMfilter.gt.0) then 
		 call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,filteredEx,7)
	     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,filteredEy,8)
	     call SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,filteredEz,9)
	 end if	 
     !Now Send-Recieve particle data
     call SendRecvPrtlSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)
     call SendRecvTestPrtlSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)

     deallocate(segment_borders_send,segment_borders_recv)

     !update all other relevant auxiliary arrays and varaiables local to this domain 
     mx=mx_new
     my=my_new
     mz=mz_new     
     ymax=my-2
     ylen=ymax-3

     deallocate(Jx,Jy,Jz)
     allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz))
     Jx=0.0_psn!reset the current 
     Jy=0.0_psn
     Jz=0.0_psn 
     deallocate(buff_rJx,buff_rJy,buff_rJz,buff_lJx,buff_lJy,buff_lJz)
     allocate(buff_lJx(3,my,mz),buff_rJx(3,my,mz),buff_lJy(3,my,mz),buff_rJy(3,my,mz),buff_lJz(3,my,mz),buff_rJz(3,my,mz)) 
#ifndef twoD	 
     deallocate(buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz)
	 allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))    
#endif

     call ReshapeShortMoverFldArr
	 	 
	 call ExchangeYZEdgeField(Ex,Ey,Ez)
     call ExchangeYZEdgeField(Bx,By,Bz) 
     call ExchangeZXEdgeField(Ex,Ey,Ez)
     call ExchangeZXEdgeField(Bx,By,Bz)
#ifndef twoD      
     call ExchangeXYEdgeField(Ex,Ey,Ez)
     call ExchangeXYEdgeField(Bx,By,Bz)
#endif 
if(nMoverEMfilter.gt.0) then 
	 call ExchangeYZEdgeField(filteredEx,filteredEy,filteredEz)
     call ExchangeZXEdgeField(filteredEx,filteredEy,filteredEz)
#ifndef twoD      
     call ExchangeXYEdgeField(filteredEx,filteredEy,filteredEz)
#endif 
end if


     yborders=yborders_new !now update the x-borders   
end subroutine SetNewYborders

subroutine DetermineLocalYBorders(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)
     integer :: i,procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv     
     integer, dimension(bmost_send:tmost_send+1)::segment_borders_send
     integer, dimension(bmost_recv:tmost_recv+1)::segment_borders_recv
     ! local boundaries of all segments in old data 
     do i=bmost_send,tmost_send
          segment_borders_send(i)=max(yborders(procyind_this),yborders_new(i))-yborders(procyind_this)+3
          segment_borders_send(i+1)=min(yborders(procyind_this+1),yborders_new(i+1))-yborders(procyind_this)+3 !can be optimised 
     end do
     !local borders of all segments in new data
     do i=bmost_recv,tmost_recv
          segment_borders_recv(i)=max(yborders_new(procyind_this),yborders(i))-yborders_new(procyind_this)+3
          segment_borders_recv(i+1)=min(yborders_new(procyind_this+1),yborders(i+1))-yborders_new(procyind_this)+3  
     end do 
end subroutine DetermineLocalYBorders
subroutine SendRecvFldSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv,Fld,tag)
     integer :: procyind_this,i,j,k,tag
     integer :: bmost_send,tmost_send,bmost_recv,tmost_recv
     integer, dimension(bmost_send:tmost_send+1)::segment_borders_send,borders_send
     integer, dimension(bmost_recv:tmost_recv+1)::segment_borders_recv
     real(psn),dimension(:,:,:), allocatable :: Fld

     i=procxind(proc)
	 k=proczind(proc)
      

    if((procyind_this.ge.bmost_send).and.(procyind_this.le.tmost_send))  then !copy the data that already exist on this proc
          F0(:,segment_borders_recv(procyind_this):segment_borders_recv(procyind_this+1)-1,:)=Fld(:,segment_borders_send(procyind_this):segment_borders_send(procyind_this+1)-1,:)
    end if

     do j=min(bmost_send,bmost_recv),max(tmost_send,tmost_recv)
          if(j.eq.procyind_this) cycle
          if((j.ge.bmost_send).and.(j.le.tmost_send)) then 
               call SendFldToNewProcY(proc_grid(i,j,k),segment_borders_send(j),segment_borders_send(j+1)-1,Fld,tag)
          else if((j.ge.bmost_recv).and.(j.le.tmost_recv)) then 
               call RecvFldFromOldProcY(proc_grid(i,j,k),segment_borders_recv(j),segment_borders_recv(j+1)-1,F0,tag,mx_new,my_new,mz_new) !Recieve data from other proc
          end if
     end do
     
	 deallocate(Fld)
     call move_alloc(F0,Fld)
     allocate(F0(mx_new,my_new,mz_new))     
end subroutine SendRecvFldSegmentsY  


subroutine SendRecvPrtlSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)
	 implicit none 
	 integer :: procyind_this,i,j,k,n,psend_size,pcount,incoming_np,precv_size
     integer :: bmost_send,tmost_send,bmost_recv,tmost_recv
     type(particle), dimension(:), allocatable :: psend,precv
     integer, dimension(bmost_send:tmost_send+1)::segment_borders_send,extent_send
     integer, dimension(bmost_recv:tmost_recv+1)::segment_borders_recv
     integer, dimension(bmost_send:tmost_send):: segment_np_send,prtl_start_ind
     integer, dimension(bmost_recv:tmost_recv):: segment_np_recv

     
     !scan through the particle data to get the size of the data to be sent to each proc
     extent_send=segment_borders_send ! an extended border is needed to include corner and edge particles
     extent_send(tmost_send+1)=extent_send(tmost_send+1)+1 
     extent_send(bmost_send)=extent_send(bmost_send)-1
     segment_np_send=0
     do n=1,used_prtl_arr_size
          if(qp(n).eq.0) cycle
          i=bmost_send
          do while(i.le.tmost_send) !this can be optimised, if the need be (if large number segments break away from one proc)
               if((yp(n).ge.extent_send(i)).and.(yp(n).lt.extent_send(i+1))) then
                    segment_np_send(i)=segment_np_send(i)+1
                    exit
               end if 
               i=i+1
          end do 
     end do 

     !find the offset value for each proc
     pcount=1
     prtl_start_ind=1
     psend_size=0
     do i=bmost_send,tmost_send
          if(i.eq.procyind_this) cycle !becuase particle data is not sent to self
          prtl_start_ind(i)=pcount          
          pcount=pcount+segment_np_send(i)
          psend_size=psend_size+segment_np_send(i)          
     end do
     
     allocate(psend(psend_size)) !all outliers would be loaded in psend
     np=np-psend_size !update number of particles in this proc

     !send the size of particle data to all relevant proc
     i=procxind(proc)
	 k=proczind(proc)

     segment_np_recv=0
     incoming_np=0

     do j=min(bmost_send,bmost_recv),max(tmost_send,tmost_recv)
          if(j.eq.procyind_this) cycle 
          if((j.ge.bmost_send).and.(j.le.tmost_send)) then 
               call SendPrtlSizeToNewProc(proc_grid(i,j,k),segment_np_send(j))
          else if((j.ge.bmost_recv).and.(j.le.tmost_recv)) then 
              call RecvPrtlSizeFromOldProc(proc_grid(i,j,k),segment_np_recv(j))
               incoming_np=incoming_np+segment_np_recv(j)     
          end if
     end do
     
     !now load outliers in all segments in psend array, ready to be sent in parts to it's new proc 
     segment_np_send=0 !reinitialise to recount particle, but this time load them in ptemp
     do n=1,used_prtl_arr_size
          if(qp(n).eq.0) cycle
          j=bmost_send
          do while(j.le.tmost_send) 
               if((yp(n).ge.extent_send(j)).and.(yp(n).lt.extent_send(j+1))) then                    
                    if(j.ne.procyind_this) then 
                       yp(n)=yp(n)-segment_borders_send(j)   
					   call LoadPrtl(psend,psend_size,prtl_start_ind(j)+segment_np_send(j),n)  
                       segment_np_send(j)=segment_np_send(j)+1
                       call DeletePrtl(n) !delete particle if the particle is leaving this proc                       
                    else !if particle remains on the same proc, make sure it's x-cordinate is updated
                       yp(n)=yp(n)+yborders(procyind_this)-yborders_new(procyind_this)                     
                    end if
                    exit
               end if 
               j=j+1
          end do 
     end do 
     
     
     
     !Now send the particles to their new proc
     np=np+incoming_np !update number of particles 
     if(np.gt.prtl_arr_size) call ReshapePrtlArr(int(1.05*np)) !make sure that there is enough place to accomodate all new particles

     precv_size=maxval(segment_np_recv)
     allocate(precv(precv_size))

     do j=min(bmost_send,bmost_recv),max(tmost_send,tmost_recv)
          if(j.eq.procyind_this) cycle 
          if((j.ge.bmost_send).and.(j.le.tmost_send)) then 
              call SendPrtlToNewProc(proc_grid(i,j,k),prtl_start_ind(j),prtl_start_ind(j)+segment_np_send(j)-1,segment_np_send(j),psend,psend_size)     
          else if((j.ge.bmost_recv).and.(j.le.tmost_recv)) then 
               call RecvPrtlFromOldProc(proc_grid(i,j,k),segment_np_recv(j),segment_np_recv(j),precv,precv_size,0,segment_borders_recv(j))
          end if
     end do
     
     deallocate(psend) !free the memory that was used to send-recv particles      
     deallocate(precv)
     used_prtl_arr_size=prtl_arr_size
end subroutine SendRecvPrtlSegmentsY


subroutine SendRecvTestPrtlSegmentsY(procyind_this,bmost_send,tmost_send,bmost_recv,tmost_recv,segment_borders_send,segment_borders_recv)
	 implicit none 
	 integer :: procyind_this,i,j,k,n,psend_size,pcount,incoming_np,precv_size
     integer :: bmost_send,tmost_send,bmost_recv,tmost_recv
     type(particle), dimension(:), allocatable :: psend,precv
     integer, dimension(bmost_send:tmost_send+1)::segment_borders_send,extent_send
     integer, dimension(bmost_recv:tmost_recv+1)::segment_borders_recv
     integer, dimension(bmost_send:tmost_send):: segment_np_send,prtl_start_ind
     integer, dimension(bmost_recv:tmost_recv):: segment_np_recv

     
     !scan through the particle data to get the size of the data to be sent to each proc
     extent_send=segment_borders_send ! an extended border is needed to include corner and edge particles
     extent_send(tmost_send+1)=extent_send(tmost_send+1)+1 
     extent_send(bmost_send)=extent_send(bmost_send)-1
     segment_np_send=0
     do n=1,used_test_prtl_arr_size
          if(qtp(n).eq.0) cycle
          i=bmost_send
          do while(i.le.tmost_send) !this can be optimised, if the need be (if large number segments break away from one proc)
               if((xtp(n).ge.extent_send(i)).and.(xtp(n).lt.extent_send(i+1))) then
                    segment_np_send(i)=segment_np_send(i)+1
                    exit
               end if 
               i=i+1
          end do 
     end do 

     !find the offset value for each proc
     pcount=1
     prtl_start_ind=1
     psend_size=0
     do i=bmost_send,tmost_send
          if(i.eq.procyind_this) cycle !becuase particle data is not sent to self
          prtl_start_ind(i)=pcount          
          pcount=pcount+segment_np_send(i)
          psend_size=psend_size+segment_np_send(i)          
     end do
     
     allocate(psend(psend_size)) !all outliers would be loaded in psend
     ntp=ntp-psend_size !update number of particles in this proc

     !send the size of particle data to all relevant proc
     i=procxind(proc)
	 k=proczind(proc)

     segment_np_recv=0
     incoming_np=0

     do j=min(bmost_send,bmost_recv),max(tmost_send,tmost_recv)
          if(j.eq.procyind_this) cycle 
          if((j.ge.bmost_send).and.(j.le.tmost_send)) then 
               call SendPrtlSizeToNewProc(proc_grid(i,j,k),segment_np_send(j))
          else if((j.ge.bmost_recv).and.(j.le.tmost_recv)) then 
              call RecvPrtlSizeFromOldProc(proc_grid(i,j,k),segment_np_recv(j))
               incoming_np=incoming_np+segment_np_recv(j)     
          end if
     end do

     
     !now load outliers in all segments in psend array, ready to be sent in parts to it's new proc 
     segment_np_send=0 !reinitialise to recount particle, but this time load them in ptemp
     do n=1,used_test_prtl_arr_size
          if(qtp(n).eq.0) cycle
          j=bmost_send
          do while(j.le.tmost_send) 
               if((xtp(n).ge.extent_send(j)).and.(xtp(n).lt.extent_send(j+1))) then                    
                    if(j.ne.procyind_this) then 
                       ytp(n)=ytp(n)-segment_borders_send(i)   
					   call LoadTestPrtl(psend,psend_size,prtl_start_ind(j)+segment_np_send(j),n)  
                       segment_np_send(j)=segment_np_send(j)+1
                       call DeleteTestPrtl(n) !delete particle if the particle is leaving this proc                       
                    else !if particle remains on the same proc, make sure it's x-cordinate is updated
                       xtp(n)=xtp(n)+yborders(procyind_this)-yborders_new(procyind_this)                     
                    end if
                    exit
               end if 
               j=j+1
          end do 
     end do 
     
     
     
     !Now send the particles to their new proc
     ntp=ntp+incoming_np !update number of particles 
     if(ntp.gt.test_prtl_arr_size) call ReshapeTestPrtlArr(int(1.05*ntp)) !make sure that there is enough place to accomodate all new particles

     precv_size=maxval(segment_np_recv)
     allocate(precv(precv_size))

     do j=min(bmost_send,bmost_recv),max(tmost_send,tmost_recv)
          if(j.eq.procyind_this) cycle 
          if((j.ge.bmost_send).and.(j.le.tmost_send)) then 
              call SendPrtlToNewProc(proc_grid(i,j,k),prtl_start_ind(j),prtl_start_ind(j)+segment_np_send(j)-1,segment_np_send(j),psend,psend_size)     
          else if((j.ge.bmost_recv).and.(j.le.tmost_recv)) then 
               call RecvTestPrtlFromOldProc(proc_grid(i,j,k),segment_np_recv(j),segment_np_recv(j),precv,precv_size,0,segment_borders_recv(j))
          end if
     end do
     
     deallocate(psend) !free the memory that was used to send-recv particles      
     deallocate(precv)
	 used_test_prtl_arr_size=test_prtl_arr_size	 
end subroutine SendRecvTestPrtlSegmentsY




!NOTE ::: The following subroutines are outdated. They were not update for structure of arrays. Currently not in use 
!------------------------------------------------------------------------------------------------------------------------
!subroutine to swap physical regions between two proc proc1 and proc2
!------------------------------------------------------------------------------------------------------------------------
! subroutine SwapDomain(proc1,proc2)
!      integer :: proc1,proc2
!      if(proc1.eq.proc2) return
!      if(proc.eq.proc1) call SwapSendRecvDomainData(proc2)
!      if(proc.eq.proc2) call SwapSendRecvDomainData(proc1)
!      call SwapUpdateNeighbourList(proc1,proc2)
! end subroutine SwapDomain
!
! subroutine SwapSendRecvDomainData(SwapProc)
!      integer :: SwapProc,new_prtl_arr_size,new_np,i
!      type(particle), dimension(:), allocatable :: SwapPrtl
!      !First Exchange fields and update all fld and related variables
!      call Swap_SendRecvFldSize(SwapProc,mx_new,my_new,mz_new)
!      !print *,'Before swap at',proc,'Ex(1,1,1) is',Ex(1,1,1)
!      call Swap_SendRecvFld(SwapProc,Ex,mx_new,my_new,mz_new)
!      !print *,'After swap at',proc,'Ex(1,1,1) is',Ex(1,1,1)
!
!      call Swap_SendRecvFld(SwapProc,Ey,mx_new,my_new,mz_new)
!      call Swap_SendRecvFld(SwapProc,Ez,mx_new,my_new,mz_new)
!      call Swap_SendRecvFld(SwapProc,Bx,mx_new,my_new,mz_new)
!      call Swap_SendRecvFld(SwapProc,By,mx_new,my_new,mz_new)
!      call Swap_SendRecvFld(SwapProc,Bz,mx_new,my_new,mz_new)
!      deallocate(Jx,Jy,Jz,F0)
!      !domain decomposition is not done in z direction
!      mx=mx_new
!      my=my_new
!      mz=mz_new
!      xmax=mx-2
!      ymax=my-2
!      xlen=xmax-3
!      ylen=ymax-3
!
!      allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz),F0(mx,my,mz))
!      Jx=0
!      Jy=0
!      Jz=0
!      deallocate(buff_tJx,buff_tJy,buff_tJz,buff_bJx,buff_bJy,buff_bJz)
!      deallocate(buff_lJx,buff_lJy,buff_lJz,buff_rJx,buff_rJy,buff_rJz)
!
!      allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz))
!      allocate(buff_lJx(3,my,mz),buff_rJx(3,my,mz),buff_lJy(3,my,mz),buff_rJy(3,my,mz),buff_lJz(3,my,mz),buff_rJz(3,my,mz))
!
!      !call Swap_SendRecvProcXYind(SwapProc) ! call this subroutine to exchage global location of the physical domain when savedata is modified to use only global x and y range instead of procx(y)ind
!
!      !Now exchange particles
!      !print *,'before swap at',proc,'p(1)x is',p(1)%x
!      call Swap_SendRecvPrtlSize(SwapProc,new_prtl_arr_size,new_np)
!      allocate(SwapPrtl(new_prtl_arr_size))
!      call Swap_SendRecvPrtl(SwapProc,SwapPrtl,new_prtl_arr_size)
!      call move_alloc(SwapPrtl,p)
!      !print *,'after swap at',proc,'p(1)x is',p(1)%x
!      np=new_np
!      prtl_arr_size=new_prtl_arr_size
!      !update the rank of neighbouring processors
!      call Swap_SendRecvNeighbourProc(SwapProc,lproc)
!      call Swap_SendRecvNeighbourProc(SwapProc,rproc)
!      call Swap_SendRecvNeighbourProc(SwapProc,tproc)
!      call Swap_SendRecvNeighbourProc(SwapProc,bproc)
!
!      !update the list of free slots
!      deallocate(pfree)
!      allocate(pfree(prtl_arr_size))
!      pfree=0
!      pfree_ind=0
!      do i=1,prtl_arr_size
!           if(p(i)%q.eq.0) then
!                pfree_ind=pfree_ind+1
!                pfree(pfree_ind)=i
!           end if
!      end do
!
! end subroutine SwapSendRecvDomainData
! subroutine SwapUpdateNeighbourList(proc1,proc2)
!      integer :: proc1,proc2
!      call SwapUpdateNeighbourList1(lproc,proc1,proc2)
!      call SwapUpdateNeighbourList1(rproc,proc1,proc2)
!      call SwapUpdateNeighbourList1(tproc,proc1,proc2)
!      call SwapUpdateNeighbourList1(bproc,proc1,proc2)
! end subroutine SwapUpdateNeighbourList
! subroutine SwapUpdateNeighbourList1(ngb_list,proc1,proc2)
!      integer :: proc1,proc2
!      integer :: ngb_list
!      if(ngb_list.eq.proc1) then
!            ngb_list=proc2
!      else if(ngb_list.eq.proc2) then
!            ngb_list=proc1
!      end if
! end subroutine SwapUpdateNeighbourList1

!--------------------- End of swapping proc. subroutines ------------------------------------------------------------------

end module loadbalance
