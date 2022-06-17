module cyl_loadbalance
    use parameters
    use vars
	use cyl_vars
	use cyl_comm_fldprtl
	use loadbalance
	use movdep
contains 
    subroutine LoadBalanceR
         integer :: i,j,k,dXdomain,imin, n
         real, dimension(0:nSubDomainsX-1) :: Max_TotalTime_column
         real :: time_temp, f

          if(modulo(t,load_balance_period).ne.0) return
		  
!           exectime_grid=0
!           dXdomain=xborders(procxind(proc)+1)-xborders(procxind(proc))
!    		  exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/dXdomain
!           call BcastExecTimeAll(exectime_grid)
!
!           Max_TotalTime_column=0
!           do i=0,nSubDomainsX-1
!               do j=0,nSubDomainsY-1
! #ifdef twoD
!                   do k=0,0
! #else
! 			   do k=0,nSubDomainsZ-1
! #endif
!                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(i)) Max_TotalTime_column(i)=exectime_grid(i,j,k)
! 			   end do
!               end do
!           end do

          call CalcNumPrtlX	  
		  call ReduceNumPrtlX(npx,-1,nx,nx+2)
  
		  if(inc_axis) npx(-1:3)=0
		  call NumPrtlX_CDF
          
          call DetermineNewBordersR

          call CommitNewBordersX
    end subroutine LoadBalanceR
	
	subroutine BalanceDomainSizeR(Den)
		integer :: n
		real(psn), external :: Den 
		do n=1,ny ! hsould be nx
		   npx(n)=n*Den(n-1.0_psn,0.0_psn,0.0_psn)
	    end do
		 
	    if(inc_axis) npx(1:4)=0
	    call NumPrtlX_CDF
		
        call DetermineNewBordersR
        
		mx=xborders_new(procxind(proc)+1)-xborders_new(procxind(proc))+5
        xmax=mx-2
        xlen=xmax-3
        deallocate(Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz,F0)
		allocate(Ex(mx,my,mz),Ey(mx,my,mz),Ez(mx,my,mz))
        allocate(Bx(mx,my,mz),By(mx,my,mz),Bz(mx,my,mz))
        allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz))
        allocate(F0(mx,my,mz))
		Ex=0.0_psn;Ey=0.0_psn;Ez=0.0_psn;Bx=0.0_psn;By=0.0_psn;Bz=0.0_psn
		Jx=0.0_psn;Jy=0.0_psn;Jz=0.0_psn;F0=0.0_psn;
        if(nMoverEMfilter.gt.0) then 
		   deallocate(FilteredEx,FilteredEy,FilteredEz)
		   allocate(FilteredEx(mx,my,mz),FilteredEy(mx,my,mz),FilteredEz(mx,my,mz))
		   FilteredEx=0;FilteredEy=0;FilteredEz=0;
        end if 
        deallocate(buff_tJx,buff_tJy,buff_tJz,buff_bJx,buff_bJy,buff_bJz)
        allocate(buff_bJx(mx,3,mz),buff_tJx(mx,3,mz),buff_bJy(mx,3,mz),buff_tJy(mx,3,mz),buff_bJz(mx,3,mz),buff_tJz(mx,3,mz)) 
#ifndef twoD	 
        deallocate(buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz)
     	allocate(buff_dJx(mx,my,3),buff_uJx(mx,my,3),buff_dJy(mx,my,3),buff_uJy(mx,my,3),buff_dJz(mx,my,3),buff_uJz(mx,my,3))    
#endif
		
		call ReshapeShortMoverFldArr
		xborders=xborders_new
	end subroutine BalanceDomainSizeR
	
	
	!Adjust the x-boundaries to redistribute load; do not include leftmost proc if the axis is in the domain 
	subroutine DetermineNewBordersR
	    integer :: imin, n, i
		real :: f
		imin=0
	    if(inc_axis) imin=1

        xborders_new(0)=xborders(0)
	    if(inc_axis) xborders_new(1)=xborders(1)
        do i=imin,nSubDomainsX-2
		  f=real(i+1)/nSubDomainsX
		  if(inc_axis) f=real(i)/(nSubDomainsX-1)
		  call ValInd(npx,f,n)
            xborders_new(i+1)=max(xborders_new(i)+4,fsave_ratio*ceiling(real(n)/fsave_ratio))
            xborders_new(i+1)=min(xborders_new(i+1),xborders(nSubDomainsX)-max(fsave_ratio,4)*(nSubDomainsX-i-1)) !to set maximum limit on domain enlargement		
        end do
        xborders_new(nSubDomainsX)=xborders(nSubDomainsX)
	    rborders=xborders_new-0.5_psn
	end subroutine DetermineNewBordersR
		
	
	subroutine ValInd(Arr,val,index)
	     real(dbpsn),intent(in), dimension(-1:nx):: Arr
	     real ,intent(in) :: val
	     integer :: index,imin,imax,imid,imid_prev
	     logical :: search
     
	     imin=-1
	     imax=nx
		 imid=(imin+imax)/2
	     search=.true. 
		 do while(search)
			  imid_prev=imid 
	          if(Arr(imid).ge.val) imax=imid
	          if(Arr(imid).le.val) imin=imid
			  imid=(imin+imax)/2
	          if(imin+1.ge.imax) search=.false.
			  if(imid.eq.imid_prev) search=.false. !to ensure exit from the loop if Arr(imin)=Arr(imax) and imin+1<imax 		  
	     end do 
	     index=imin	 
	end subroutine ValInd 
	
!     subroutine LoadBalanceR
!          integer :: i,j,k,dXdomain,imin
!          real, dimension(0:nSubDomainsX-1) :: Max_TotalTime_column
!          real :: time_temp
!
!           if(modulo(t,load_balance_period).ne.0) return
!           exectime_grid=0
!           dXdomain=xborders(procxind(proc)+1)-xborders(procxind(proc))
!    		  exectime_grid(procxind(proc),procyind(proc),proczind(proc))=real(max(np,10))/dXdomain
!           call BcastExecTimeAll(exectime_grid)
!
!           Max_TotalTime_column=0
!           do i=0,nSubDomainsX-1
!               do j=0,nSubDomainsY-1
! #ifdef twoD
!                   do k=0,0
! #else
! 			   do k=0,nSubDomainsZ-1
! #endif
!                       if(exectime_grid(i,j,k).gt.Max_TotalTime_column(i)) Max_TotalTime_column(i)=exectime_grid(i,j,k)
! 			   end do
!               end do
!           end do
!
!           !Adjust the x-boundaries to redistribute load; do not include leftmost proc if the axis is in the domain
!           time_temp=0
! 		  imin=0
! 		  if(inc_axis) imin=1
!           do i=imin,nSubDomainsX-1
!                time_temp=time_temp+1.0_psn/Max_TotalTime_column(i)
!           end do
! 		  time_temp=(xborders(nSubDomainsX)-xborders(imin))/time_temp
!
!
!           xborders_new(0)=xborders(0)
! 		  if(inc_axis) xborders_new(1)=xborders(1)
!           do i=imin,nSubDomainsX-2
!               xborders_new(i+1)=xborders_new(i)+max(4,fsave_ratio*ceiling(time_temp/(Max_TotalTime_column(i)*fsave_ratio)))
!               xborders_new(i+1)=min(xborders_new(i+1),xborders(nSubDomainsX)-max(fsave_ratio,4)*(nSubDomainsX-i-1)) !to set maximum limit on domain enlargement
!           end do
!           xborders_new(nSubDomainsX)=xborders(nSubDomainsX)
! 		  rborders=xborders_new-0.5_psn
!           call CommitNewBordersX
!     end subroutine LoadBalanceR
	
end module cyl_loadbalance