module cyl_comm_fldprtl
	use parameters
	use vars
	use cyl_vars
	use communication
	use mem_prtl
	implicit none
contains
	
	subroutine ExchangePrtlAxis		
		if(procxind.eq.0) call SortPrtlOut_Axis
		call SendRecvPrtlSize_Axis
		call SendRecvPrtl_Axis
		if(procxind.eq.0) call MergePrtlIn_Axis 
	end subroutine ExchangePrtlAxis
	
	subroutine SendRecvPrtlSize_Axis
        integer :: stat(MPI_STATUS_SIZE)     
        integer :: i, mpi_err, max_size  
		
		if(inc_axis.eqv..false.) return
	    if(nSubDomainsX.lt.2) return
		if(procxind.eq.0) then 
			do i=0,nSubDomainsY-1				
		       call MPI_RECV(rinp_count_axis(i),1,MPI_INTEGER,proc_grid(1,i,proczind),0,MPI_COMM_WORLD,stat,mpi_err)
		       call MPI_RECV(rintp_count_axis(i),1,MPI_INTEGER,proc_grid(1,i,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
		    end do 
			do i=0,nSubDomainsY-1		
		       call MPI_SEND(rpcross_axis(i),1,MPI_INTEGER,proc_grid(1,i,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
		       call MPI_SEND(rtpcross_axis(i),1,MPI_INTEGER,proc_grid(1,i,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		    end do 
			!change size of the transfer arrays if needed
			max_size=maxval(rinp_count_axis)+maxval(rintp_count_axis)
			if(max_size.gt.rinp_size) then 
				deallocate(rinp_axis, rinp)
				max_size=int(1.1*max_size)+100
				allocate(rinp_axis(max_size,0:nSubDomainsY-1), rinp(nSubDomainsY*max_size))
				rinp_size=max_size
			end if
		end if
		
		if(procxind.eq.1) then			
			 call MPI_SEND(lpcross,1,MPI_INTEGER,proc_grid(0,0,proczind),0,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_SEND(ltpcross,1,MPI_INTEGER,proc_grid(0,0,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(linp_count,1,MPI_INTEGER,proc_grid(0,0,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(lintp_count,1,MPI_INTEGER,proc_grid(0,0,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
			 call UpdateTransferInSize
		end if	
			
	end subroutine SendRecvPrtlSize_Axis 
	
	subroutine SendRecvPrtl_Axis
        integer :: stat(MPI_STATUS_SIZE)     
        integer :: i, mpi_err  
		
	    if(nSubDomainsX.lt.2) return
		if(procxind.eq.0) then 
			do i=0,nSubDomainsY-1
		       call MPI_RECV(rinp_axis(1:rinp_count_axis(i)+rintp_count_axis(i),i),rinp_count_axis(i)+rintp_count_axis(i),mpi_prtltype,proc_grid(1,i,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
		    end do 
			do i=0,nSubDomainsY-1
		       call MPI_SEND(routp_axis(1:rpcross_axis(i)+rtpcross_axis(i),i),rpcross_axis(i)+rtpcross_axis(i),mpi_prtltype,proc_grid(1,i,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		    end do 
			call UpdateTransferOutSizeAxis
		end if
		
		if(procxind.eq.1) then
			 call MPI_SEND(loutp,lcross,mpi_prtltype,proc_grid(0,0,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(linp(1:linp_count+lintp_count),linp_count+lintp_count,mpi_prtltype,proc_grid(0,0,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		     call UpdateTransferOutSize
		end if
	end subroutine SendRecvPrtl_Axis
	
	subroutine UpdateTransferOutSizeAxis
		integer :: max_size
		max_size=maxval(rpcross_axis)+maxval(rtpcross_axis)
		if(routp_size.lt.(2*max_size+100000)) then
			max_size=int(2.5*max_size+100000)
			deallocate(routp_axis,routp)
			allocate(routp_axis(max_size,0:nSubDomainsY-1), routp(nSubDomainsY*max_size))
		end if  
	end subroutine UpdateTransferOutSizeAxis
	
	subroutine MergePrtlIn_Axis
		integer :: n, off
				
		do n=0,nSubDomainsY-1
			rinp_axis(1:rinp_count_axis(n)+rintp_count_axis(n),n)%y=rinp_axis(1:rinp_count_axis(n)+rintp_count_axis(n),n)%y+yborders(n)
		end do 
					
		off=0
		do n=0,nSubDomainsY-1 !particles
			rinp(1+off:rinp_count_axis(n)+off)=rinp_axis(1:rinp_count_axis(n),n)
			off=off+rinp_count_axis(n)
		end do 
		rinp_count=off
		do n=0,nSubDomainsY-1 !test particles
			rinp(1+off+rinp_count_axis(n):rinp_count_axis(n)+rintp_count_axis(n)+off)=rinp_axis(1+rinp_count_axis(n):rinp_count_axis(n)+rintp_count_axis(n),n)
			off=off+rintp_count_axis(n)
		end do 
		rintp_count=off-rinp_count
	end subroutine MergePrtlIn_Axis
	
	subroutine SortPrtlOut_Axis
		integer :: n, ind
		
		rpcross_axis=0
		rtpcross_axis=0
		do n=1,rpcross 
			ind=find_right_subdomain_ind(routp(n)%y-3)
			rpcross_axis(ind)=rpcross_axis(ind)+1
			routp(n)%y=routp(n)%y-yborders(ind)
			routp_axis(rpcross_axis(ind),ind)=routp(n)
		end do 
		do n=rpcross+1,rcross
			ind=find_right_subdomain_ind(routp(n)%y-3)
			rtpcross_axis(ind)=rtpcross_axis(ind)+1
			routp(n)%y=routp(n)%y-yborders(ind)
			routp_axis(rpcross_axis(ind)+rtpcross_axis(ind),ind)=routp(n)
		end do
		
	end subroutine SortPrtlOut_Axis
	
	integer function find_right_subdomain_ind(y)
		real(psn) :: y
		integer   :: imin,imax,imid
		logical   :: search
		search=.true.
		imin=0
		imax=nSubDomainsY-1
		imid=(imin+imax)/2
		if(y.le.yborders(1)) then 
			imid=0
			search=.false.
		end if 
		if(y.ge.yborders(nSubDomainsY-1)) then 
			imid=nSubDomainsY-1
			search=.false.
		end if 	
		
		do while(search)
           if(y.ge.yborders(imid)) imin=imid
           if(y.le.yborders(imid)) imax=imid
		   imid=(imin+imax)/2
		   if(imin+1.ge.imax) search=.false.
		end do 
		find_right_subdomain_ind=imid
	end function find_right_subdomain_ind
	
    
    subroutine ExchangeYZEdgeField_Axis(Fldx,Fldy,Fldz)
         real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz 
         integer :: i,off,my_short,dcount1,dcount2,mpi_err
         integer :: stat(MPI_STATUS_SIZE)
         dcount1=3*my*mz
         dcount2=2*my*mz
 		 if(inc_axis.eqv..false.) return
	     if(nSubDomainsX.lt.2) return
 		 if(procxind.eq.0) then
             do i=0,nSubDomainsY-1
				 off=yborders(i)+2
				 if(i.eq.0) off=0
				 my_short=yborders(i+1)-yborders(i) 
				 if(i.eq.0) my_short=my_short+2 
				 if(i.eq.nSubDomainsY-1) my_short=my_short+3
		         dcount1=3*my_short*mz
		         dcount2=2*my_short*mz			 
		         call MPI_RECV(Fldx(mx-2:mx,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind),1,MPI_COMM_WORLD,stat,mpi_err)				 
				 call MPI_RECV(Fldy(mx-2:mx,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
				 call MPI_RECV(Fldz(mx-2:mx,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Fldx(mx-4:mx-3,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(1,i,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Fldy(mx-4:mx-3,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(1,i,proczind),5,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Fldz(mx-4:mx-3,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(1,i,proczind),6,MPI_COMM_WORLD,stat,mpi_err)
			 end do 
		 end if 
		 if(procxind.eq.1) then
	         off=2 
			 if(procyind.eq.0) off=0
			 my_short=my-5 
			 if(procyind.eq.0) my_short=my_short+2
			 if(procyind.eq.nSubDomainsY-1) my_short=my_short+3
			 dcount1=3*my_short*mz
	         dcount2=2*my_short*mz		 
		     call MPI_SEND(Fldx(3:5,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_SEND(Fldy(3:5,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_SEND(Fldz(3:5,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_RECV(Fldx(1:2,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(0,0,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(Fldy(1:2,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(0,0,proczind),5,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(Fldz(1:2,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(0,0,proczind),6,MPI_COMM_WORLD,stat,mpi_err)
		 end if              
    end subroutine ExchangeYZEdgeField_Axis
	
    subroutine ExchangeYZEdgeCurrent_Axis
	      integer :: i,off,my_right,dcount1,mpi_err
	      integer :: stat(MPI_STATUS_SIZE)
	        
		  if(inc_axis.eqv..false.) return
		  if(nSubDomainsX.lt.2) return
		  if(procxind.eq.0) then
			  buff_rJx=0.0_psn; buff_rJy=0.0_psn; buff_rJz=0.0_psn;
			  buff_lJx=0.0_psn; buff_lJy=0.0_psn; buff_lJz=0.0_psn;
	          do i=0,nSubDomainsY-1
				 off=yborders(i)
				 my_right=yborders(i+1)-yborders(i)+5
				 dcount1=3*my_right*mz
		         call MPI_RECV(buff_lJx(:,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
				 call MPI_RECV(buff_lJy(:,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
				 call MPI_RECV(buff_lJz(:,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Jx(mx-2:mx,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Jy(mx-2:mx,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),5,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(Jz(mx-2:mx,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind),6,MPI_COMM_WORLD,stat,mpi_err)	
				 buff_rJx(:,off+1:off+my_right,:)=buff_rJx(:,off+1:off+my_right,:)+buff_lJx(:,off+1:off+my_right,:)
				 buff_rJy(:,off+1:off+my_right,:)=buff_rJy(:,off+1:off+my_right,:)+buff_lJy(:,off+1:off+my_right,:)
				 buff_rJz(:,off+1:off+my_right,:)=buff_rJz(:,off+1:off+my_right,:)+buff_lJz(:,off+1:off+my_right,:)
				 buff_lJx(:,off+1:off+my_right,:)=0.0_psn
				 buff_lJy(:,off+1:off+my_right,:)=0.0_psn
				 buff_lJz(:,off+1:off+my_right,:)=0.0_psn
			 end do 
		 end if 
		 if(procxind.eq.1) then
			 dcount1=3*my*mz
		     call MPI_SEND(Jx(1:3,:,:),dcount1,mpi_psn,proc_grid(0,0,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_SEND(Jy(1:3,:,:),dcount1,mpi_psn,proc_grid(0,0,proczind),2,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_SEND(Jz(1:3,:,:),dcount1,mpi_psn,proc_grid(0,0,proczind),3,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_RECV(buff_lJx,dcount1,mpi_psn,proc_grid(0,0,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(buff_lJy,dcount1,mpi_psn,proc_grid(0,0,proczind),5,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(buff_lJz,dcount1,mpi_psn,proc_grid(0,0,proczind),6,MPI_COMM_WORLD,stat,mpi_err)
			 if(procyind.ne.0) then ! avoid repition in ghost cells during the exchange of XZ and XY edges
				 buff_lJx(:,1:2,:)=0.0_psn
				 buff_lJy(:,1:2,:)=0.0_psn
				 buff_lJz(:,1:2,:)=0.0_psn 
			 end if
			 if(procyind.ne.nSubDomainsY-1) then 
				 buff_lJx(:,my-2:my,:)=0.0_psn
				 buff_lJy(:,my-2:my,:)=0.0_psn
				 buff_lJz(:,my-2:my,:)=0.0_psn
			 end if
		 end if 
         
    end subroutine ExchangeYZEdgeCurrent_Axis
	
    subroutine ExchangeYZEdgeCurrent1_Axis(J1)
         integer :: stat(MPI_STATUS_SIZE)
         integer :: dcount1,mpi_err, i,off,my_short
         real(psn),dimension(mx,my,mz) :: J1
		 
		 if(nSubDomainsX.lt.2) return
 		 if(procxind.eq.0) then
             do i=0,nSubDomainsY-1
				 off=yborders(i)+2
				 if(i.eq.0) off=0
				 my_short=yborders(i+1)-yborders(i) 
				 if(i.eq.0) my_short=my_short+2 
				 if(i.eq.nSubDomainsY-1) my_short=my_short+3
		         dcount1=1*my_short*mz			 
		         call MPI_RECV(J1(mx-2:mx-2,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind),1,MPI_COMM_WORLD,stat,mpi_err)				 
		         call MPI_SEND(J1(mx-3:mx-3,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
			 end do 
		 end if 
		 if(procxind.eq.1) then
	         off=2 
			 if(procyind.eq.0) off=0
			 my_short=my-5 
			 if(procyind.eq.0) my_short=my_short+2
			 if(procyind.eq.nSubDomainsY-1) my_short=my_short+3
			 dcount1=1*my_short*mz		 
		     call MPI_SEND(J1(3:3,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind),1,MPI_COMM_WORLD,stat,mpi_err)
		     call MPI_RECV(J1(2:2,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind),4,MPI_COMM_WORLD,stat,mpi_err)
		 end if 
    end subroutine ExchangeYZEdgeCurrent1_Axis
	
end module cyl_comm_fldprtl