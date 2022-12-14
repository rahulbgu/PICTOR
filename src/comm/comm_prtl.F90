! This module is designed to facilitate transfer of particles between MPI proc.
module comm_prtl 
     use parameters 
     use vars
     use mem_prtl
	 use communication
	 use loadprtlout
     implicit none
contains
           
     subroutine AppendParticles(recv)
		  type(ngbr_list), dimension(:) ::  recv
          integer :: n, m , i 
		  
		  !print*,proc,'np_recv',np_recv,'np',np
		  
		  if(used_prtl_arr_size + np_recv .gt. prtl_arr_size) call ReshapePrtlArr( used_prtl_arr_size + np_recv + 1000000  )
		  
		  do n = 1,size(recv)
			  do m = 1, recv(n)%num_ngbrs
					  
					  call MPI_WAIT(recv(n)%ngbr(m)%mpi_req_prtl,MPI_STATUS_IGNORE)
					  
					  call InsertPrtls(recv(n)%ngbr(m)%p, recv(n)%ngbr(m)%xshift, recv(n)%ngbr(m)%yshift, recv(n)%ngbr(m)%zshift, recv(n)%ngbr(m)%pcount )
					  
			  end do 
		  end do
		  
		       
     end subroutine AppendParticles
	 
	 subroutine ExchangePrtlCPU
		  call ExchangePrtl(qp,xp,yp,zp,up,vp,wp,var1p,flvp,tagp,1,used_prtl_arr_size)
	 end subroutine ExchangePrtlCPU

     subroutine ExchangePrtl(q,x,y,z,u,v,w,var1,flv,tag,min_ind,max_ind)
		  real(psn), dimension(:) :: q,x,y,z,u,v,w,var1
		  integer, dimension(:) :: flv,tag
		  integer :: min_ind, max_ind
		  
		  call StartTimer(9)
          
		  call reset_pcount(ngbr_send)
		  call reset_pcount(ngbr_recv)

		  call LoadPrtlOutliers(q,x,y,z,u,v,w,var1,flv,tag,min_ind,max_ind)
		  
		  call SendPrtlSize(ngbr_send)
		  call RecvPrtlSize(ngbr_recv)
		  
		  if(modulo(t,prtl_reorder_period).eq.0) call update_max_nprtl_out(ngbr_send,max_nprtl_out)
		  
		  call SendPrtl(ngbr_send)
		  
		  call SubstractOutgoingPrtlCount(ngbr_send)
		 
		  call RecvPrtl(ngbr_recv)
		  
		 

		  
          
		  call StopTimer(9)

     end subroutine ExchangePrtl
	 
	 !wait for prtl sends to compelete, so that the buff can be used safely in the next cycle
	 subroutine WaitPrtlSendRecvComplete
	 	  call WaitPrtlComm(ngbr_send)
	 	  call WaitPrtlSizeComm(ngbr_send)
	 end subroutine WaitPrtlSendRecvComplete
	 
	 
	 subroutine SubstractOutgoingPrtlCount(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: m, n, SumPcount 
		 SumPcount = 0 
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 SumPcount = SumPcount + list(n)%ngbr(m)%pcount
			 end do 
		 end do 
		 np = np - SumPcount 	 
	 end subroutine SubstractOutgoingPrtlCount 
	 

 !----------------------------------------------------------------------------------------------
 ! The following subroutines are used for tranferring particles between MPI proc. 
 !----------------------------------------------------------------------------------------------
 
      subroutine SendPrtlSize(send)
 		  type(ngbr_list), dimension(:) :: send	  
           integer :: n, m , ind
 		  do n =1,size(send)
 			  do m=1,send(n)%num_ngbrs
			  
 				  send(n)%ngbr(m)%mpi_req_prtl_size = MPI_REQUEST_NULL
 			  	  call MPI_ISEND(send(n)%ngbr(m)%pcount,1,MPI_INTEGER,send(n)%ngbr(m)%proc, send(n)%ngbr(m)%mpi_tag_offset, MPI_COMM_WORLD, send(n)%ngbr(m)%mpi_req_prtl_size)
							  			  
 			  end do
 		  end do                       
      end subroutine SendPrtlSize
 
      subroutine RecvPrtlSize(recv)
 		  type(ngbr_list), dimension(:) :: recv
           integer :: n, m 
	  
 		  do n =1,size(recv)
 			  do m=1,recv(n)%num_ngbrs
			  
 				  recv(n)%ngbr(m)%pcount = 0 
 				  recv(n)%ngbr(m)%mpi_req_prtl_size = MPI_REQUEST_NULL
 				  call MPI_IRECV(recv(n)%ngbr(m)%pcount,1,MPI_INTEGER,recv(n)%ngbr(m)%proc, recv(n)%ngbr(m)%mpi_tag_offset, MPI_COMM_WORLD, recv(n)%ngbr(m)%mpi_req_prtl_size)
				  
 			  end do
 		  end do                       
      end subroutine RecvPrtlSize
 
      subroutine SendPrtl(send)
 	      type(ngbr_list), dimension(:) :: send
 	      integer :: n, m 
	  
 		  do n =1,size(send)
 			  do m=1,send(n)%num_ngbrs					  

 				  send(n)%ngbr(m)%mpi_req_prtl = MPI_REQUEST_NULL
 				  if(send(n)%ngbr(m)%pcount.gt.0) then 
 					  call MPI_ISEND(send(n)%ngbr(m)%p,send(n)%ngbr(m)%pcount,mpi_prtltype,send(n)%ngbr(m)%proc, send(n)%ngbr(m)%mpi_tag_offset, MPI_COMM_WORLD, send(n)%ngbr(m)%mpi_req_prtl)
 				  end if 
				  				  
 			  end do 
 		  end do 
               
      end subroutine SendPrtl
 
      subroutine RecvPrtl(recv)
 	      type(ngbr_list), dimension(:) :: recv	  
 	      integer :: n, m
	  
 		  np_recv = 0 
	  
 		  do n=1,size(recv)
 			  do m=1,recv(n)%num_ngbrs
		
 			          !Make sure that the size recv is complete
 					  call MPI_WAIT(recv(n)%ngbr(m)%mpi_req_prtl_size,MPI_STATUS_IGNORE)
				  
 					  np_recv = np_recv + recv(n)%ngbr(m)%pcount
				  
 					  !make sure that there is enough space in the buffer
 					  if(size(recv(n)%ngbr(m)%p).lt.recv(n)%ngbr(m)%pcount) then 
 						  deallocate(recv(n)%ngbr(m)%p)
 						  allocate(recv(n)%ngbr(m)%p( int(1.1*recv(n)%ngbr(m)%pcount)+10000))
 					  end if
				  
 				      recv(n)%ngbr(m)%mpi_req_prtl = MPI_REQUEST_NULL
				  
 					  if(recv(n)%ngbr(m)%pcount.gt.0) then 
					  
 						  call MPI_IRECV(recv(n)%ngbr(m)%p,recv(n)%ngbr(m)%pcount,mpi_prtltype,recv(n)%ngbr(m)%proc, recv(n)%ngbr(m)%mpi_tag_offset, MPI_COMM_WORLD, recv(n)%ngbr(m)%mpi_req_prtl)
				  
 					  end if					  


 			  end do 
 		  end do 
               
      end subroutine RecvPrtl
	  
!----------------------------------------------------------------------------------------------
! Wait for send/recv comm. to complete
!----------------------------------------------------------------------------------------------		  
 
 	 subroutine WaitPrtlComm(list)
      	type(ngbr_list), dimension(:) :: list
      	integer :: n, m 
	
 	    do n =1,size(list)
 		    do m=1,list(n)%num_ngbrs
			
 				if(list(n)%ngbr(m)%pcount.gt.0) call MPI_WAIT(list(n)%ngbr(m)%mpi_req_prtl, MPI_STATUS_IGNORE)

 		    end do
 	    end do

 	 end subroutine WaitPrtlComm
 
 	 subroutine WaitPrtlSizeComm(list)
      	type(ngbr_list), dimension(:) :: list
      	integer :: n, m 
	
 	    do n =1,size(list)
 		    do m=1,list(n)%num_ngbrs
			
 				call MPI_WAIT(list(n)%ngbr(m)%mpi_req_prtl_size, MPI_STATUS_IGNORE)

 		    end do
 	    end do

 	 end subroutine WaitPrtlSizeComm
 
 
 !----------------------------------------------------------------------------------------------
 ! set/reset arrays/variables used in communicating particles
 !----------------------------------------------------------------------------------------------	
 
 	 subroutine reset_pcount(list)
 		 type(ngbr_list), dimension(:) :: list
 		 integer :: n, m 
	 
 		 do n=1, size(list)
 			 do m=1,list(n)%num_ngbrs
 				 list(n)%ngbr(m)%pcount=0
 			 end do 
 		 end do 
	 
 	 end subroutine reset_pcount
 
 	 subroutine allocate_comm_prtl_ngbrs(list)
 		 type(ngbr_list), dimension(:) :: list
 		 integer :: n, m 
	 
 		 do n=1, size(list)
 			 do m=1,list(n)%num_ngbrs				 
 				 allocate(list(n)%ngbr(m)%p(max(4*outp_arr_block_size,max_nprtl_out)))
 			 end do 
 		 end do 
	 
 	 end subroutine allocate_comm_prtl_ngbrs
	 
 !----------------------------------------------------------------------------------------------
 ! Periodically update the max. number of particles that are sent to appropiately allocate send buff memory
 !----------------------------------------------------------------------------------------------
 
 	subroutine update_max_nprtl_out(list,npmax)
		type(ngbr_list), dimension(:) :: list
		integer :: npmax
		integer :: n,m
		
		npmax = 0
		do n=1, size(list)
			do m=1,list(n)%num_ngbrs 
				npmax = max(npmax,list(n)%ngbr(m)%pcount)
			end do 
		end do
		npmax = int(1.1*npmax) + 102400
	end subroutine update_max_nprtl_out 	 
	 
 
	 	 
end module comm_prtl