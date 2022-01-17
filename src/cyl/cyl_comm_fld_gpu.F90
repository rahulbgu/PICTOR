module cyl_comm_fld_gpu
	use parameters
	use vars
	use cyl_vars
	use communication
	implicit none
contains
    subroutine ExchangeYZEdgeField_Axis_GPU(lFld_send,lFld_recv,rFld_send,rFld_recv)
		 real(psn), dimension(3,my,mz)    :: lFld_send,rFld_recv
		 real(psn), dimension(2,my,mz)    :: rFld_send,lFld_recv
         integer :: i,off,my_short,dcount1,dcount2,mpi_err
         integer :: stat(MPI_STATUS_SIZE)
         dcount1=3*my*mz
         dcount2=2*my*mz
 		 if(inc_axis.eqv..false.) return
	     if(nSubDomainsX.lt.2) return
 		 if(procxind(proc).eq.0) then
             do i=0,nSubDomainsY-1
				 off=yborders(i)+2
				 if(i.eq.0) off=0
				 my_short=yborders(i+1)-yborders(i) 
				 if(i.eq.0) my_short=my_short+2 
				 if(i.eq.nSubDomainsY-1) my_short=my_short+3
		         dcount1=3*my_short*mz
		         dcount2=2*my_short*mz			 
		         call MPI_RECV(rFld_recv(:,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(1,i,proczind(proc)),1,MPI_COMM_WORLD,stat,mpi_err)				 
		         call MPI_SEND(rFld_send(:,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(1,i,proczind(proc)),4,MPI_COMM_WORLD,stat,mpi_err)
			 end do 
		 end if 
		 if(procxind(proc).eq.1) then
	         off=2 
			 if(procyind(proc).eq.0) off=0
			 my_short=my-5 
			 if(procyind(proc).eq.0) my_short=my_short+2
			 if(procyind(proc).eq.nSubDomainsY-1) my_short=my_short+3
			 dcount1=3*my_short*mz
	         dcount2=2*my_short*mz	 
		     call MPI_SEND(lFld_send(:,off+1:off+my_short,:),dcount1,mpi_psn,proc_grid(0,0,proczind(proc)),1,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(lFld_recv(:,off+1:off+my_short,:),dcount2,mpi_psn,proc_grid(0,0,proczind(proc)),4,MPI_COMM_WORLD,stat,mpi_err)
		 end if              
    end subroutine ExchangeYZEdgeField_Axis_GPU
	
    subroutine ExchangeYZEdgeCurrent_Axis_GPU(lFld_send,lFld_recv,rFld_send,rFld_recv)
	      real(psn), dimension(3,my,mz)    :: lFld_send,rFld_recv,rFld_send,lFld_recv
	      integer :: i,off,my_right,dcount1,mpi_err
	      integer :: stat(MPI_STATUS_SIZE)
	        
		  if(inc_axis.eqv..false.) return
		  if(nSubDomainsX.lt.2) return
		  if(procxind(proc).eq.0) then
              rFld_recv=0.0_psn 
			  lFld_recv=0.0_psn
	          do i=0,nSubDomainsY-1
				 off=yborders(i)
				 my_right=yborders(i+1)-yborders(i)+5
				 dcount1=3*my_right*mz
		         call MPI_RECV(lFld_recv(:,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind(proc)),1,MPI_COMM_WORLD,stat,mpi_err)
		         call MPI_SEND(rFld_send(:,off+1:off+my_right,:),dcount1,mpi_psn,proc_grid(1,i,proczind(proc)),4,MPI_COMM_WORLD,stat,mpi_err)	
				 rFld_recv(:,off+1:off+my_right,:)=rFld_recv(:,off+1:off+my_right,:)+lFld_recv(:,off+1:off+my_right,:)
				 lFld_recv(:,off+1:off+my_right,:)=0.0_psn
			 end do 
		 end if 
		 if(procxind(proc).eq.1) then
			 dcount1=3*my*mz
		     call MPI_SEND(lFld_send(1:3,:,:),dcount1,mpi_psn,proc_grid(0,0,proczind(proc)),1,MPI_COMM_WORLD,stat,mpi_err)
			 call MPI_RECV(lFld_recv,dcount1,mpi_psn,proc_grid(0,0,proczind(proc)),4,MPI_COMM_WORLD,stat,mpi_err)
			 if(procyind(proc).ne.0) then ! avoid repetition in ghost cells during the exchange of XZ and XY edges
				 lFld_recv(:,1:2,:)=0.0_psn
			 end if
			 if(procyind(proc).ne.nSubDomainsY-1) then 
				 lFld_recv(:,my-2:my,:)=0.0_psn
			 end if
		 end if 
         
    end subroutine ExchangeYZEdgeCurrent_Axis_GPU

end module cyl_comm_fld_gpu