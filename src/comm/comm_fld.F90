module comm_fld
	use parameters
	use vars
    use communication
	 
    implicit none
                                        
contains 
		 
!----------------------------------------------------------------------------------------------
! Exchange fld data
!----------------------------------------------------------------------------------------------	 
	 
	 
	 subroutine UpdateCurrentAllFaces(Fldx,Fldy,Fldz)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz

! 		 call SendRecvCurrentYZ(Fldx,Fldy,Fldz)
! 		 call SendRecvCurrentZX(Fldx,Fldy,Fldz)
! #ifndef twoD
!        call SendRecvCurrentXY(Fldx,Fldy,Fldz)
! #endif


		 
		 if(indepLBaxis.eq.0) call SendRecvCurrentYZ(Fldx,Fldy,Fldz)
		 if(indepLBaxis.eq.1) call SendRecvCurrentZX(Fldx,Fldy,Fldz)
#ifndef twoD
         if(indepLBaxis.eq.2) call SendRecvCurrentXY(Fldx,Fldy,Fldz)
#endif


		 if(indepLBaxis.ne.0) call SendRecvCurrentYZ(Fldx,Fldy,Fldz)
		 if(indepLBaxis.ne.1) call SendRecvCurrentZX(Fldx,Fldy,Fldz)
#ifndef twoD
         if(indepLBaxis.ne.2) call SendRecvCurrentXY(Fldx,Fldy,Fldz)
#endif



	 end subroutine UpdateCurrentAllFaces
	 
	 subroutine SendRecvCurrentYZ(Fldx,Fldy,Fldz)
		 real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
		 call WaitCommBuffFlds(ngbr_recv(1:2))
		 call WaitCommBuffFlds(ngbr_send(1:2))
		 call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_recv(1:2))
		 call SendBuffFlds(ngbr_recv(1:2))
		 call RecvBuffFlds(ngbr_send(1:2))
 	     call WaitCommBuffFlds(ngbr_send(1:2))
		 call AddBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_send(1:2))
	
	end subroutine SendRecvCurrentYZ
	
	subroutine SendRecvCurrentZX(Fldx,Fldy,Fldz)
		 real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
		 call WaitCommBuffFlds(ngbr_recv(3:4))
		 call WaitCommBuffFlds(ngbr_send(3:4))
		 call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_recv(3:4))
		 call SendBuffFlds(ngbr_recv(3:4))
		 call RecvBuffFlds(ngbr_send(3:4))
 	     call WaitCommBuffFlds(ngbr_send(3:4))
		 call AddBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_send(3:4))
	
	end subroutine SendRecvCurrentZX
	
	subroutine SendRecvCurrentXY(Fldx,Fldy,Fldz)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
         call WaitCommBuffFlds(ngbr_recv(5:6))
         call WaitCommBuffFlds(ngbr_send(5:6))
         call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_recv(5:6))
		 call SendBuffFlds(ngbr_recv(5:6))
		 call RecvBuffFlds(ngbr_send(5:6))
 	     call WaitCommBuffFlds(ngbr_send(5:6))
		 call AddBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_send(5:6))
		
	end subroutine SendRecvCurrentXY
		 
	 
	 subroutine SendRecvFlds(Fldx,Fldy,Fldz)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz

! #ifndef twoD
! 	     call SendRecvFldsXY(Fldx,Fldy,Fldz)
! #endif
!  		 call SendRecvFldsZX(Fldx,Fldy,Fldz)
!        call SendRecvFldsYZ(Fldx,Fldy,Fldz)


		 if(indepLBaxis.ne.0) call SendRecvFldsYZ(Fldx,Fldy,Fldz)
		 if(indepLBaxis.ne.1) call SendRecvFldsZX(Fldx,Fldy,Fldz)
#ifndef twoD
		 if(indepLBaxis.ne.2) call SendRecvFldsXY(Fldx,Fldy,Fldz)
#endif

	 	 if(indepLBaxis.eq.0) call SendRecvFldsYZ(Fldx,Fldy,Fldz)
		 if(indepLBaxis.eq.1) call SendRecvFldsZX(Fldx,Fldy,Fldz)
#ifndef twoD
		 if(indepLBaxis.eq.2) call SendRecvFldsXY(Fldx,Fldy,Fldz)
#endif
	 
	 end subroutine SendRecvFlds
	 
	 subroutine SendRecvFldsYZ(Fldx,Fldy,Fldz)
		 real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
		 call WaitCommBuffFlds(ngbr_send(1:2))
		 call WaitCommBuffFlds(ngbr_recv(1:2))
		 call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_send(1:2))
		 call SendBuffFlds(ngbr_send(1:2))
		 call RecvBuffFlds(ngbr_recv(1:2))
 	     call WaitCommBuffFlds(ngbr_recv(1:2))
		 call CopyBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_recv(1:2))
		 
	 end subroutine SendRecvFldsYZ
	 
	 subroutine SendRecvFldsZX(Fldx,Fldy,Fldz)
		 real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
		 call WaitCommBuffFlds(ngbr_send(3:4))
		 call WaitCommBuffFlds(ngbr_recv(3:4))
		 call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_send(3:4))
		 call SendBuffFlds(ngbr_send(3:4))
		 call RecvBuffFlds(ngbr_recv(3:4))
 	     call WaitCommBuffFlds(ngbr_recv(3:4))
		 call CopyBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_recv(3:4))
		 
	 end subroutine SendRecvFldsZX
	 
	 subroutine SendRecvFldsXY(Fldx,Fldy,Fldz)
		 real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 
         call WaitCommBuffFlds(ngbr_send(5:6))
         call WaitCommBuffFlds(ngbr_recv(5:6))
         call CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,ngbr_send(5:6))
		 call SendBuffFlds(ngbr_send(5:6))
		 call RecvBuffFlds(ngbr_recv(5:6))
 	     call WaitCommBuffFlds(ngbr_recv(5:6))
		 call CopyBuffToFldsNgbrs(Fldx,Fldy,Fldz,ngbr_recv(5:6))
		 
	 end subroutine SendRecvFldsXY
	 
		 
	 !---------------------------------------------------------------------------------------
	 ! Wait for Recv operaton to complete
	 !---------------------------------------------------------------------------------------
	 
	 subroutine WaitCommBuffFlds(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m 
		 
		 do n =1, size(list)
			 do m=1,list(n)%num_ngbrs				 
				  call MPI_WAIT(list(n)%ngbr(m)%mpi_req_fldx,MPI_STATUS_IGNORE)
				  call MPI_WAIT(list(n)%ngbr(m)%mpi_req_fldy,MPI_STATUS_IGNORE)
				  call MPI_WAIT(list(n)%ngbr(m)%mpi_req_fldz,MPI_STATUS_IGNORE)
			 end do 
		 end do 
		 		 
	 end subroutine WaitCommBuffFlds
		 
	 !---------------------------------------------------------------------------------------
	 ! Send and Recv buffer memory and return immediately (non-blocking)
	 !---------------------------------------------------------------------------------------
	 	
	 subroutine SendBuffFlds(list)
		 type(ngbr_list), dimension(:) :: list
 		 integer :: n,m, dcount

		 do n =1, size(list)
			 do m=1,list(n)%num_ngbrs
				 				 
				 list(n)%ngbr(m)%mpi_req_fldx = MPI_REQUEST_NULL
				 list(n)%ngbr(m)%mpi_req_fldy = MPI_REQUEST_NULL	
				 list(n)%ngbr(m)%mpi_req_fldz = MPI_REQUEST_NULL					 
				 dcount = (list(n)%ngbr(m)%ind(2) - list(n)%ngbr(m)%ind(1) +1)*(list(n)%ngbr(m)%ind(4) - list(n)%ngbr(m)%ind(3) +1)*(list(n)%ngbr(m)%ind(6) - list(n)%ngbr(m)%ind(5) +1)
				 call MPI_ISEND( list(n)%ngbr(m)%Fldx, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +1 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldx )
				 call MPI_ISEND( list(n)%ngbr(m)%Fldy, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +2 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldy )
				 call MPI_ISEND( list(n)%ngbr(m)%Fldz, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +3 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldz )
				 
			 end do
		 end do

	 end subroutine SendBuffFlds

	 subroutine RecvBuffFlds(list)
		 type(ngbr_list), dimension(:) :: list
 		 integer :: n,m, dcount

		 do n =1, size(list)
			 do m=1,list(n)%num_ngbrs				 
				 
				 list(n)%ngbr(m)%mpi_req_fldx = MPI_REQUEST_NULL
				 list(n)%ngbr(m)%mpi_req_fldy = MPI_REQUEST_NULL	
				 list(n)%ngbr(m)%mpi_req_fldz = MPI_REQUEST_NULL	
				 dcount = (list(n)%ngbr(m)%ind(2) - list(n)%ngbr(m)%ind(1) +1)*(list(n)%ngbr(m)%ind(4) - list(n)%ngbr(m)%ind(3) +1)*(list(n)%ngbr(m)%ind(6) - list(n)%ngbr(m)%ind(5) +1)
				 call MPI_IRECV( list(n)%ngbr(m)%Fldx, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +1 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldx )
				 call MPI_IRECV( list(n)%ngbr(m)%Fldy, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +2 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldy )
				 call MPI_IRECV( list(n)%ngbr(m)%Fldz, dcount, mpi_psn, list(n)%ngbr(m)%proc, list(n)%ngbr(m)%mpi_tag_offset +3 , MPI_COMM_WORLD, list(n)%ngbr(m)%mpi_req_fldz )
	
			 end do
		 end do

	 end subroutine RecvBuffFlds
	 	 
	 !---------------------------------------------------------------------------------------
	 ! Copy To/From buffer memory 
	 !---------------------------------------------------------------------------------------
	 
	 subroutine CopyFldsToBuffNgbrs(Fldx,Fldy,Fldz,list)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 if(list(n)%ngbr(m)%proc.eq.MPI_PROC_NULL) cycle
				 call CopyFldToBuff(Fldx,list(n)%ngbr(m)%Fldx,list(n)%ngbr(m)%ind)
				 call CopyFldToBuff(Fldy,list(n)%ngbr(m)%Fldy,list(n)%ngbr(m)%ind)
				 call CopyFldToBuff(Fldz,list(n)%ngbr(m)%Fldz,list(n)%ngbr(m)%ind)
			 end do
		 end do  	 
	 end subroutine CopyFldsToBuffNgbrs
	 
	 subroutine CopyBuffToFldsNgbrs(Fldx,Fldy,Fldz,list)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m 
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 if(list(n)%ngbr(m)%proc.eq.MPI_PROC_NULL) cycle
				 call CopyBuffToFld(Fldx,list(n)%ngbr(m)%Fldx,list(n)%ngbr(m)%ind)
				 call CopyBuffToFld(Fldy,list(n)%ngbr(m)%Fldy,list(n)%ngbr(m)%ind)
				 call CopyBuffToFld(Fldz,list(n)%ngbr(m)%Fldz,list(n)%ngbr(m)%ind)
			 end do
		 end do  	 
	 end subroutine CopyBuffToFldsNgbrs
		 
	 
	 subroutine CopyFldToBuff(Fld,buff,ind)
		 real(psn), dimension(:,:,:) :: Fld, buff 
		 integer, dimension(:) :: ind
		 integer :: sx,sy,sz
		 sx=ind(2)-ind(1) +1 
		 sy=ind(4)-ind(3) +1
		 sz=ind(6)-ind(5) +1 
		 buff(1:sx,1:sy,1:sz) = Fld(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) 
	 end subroutine CopyFldToBuff
	 
	 subroutine CopyBuffToFld(Fld,buff,ind)
		 real(psn), dimension(:,:,:) :: Fld, buff 
		 integer, dimension(:) :: ind
		 integer :: sx,sy,sz
		 sx=ind(2)-ind(1) +1 
		 sy=ind(4)-ind(3) +1
		 sz=ind(6)-ind(5) +1 
		 Fld(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) =  buff(1:sx,1:sy,1:sz)
	 end subroutine CopyBuffToFld
	 
	 !---------------------------------------------------------------------------------------
	 ! Add buffer memory array to the main fld arrays
	 !---------------------------------------------------------------------------------------
	 
	 subroutine AddBuffToFldsNgbrs(Fldx,Fldy,Fldz,list)
         real(psn), dimension(:,:,:) :: Fldx,Fldy,Fldz
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m 
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 if(list(n)%ngbr(m)%proc.eq.MPI_PROC_NULL) cycle
				 call AddBuffToFld(Fldx,list(n)%ngbr(m)%Fldx,list(n)%ngbr(m)%ind)
				 call AddBuffToFld(Fldy,list(n)%ngbr(m)%Fldy,list(n)%ngbr(m)%ind)
				 call AddBuffToFld(Fldz,list(n)%ngbr(m)%Fldz,list(n)%ngbr(m)%ind)
			 end do 
		 end do 	 
	 end subroutine AddBuffToFldsNgbrs
	 
	 subroutine AddBuffToFld(Fld,buff,ind)
		 real(psn), dimension(:,:,:) :: Fld, buff 
		 integer, dimension(:), intent(IN) :: ind
		 integer :: sx,sy,sz
		 sx=ind(2)-ind(1) +1 
		 sy=ind(4)-ind(3) +1
		 sz=ind(6)-ind(5) +1 
		 Fld(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) = Fld(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6)) + buff(1:sx,1:sy,1:sz)
	 end subroutine AddBuffToFld
	 
     
end module comm_fld