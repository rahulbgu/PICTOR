module comm_loadbalance
     use parameters
     use vars
     use communication
	 use memory
	 implicit none 
     integer :: InsertAt
contains 
     !subroutines in this module help communicate data for the purpose of balancing overall load 
     
!-----------------------------------------------------------------------------------------------
! this following subroutine are used in changing the x-boundaries 
!-----------------------------------------------------------------------------------------------     
     subroutine SendFldToNewProcX(toproc,ledge,redge,Fld,tag)
          integer :: toproc,ledge,redge,tag
          real(psn), dimension(mx,my,mz):: Fld
          integer :: mpi_err
          call MPI_SEND(Fld(ledge:redge,:,:),(redge-ledge+1)*my*mz,mpi_psn,toproc,tag,MPI_COMM_WORLD,mpi_err)
     end subroutine SendFldToNewProcX     
     subroutine RecvFldFromOldProcX(fromproc,ledge,redge,Fld,tag,mx_new,my_new,mz_new)
          integer :: fromproc,ledge,redge,tag
          integer :: mx_new,my_new,mz_new
          real(psn), dimension(mx_new,my_new,mz_new) :: Fld
          integer, dimension(MPI_STATUS_SIZE) :: stat
          integer :: mpi_err
          call MPI_RECV(Fld(ledge:redge,:,:),(redge-ledge+1)*my*mz,mpi_psn,fromproc,tag,MPI_COMM_WORLD,stat,mpi_err)
     end subroutine RecvFldFromOldProcX
     subroutine SendFldToNewProcY(toproc,bedge,tedge,Fld,tag)
          integer :: toproc,bedge,tedge,tag
          real(psn), dimension(mx,my,mz):: Fld
          integer :: mpi_err
          call MPI_SEND(Fld(:,bedge:tedge,:),(tedge-bedge+1)*mx*mz,mpi_psn,toproc,tag,MPI_COMM_WORLD,mpi_err)
     end subroutine SendFldToNewProcY     
     subroutine RecvFldFromOldProcY(fromproc,bedge,tedge,Fld,tag,mx_new,my_new,mz_new)
          integer :: fromproc,bedge,tedge,tag
          integer :: mx_new,my_new,mz_new
          real(psn), dimension(mx_new,my_new,mz_new) :: Fld
          integer, dimension(MPI_STATUS_SIZE) :: stat
          integer :: mpi_err
          call MPI_RECV(Fld(:,bedge:tedge,:),(tedge-bedge+1)*mx*mz,mpi_psn,fromproc,tag,MPI_COMM_WORLD,stat,mpi_err)
     end subroutine RecvFldFromOldProcY
     subroutine SendPrtlSizeToNewProc(proc_no,size)
          integer :: proc_no,size
          integer :: mpi_err
          call MPI_SEND(size,1,MPI_INTEGER,proc_no,1,MPI_COMM_WORLD,mpi_err)
     end subroutine SendPrtlSizeToNewProc
     subroutine RecvPrtlSizeFromOldProc(proc_no,size)
          integer :: proc_no,size
          integer, dimension(MPI_STATUS_SIZE) :: stat  
          integer :: mpi_err   
          call MPI_RECV(size,1,MPI_INTEGER,proc_no,1,MPI_COMM_WORLD,stat,mpi_err)
     end subroutine RecvPrtlSizeFromOldProc
     subroutine SendPrtlToNewProc(proc_no,indi,indf,pcount,prtl_arr,arr_dim)
          integer :: indi,indf,proc_no,pcount,arr_dim,i
          type(particle), dimension(arr_dim) :: prtl_arr
		  integer :: mpi_err	  
          call MPI_SEND(prtl_arr(indi:indf),pcount,mpi_prtltype,proc_no,0,MPI_COMM_WORLD,mpi_err)
     end subroutine SendPrtlToNewProc

	 subroutine RecvPrtlFromOldProc(proc_no,indf,pcount,p1,arr_dim,xshift,yshift)
          integer :: proc_no,indf,pcount,arr_dim,xshift,yshift,n
          type(particle), dimension(arr_dim) :: p1
          integer, dimension(MPI_STATUS_SIZE) :: stat     
          integer :: mpi_err
		  
          call MPI_RECV(p1(1:indf),pcount,mpi_prtltype,proc_no,0,MPI_COMM_WORLD,stat,mpi_err)
          do n=1,indf
		       p1(n)%x=p1(n)%x+xshift
		       p1(n)%y=p1(n)%y+yshift
			   call InsertParticleAt(used_prtl_arr_size+1,p1(n)%x,p1(n)%y,p1(n)%z,p1(n)%u,p1(n)%v,p1(n)%w,p1(n)%q,p1(n)%tag,p1(n)%flv,p1(n)%var1)
			   used_prtl_arr_size = used_prtl_arr_size + 1 				   
          end do
     end subroutine RecvPrtlFromOldProc

! 	 subroutine RecvPrtlFromOldProc(proc_no,indf,pcount,p1,arr_dim,xshift,yshift)
!           integer :: proc_no,indf,pcount,arr_dim,xshift,yshift,n
!           type(particle), dimension(arr_dim) :: p1
!           integer, dimension(MPI_STATUS_SIZE) :: stat
!           integer :: mpi_err
! 		  integer :: InsertAt,LastCallStep
! 		  save InsertAt,LastCallStep
! 		  data LastCallStep /-1/
! 		  if(t.ne.LastCallStep) then
! 			  InsertAt=1
! 			  LastCallStep=t
! 		  end if
!           call MPI_RECV(p1(1:indf),pcount,mpi_prtltype,proc_no,0,MPI_COMM_WORLD,stat,mpi_err)
!           do n=1,indf
! 		       p1(n)%x=p1(n)%x+xshift
! 		       p1(n)%y=p1(n)%y+yshift
! 			   do while(qp(InsertAt).ne.0)
! 				   InsertAt=InsertAt+1
! 			   end do
! 			   call InsertParticleAt(InsertAt,p1(n)%x,p1(n)%y,p1(n)%z,p1(n)%u,p1(n)%v,p1(n)%w,p1(n)%q,p1(n)%tag,p1(n)%flv,p1(n)%var1)
!           end do
!      end subroutine RecvPrtlFromOldProc
	 
	 
     subroutine RecvTestPrtlFromOldProc(proc_no,indf,pcount,p1,arr_dim,xshift,yshift)
          integer :: proc_no,indf,pcount,arr_dim,xshift,yshift,n
          type(particle), dimension(arr_dim) :: p1
          integer, dimension(MPI_STATUS_SIZE) :: stat     
          integer :: mpi_err
		  integer :: InsertAt,LastCallStep 
		  save InsertAt,LastCallStep 
		  data LastCallStep /-1/
		  if(t.ne.LastCallStep) then 
			  InsertAt=1 
			  LastCallStep=t
		  end if
          call MPI_RECV(p1(1:indf),pcount,mpi_prtltype,proc_no,0,MPI_COMM_WORLD,stat,mpi_err)
          do n=1,indf
		       p1(n)%x=p1(n)%x+xshift
		       p1(n)%y=p1(n)%y+yshift
			   do while(qtp(InsertAt).ne.0)
				   InsertAt=InsertAt+1
			   end do 
			   call InsertTestParticleAt(InsertAt,p1(n)%x,p1(n)%y,p1(n)%z,p1(n)%u,p1(n)%v,p1(n)%w,p1(n)%q,p1(n)%tag,p1(n)%flv,p1(n)%var1)				   
          end do
     end subroutine RecvTestPrtlFromOldProc
     

!-----------------------------------------------------------------------------------------------------
!Load balance incase of Homogeneous case
!-----------------------------------------------------------------------------------------------------
subroutine BcastExecTimeAll(arr)
#ifdef twoD	
     real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: arr
#else	 
	 real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: arr
#endif	 
     integer :: mpi_err, size
	 size=nSubDomainsX*nSubDomainsY*nSubDomainsZ
#ifdef twoD 
     size=nSubDomainsX*nSubDomainsY
#endif	 
     call MPI_ALLREDUCE(MPI_IN_PLACE,arr,size,MPI_REAL,mpi_sum,MPI_COMM_WORLD,mpi_err) !synch the matrix on all proc
end subroutine BcastExecTimeAll

subroutine ReduceNumPrtlX(arr, i1, i2, size)
	integer  :: i1, i2, size, mpi_err
	real(dbpsn), dimension(i1:i2) :: arr
	call MPI_ALLREDUCE(MPI_IN_PLACE,arr(i1:i2),size,MPI_REAL8,mpi_sum,MPI_COMM_WORLD,mpi_err)
end subroutine ReduceNumPrtlX


     
end module comm_loadbalance
          