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
     subroutine ResetInsertIndex
		 InsertAt=1 
	 end subroutine ResetInsertIndex 
	 subroutine RecvPrtlFromOldProc(proc_no,indf,pcount,p1,arr_dim,xshift,yshift)
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
			   do while(qp(InsertAt).ne.0)
				   InsertAt=InsertAt+1
			   end do 
			   call InsertParticleAt(InsertAt,p1(n)%x,p1(n)%y,p1(n)%z,p1(n)%u,p1(n)%v,p1(n)%w,p1(n)%q,p1(n)%tag,p1(n)%flv,p1(n)%var1)				   
          end do
     end subroutine RecvPrtlFromOldProc
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
     
!-------------------------------------------------------------------------------------------------------
! subroutines used for swapping physical domain
!-------------------------------------------------------------------------------------------------------
subroutine Swap_SendRecvFldSize(proc_no,mx_new,my_new,mz_new)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: proc_no,mx_new,my_new,mz_new
     integer, dimension(3) :: SendFldSize,RecvFldSize
     integer :: mpi_err
     SendFldSize(1)=mx
     SendFldSize(2)=my
     SendFldSize(3)=mz
     call MPI_SENDRECV(SendFldSize,3,MPI_INTEGER,proc_no,proc,RecvFldSize,3,MPI_INTEGER,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)
     mx_new=RecvFldSize(1)
     my_new=RecvFldSize(2)
     mz_new=RecvFldSize(3)
end subroutine Swap_SendRecvFldSize     
subroutine Swap_SendRecvFld(proc_no,Fld,mx_new,my_new,mz_new)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: proc_no,mx_new,my_new,mz_new
     real(psn), dimension(:,:,:), allocatable :: Fld,SwapFldTemp
     integer :: mpi_err
     allocate(SwapFldTemp(mx_new,my_new,mz_new))
     call MPI_SENDRECV(Fld,mx*my*mz,mpi_psn,proc_no,proc,SwapFldTemp,mx_new*my_new*mz_new,mpi_psn,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)
     call move_alloc(SwapFldTemp,Fld)
end subroutine Swap_SendRecvFld
subroutine Swap_SendRecvProcXYind(proc_no)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: proc_no
     integer, dimension(2) :: SendXYind,RecvXYind
     integer :: mpi_err
     SendXYind(1)=procxind(proc)
     SendXYind(2)=procyind(proc)
     call MPI_SENDRECV(SendXYind,2,MPI_INTEGER,proc_no,proc,RecvXYind,2,MPI_INTEGER,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)
     procxind(proc)=RecvXYind(1)
     procyind(proc)=RecvXYind(2)
end subroutine Swap_SendRecvProcXYind

subroutine Swap_SendRecvPrtlSize(proc_no,new_prtl_arr_size,new_np)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: proc_no,new_prtl_arr_size,new_np
     integer,dimension(2) :: SendPrtlSize,RecvPrtlSize
     integer :: mpi_err
     SendPrtlSize(1)=prtl_arr_size
     SendPrtlSize(2)=np
     call MPI_SENDRECV(SendPrtlSize,2,MPI_INTEGER,proc_no,proc,RecvPrtlSize,2,MPI_INTEGER,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)
     new_prtl_arr_size=RecvPrtlSize(1)
     new_np=RecvPrtlSize(2)
end subroutine Swap_SendRecvPrtlSize
subroutine Swap_SendRecvPrtl(proc_no,SwapPrtl,new_size)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: proc_no,new_size
     type(particle), dimension(new_size) :: SwapPrtl 
     integer :: mpi_err
     call MPI_SENDRECV(p,prtl_arr_size,mpi_prtltype,proc_no,proc,SwapPrtl,new_size,mpi_prtltype,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)
end subroutine Swap_SendRecvPrtl
subroutine Swap_SendRecvNeighbourProc(proc_no,ngb_proc)
     integer :: stat(MPI_STATUS_SIZE)
     integer :: ngb_proc, proc_no,ngb_proc_temp
     integer :: mpi_err
     call MPI_SENDRECV(ngb_proc,1,MPI_INTEGER,proc_no,proc,ngb_proc_temp,1,MPI_INTEGER,proc_no,proc_no,MPI_COMM_WORLD,stat,mpi_err)     
     ngb_proc=ngb_proc_temp
end subroutine Swap_SendRecvNeighbourProc
!---------------end of subroutines used for swapping physical domain ---------------------------------

!-----------------------------------------------------------------------------------------------------
!Load balance incase of Homogeneous case
!-----------------------------------------------------------------------------------------------------
subroutine BcastExecTimeAllHomogeneous(arr)
     real, dimension(0:nSubDomainsX-1,0:nSubDomainsY-1), intent(inout) :: arr
     integer :: mpi_err
     call MPI_ALLREDUCE(MPI_IN_PLACE,arr,nproc,MPI_REAL,mpi_sum,MPI_COMM_WORLD,mpi_err) !synch the matrix on all proc
end subroutine BcastExecTimeAllHomogeneous


     
end module comm_loadbalance
          