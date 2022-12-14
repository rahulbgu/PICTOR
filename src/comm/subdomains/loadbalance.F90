module loadbalance
     use parameters
     use vars
     use mem_prtl
	 use movdep
	 use initialise
	 use mem_fld
	 use subdomains
	 use comm_fld
	 use loadprtlout
	 use comm_prtl
#ifdef gpu	 
	 use loadbalance_gpu
	 use fields_gpu
	 use particles_gpu
	 use initialise_gpu
#endif	 
	 use loadprtlout 
     implicit none 
      
contains 
	
!---------------------------------------------------------------------------
! Following subroutines are for computing and using "npx"; particle density along the x-direction (cell-by-cell)
!---------------------------------------------------------------------------	
    subroutine BalanceLoad
		if(modulo(t,load_balance_period).ne.0) return

		call AdjustIndepBorders

	end subroutine BalanceLoad 
	
	subroutine NumPrtl1D(np1d_local,pos)
		real(dbpsn), dimension(:) ::  np1d_local
		real(psn), dimension(:) :: pos
		integer :: n, ind
		np1d_local=0.01
#ifdef gpu
        !call CalcNumPrtlX_GPU( np1d_local,shift)
#else
		do n=1,used_prtl_arr_size
			if(flvp(n).eq.0) cycle
			ind=pos(n)
			np1d_local(ind)= np1d_local(ind)+1
		end do
#endif		
	end subroutine NumPrtl1D
	
	subroutine NumPrtl1D_Global(np1d,np1d_local,bsize,borders,procind)
		integer :: bsize, procind
		integer, dimension(0:bsize) :: borders
		real(dbpsn), dimension(:) :: np1d, np1d_local
		integer :: i,j, jmax
		
		np1d = 0
		jmax = borders(bsize) - borders(0)
		do i=1,size(np1d_local)
			j= i -3 +1  + ( borders(procind) - borders(0) ) 
			if(j.ge.1 .and. j.le.jmax) np1d(j) = np1d_local(i)		
		end do 
	end subroutine NumPrtl1D_Global
	
	
    !--------------------------------------------------------------------------
	! Determine new borders to balance the load
	!--------------------------------------------------------------------------
	subroutine AdjustIndepBorders
		real(dbpsn), dimension(:), allocatable :: np1d_local, np1d
		integer, dimension(:), allocatable :: borders_new
		integer :: proc_count
		
		proc_count = ProcGrid(iproc,jproc)%count	
		allocate(np1d_local( ProcGrid(iproc,jproc)%borders(kproc+1) - ProcGrid(iproc,jproc)%borders(kproc) +5 ))
		allocate(np1d( ProcGrid(iproc,jproc)%borders(proc_count) - ProcGrid(iproc,jproc)%borders(0)))
		
		if(indepLBaxis.eq.0) call NumPrtl1D(np1D_local,xp)
		if(indepLBaxis.eq.1) call NumPrtl1D(np1D_local,yp)
		if(indepLBaxis.eq.2) call NumPrtl1D(np1D_local,zp)
		
		call NumPrtl1D_Global( np1d, np1d_local, proc_count, ProcGrid(iproc,jproc)%borders, kproc )
		
		
		call MPI_ALLREDUCE(MPI_IN_PLACE,np1d, ProcGrid(iproc,jproc)%borders(proc_count) - ProcGrid(iproc,jproc)%borders(0) , MPI_DOUBLE_PRECISION, mpi_sum, comm_indepLBaxis)

		if(proc.eq.0) print*,proc,'before borders are:',ProcGrid(iproc,jproc)%borders,'time',t

		call DetermineNewBorders(ProcGrid(iproc,jproc)%count,ProcGrid(iproc,jproc)%borders,borders_new,np1d,kproc)
		
		!print*,proc,'new borders are:',borders_new,'time',t
!>>>>>-----------		
!if(procyind.eq.0) borders_new(1) = 12
		!borders_new(1) =32
		!borders_new(1) = 132!borders_new(1) +4
		
		if( borders_changed(ProcGrid(iproc,jproc)%count,ProcGrid(iproc,jproc)%borders,borders_new) ) then
			 
			call SetNewBorders(ProcGrid(iproc,jproc)%count,ProcGrid(iproc,jproc)%borders,ProcGrid(iproc,jproc)%count,borders_new)
			
			call SetSubDomainGridSize
			call SetPrtlBoundaries
			call InitAuxFld
		
		end if
				
	    call SetNgbr
        call SetSendRecvFldDomain
        call SetSendRecvPrtl
		call SetBoundaryNgbr
		
		!reset MPI groups
		call FreeProcGridCommGroup
		call SetProcGridCommGroup
		
				
		!update ghost cells
		call SendRecvFlds(Ex,Ey,Ez)
		call SendRecvFlds(Bx,By,Bz)
		
		
		
		if(proc.eq.0) print*,proc,'new borders are:',borders_new,'time',t
		
	end subroutine AdjustIndepBorders
	
	subroutine DetermineNewBorders(bsize,borders,borders_new, np1d ,procind)
		integer :: bsize, procind
		integer, dimension(0:bsize) :: borders 
		integer, dimension(:), allocatable :: borders_new
		real(dbpsn), dimension(:) :: np1d
		real(dbpsn) :: np1d_subdomain, np1d_sum
		integer :: n , k
		
		np1d_subdomain = sum(np1d)/bsize
		
		allocate(borders_new(0:bsize))
		
		np1d_sum= 0.0_dbpsn

	    borders_new(0) = borders(0)
		k=1
		do n=1,size(np1d)
			if(np1d_sum + np1d(n).ge.k*np1d_subdomain) then 
				borders_new(k) = max(borders_new(k-1)+4,borders_new(0)+n) ! +4:subdomains are atleast 4 cells wide 				
				k=k+1
			end if 
			if( n .ge. borders(bsize) - 4*(bsize-k) ) then ! to ensure that next borders defined in the loop are atleast 4 cells wide 
				borders_new(k) = borders(bsize) - 4*(bsize-k)
				k=k+1
			end if
			np1d_sum = np1d_sum + np1d(n)
		end do 
	
		borders_new(bsize)=borders(bsize)
		
		  
	end subroutine DetermineNewBorders
	
	logical function borders_changed(size,borders,borders_new)
		integer :: size
		integer, dimension(0:size) :: borders, borders_new
		integer :: n
		
		borders_changed = .false.
		do n = 1,size-1
			if(abs(borders(n)-borders_new(n)).ge.4) then ! borders are changed only if shifts are significant (>= 4)
				borders_changed = .true.
			end if 
		end do  
	
	end function borders_changed
	
	
    !--------------------------------------------------------------------------
	! Exchange data and resize arrays to set new borders along axisProcGrid direction of the ProcGrid array
	!--------------------------------------------------------------------------
	subroutine SetNewBorders(size,borders,size_new,borders_new)
		integer :: size, size_new
		integer, dimension(0:size) :: borders
		integer, dimension(0:size_new) :: borders_new
				
		! deallocate older send/recv buff, they will be reallocated later 
	    if(allocated(ngbr_send)) then 
			call DeallocateNgbrList(ngbr_send)
			deallocate(ngbr_send)
	    end if 
	    if(allocated(ngbr_recv)) then 
			call DeallocateNgbrList(ngbr_recv)
			deallocate(ngbr_recv)
	    end if 
		
		!determine proc. and domains(edges/shifts) for the loadbalance send/recv 
		allocate(lb_send(1),lb_recv(1))
		call SetCommSkeltonLB(size,borders,size_new,borders_new)
		
	
		!load particles in the buffer and send/recv
	    call reset_pcount(lb_send) 
		call ScanPrtlOutlierLB(indepLBaxis,lb_send(1),lb_send(1)%edges)
		
		!set prtl buffer memomry
		call allocate_comm_prtl_ngbrs(lb_recv)
		call allocate_comm_prtl_ngbrs_LB(lb_send)
		
		call reset_pcount(lb_send) 
		call LoadPrtlOutlierLB(indepLBaxis,lb_send(1),lb_send(1)%edges)
		
		
		call reset_pcount(lb_recv) 
	    
		!send/recv prtl data		
	    call SendPrtlSize(lb_send)
	    call RecvPrtlSize(lb_recv)
	    call SendPrtl(lb_send)
		call RecvPrtl(lb_recv)
		
			
		
		!Resize electric field 
		call WaitCommBuffFlds(lb_send)
		call WaitCommBuffFlds(lb_recv)
		call CopyFldsToBuffNgbrs(Ex,Ey,Ez,lb_send)		
		call SendBuffFlds(lb_send)
		call ResizeFldsLB(Ex,size,borders,borders_new,indepLBaxis,kproc)
		call ResizeFldsLB(Ey,size,borders,borders_new,indepLBaxis,kproc)
		call ResizeFldsLB(Ez,size,borders,borders_new,indepLBaxis,kproc)
		call RecvBuffFlds(lb_recv)
		call WaitCommBuffFlds(lb_recv)
		call CopyBuffToFldsNgbrs(Ex,Ey,Ez,lb_recv)


		!Resize magnetic field
		call WaitCommBuffFlds(lb_send)
		call WaitCommBuffFlds(lb_recv)		
		call CopyFldsToBuffNgbrs(Bx,By,Bz,lb_send)
		call SendBuffFlds(lb_send)
		call ResizeFldsLB(Bx,size,borders,borders_new,indepLBaxis,kproc)
		call ResizeFldsLB(By,size,borders,borders_new,indepLBaxis,kproc)
		call ResizeFldsLB(Bz,size,borders,borders_new,indepLBaxis,kproc)
		call RecvBuffFlds(lb_recv)
		call WaitCommBuffFlds(lb_recv)
		call CopyBuffToFldsNgbrs(Bx,By,Bz,lb_recv)
		
		
		!append the incoming particles into the main particle arrays
		call AppendParticles(lb_recv)
		!call MPI_BARRIER(MPI_COMM_WORLD)
				
		
		!update indep. borders
		if(indepLBaxis.eq.0) xborders = borders_new
		if(indepLBaxis.eq.1) yborders = borders_new
		if(indepLBaxis.eq.2) zborders = borders_new
		ProcGrid(iproc,jproc)%borders = borders_new
		
		!wait for the data trasfer completion in loadbalacning 
		call WaitCommBuffFlds(lb_send)
		call WaitPrtlSizeComm(lb_send)
		call WaitPrtlComm(lb_send)
				
		!deallocate the memory used in load-balancing
	    if(allocated(lb_send)) then 
			call DeallocateNgbrList(lb_send)
			deallocate(lb_send)
	    end if 
	    if(allocated(lb_recv)) then 
			call DeallocateNgbrList(lb_recv)
			deallocate(lb_recv)
	    end if 
		
	end subroutine SetNewBorders
		
	
	subroutine SetCommSkeltonLB(size,borders,size_new,borders_new)
		integer :: size, size_new
		integer, dimension(0:size) :: borders
		integer, dimension(0:size_new) :: borders_new
				
		
		call MultipleNgbr( size, borders, size_new , borders_new, ProcGrid(iproc,jproc)%procs , lb_send(1))
		call MultipleNgbr( size_new, borders_new, size, borders, ProcGrid(iproc,jproc)%procs , lb_recv(1))
		
		call SetSendRecvFldDomainLB(lb_send)
		call SetSendRecvFldDomainLB(lb_recv)
		
		call SetSendRecvPrtlOffsets(lb_send)
		call SetSendRecvPrtlOffsets(lb_recv)
				
	end subroutine SetCommSkeltonLB
	

	
	subroutine SetSendRecvFldDomainLB(list)
		type(ngbr_list), dimension(:) :: list
		integer :: m , e1 ,e2
	
		
		do m=1,list(1)%num_ngbrs
			list(1)%ngbr(m)%ind(1)=1; list(1)%ngbr(m)%ind(2)=mx; list(1)%ngbr(m)%ind(3)=1; list(1)%ngbr(m)%ind(4)=my; list(1)%ngbr(m)%ind(5)=1; list(1)%ngbr(m)%ind(6)=mz;
									
			e1 = list(1)%edges(m)
			e2 = list(1)%edges(m+1)

			if(m.ne.1)                 e1 = e1 -2
			if(m.ne.list(1)%num_ngbrs) e2 = e2 +2
			
			list(1)%ngbr(m)%ind(1+2*indepLBaxis) = e1
			list(1)%ngbr(m)%ind(2+2*indepLBaxis) = e2
			
			list(1)%ngbr(m)%mpi_tag_offset = 0 
			
		end do 
		
		call AllocCommFld(list)
	end subroutine SetSendRecvFldDomainLB
	
			
	!-------------------------------------------------------------------------
	! resize fld arrays and shift old data
	!-------------------------------------------------------------------------
	
	subroutine ResizeFldsLB(Fld,size,borders,borders_new,axis,procind)
		real(psn), dimension(:,:,:), allocatable :: Fld
		real(psn), dimension(:,:,:), allocatable :: Fld_new
		integer :: procind, axis, size
		integer, dimension(0:size) :: borders, borders_new 
		integer :: p1, p2
		integer :: sx, sy, sz, snew ! new size
		integer, dimension(6) :: ind_new, ind_old ! range of indices to copy form the old arrays
		
		!size of the new array 
		sx = mx; sy=my; sz=mz;
		snew = borders_new(procind+1) - borders_new(procind) + 5 
		if(axis.eq.0) sx = snew
		if(axis.eq.1) sy = snew
		if(axis.eq.2) sz = snew
		
		ind_old(1)=1; ind_old(2)=mx; ind_old(3)=1; ind_old(4)=my; ind_old(5)=1; ind_old(6)=mz;
		ind_new = ind_old
		
		p1 = max(borders_new(procind),borders(procind+1)) - borders(procind) +3 
		p2 = min(borders_new(procind+1),borders(procind)) - borders(procind+1) +3 
		
		ind_old(1+2*axis) = p1
		ind_old(2+2*axis) = p2 
		
		p1 = max(borders(procind),borders_new(procind+1)) - borders_new(procind) +3 
		p2 = min(borders(procind+1),borders_new(procind)) - borders_new(procind+1) +3 
		
		ind_new(1+2*axis) = p1
		ind_new(2+2*axis) = p2 
		
		allocate(Fld_new(sx,sy,sz))
		if(p2.ge.p1) Fld_new(ind_new(1):ind_new(2),ind_new(3):ind_new(4),ind_new(5):ind_new(6)) = Fld(ind_old(1):ind_old(2),ind_old(3):ind_old(4),ind_old(5):ind_old(6))
		call move_alloc(Fld_new,Fld)
		
		print*,'new size',sx,sy,sz
	end subroutine ResizeFldsLB
	
	
	
	
	!-------------------------------------------------------------------------
	! Scan and load outgoing particles in the buffer memory
	!-------------------------------------------------------------------------
	
	subroutine allocate_comm_prtl_ngbrs_LB(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m 
	 
		 do n=1, size(list)
			 do m=1,list(n)%num_ngbrs
				 allocate(list(n)%ngbr(m)%p(list(n)%ngbr(m)%pcount))
			 end do 
		 end do 
	 
	end subroutine allocate_comm_prtl_ngbrs_LB
	
	subroutine ScanPrtlOutlierLB(axis,list,edges)
		integer :: axis
		integer, dimension(:) :: edges
		type(ngbr_list) :: list
		real(psn) :: pos
		integer :: n, m, num_ngbrs
		
		num_ngbrs = list%num_ngbrs
			
		do n=1,used_prtl_arr_size
		
			if(flvp(n).eq.0) cycle
		
			if(axis.eq.0) pos = xp(n)
			if(axis.eq.1) pos = yp(n)
			if(axis.eq.2) pos = zp(n)
			
			call SendPrtlNgbrInd(num_ngbrs,edges,pos,m)
			
			list%ngbr(m)%pcount = list%ngbr(m)%pcount +1 
					
		end do 
		
	end subroutine ScanPrtlOutlierLB
	
	
	subroutine LoadPrtlOutlierLB(axis,list,edges)
		integer :: axis
		integer, dimension(:) :: edges
		type(ngbr_list) :: list
		real(psn) :: pos
		integer ::  n, m
			
		do n=1,used_prtl_arr_size
		
			if(flvp(n).eq.0) cycle
						
            call CopyTransferPrtl(qp,xp,yp,zp,up,vp,wp,var1p,flvp,tagp,list,n)
            call RemovePrtl(qp,flvp,tagp,xp,yp,zp,up,vp,wp,n) 
	
		end do 
		
		
		used_prtl_arr_size = 0 
		np = 0
		
	end subroutine LoadPrtlOutlierLB
	

end module loadbalance
