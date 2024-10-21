module subdomains
	use parameters
	use vars
	use communication
	use comm_fld
	use comm_prtl
#ifdef cyl
    use cyl_common
#endif 	
	implicit none

	integer :: num_send_indep_borders = 0
	type(MPI_Request) :: mpi_req_recv_indep_borders! = MPI_REQUEST_NULL
	type(MPI_Request), dimension(:), allocatable :: mpi_req_send_indep_borders

contains

	subroutine InitDomainSkelton
		integer :: n,m
		 
		 mpi_req_recv_indep_borders = MPI_REQUEST_NULL
		 
		 call InitProcGridPlane
		 call DetermineProcGridCount
		 call InitProcGrid
		 call SetNgbr
		 
		 call SetProcGridCommGroup
		
	 end subroutine InitDomainSkelton
	 
	 subroutine RestartDomainSkelton
		 mpi_req_recv_indep_borders = MPI_REQUEST_NULL
		 call SetNgbr
		 call SetProcGridCommGroup
	 end subroutine RestartDomainSkelton
     
  	!---------------------------------------------------------------------------------------
  	! Local limits and size of the grid
  	!--------------------------------------------------------------------------------------- 
	
     subroutine SetSubDomainGridSize	 
          mx=xborders(procxind+1)-xborders(procxind)+5
          my=yborders(procyind+1)-yborders(procyind)+5
#ifdef twoD
          mz=1 
#else
          mz=zborders(proczind+1)-zborders(proczind)+5
#endif     
     end subroutine SetSubDomainGridSize
	 
     subroutine SetPrtlBoundaries
          xmin=3.0_psn
          xmax=mx-2
          ymin=3.0_psn
          ymax=my-2
#ifdef twoD
          zmin=1.0_psn
          zmax=2.0_psn
#else          
          zmin=3.0_psn
          zmax=mz-2
#endif
          xlen=xmax-xmin
          ylen=ymax-ymin
          zlen=zmax-zmin
     end subroutine SetPrtlBoundaries
	 

 	!---------------------------------------------------------------------------------------
 	! Layout of the ProcGrid matrix 
 	!--------------------------------------------------------------------------------------- 
	 subroutine InitProcGridPlane	
#ifdef twoD
         nSubDomainsZ = 1 
#endif		 	 
		 select case (indepLBaxis)
		 	case(0)
				isizeProcGrid = nSubDomainsY
				jsizeProcGrid = nSubDomainsZ
			case(1)
				isizeProcGrid = nSubDomainsX
				jsizeProcGrid = nSubDomainsZ
			case(2)
				isizeProcGrid = nSubDomainsX
				jsizeProcGrid = nSubDomainsY
	     end select
		 		 
		 allocate(ProcGrid(0:isizeProcGrid-1,0:jsizeProcGrid-1)) 
	 end subroutine InitProcGridPlane

	 
	!---------------------------------------------------------------------------------------
	! 0 dirn of ProcGrid is along indepLBaxis and other indices are along the other two, as defined below
	!--------------------------------------------------------------------------------------- 
 	integer function axisProcGrid_To_CordAxis(axisProcGrid)
 		integer :: axisProcGrid
 		if(axisProcGrid.eq.0) axisProcGrid_To_CordAxis = indepLBaxis
 		if(axisProcGrid.eq.1) then 
 			if(indepLBaxis.eq.0) axisProcGrid_To_CordAxis = 1 
 			if(indepLBaxis.ge.1) axisProcGrid_To_CordAxis = 0 
 		end if
 		if(axisProcGrid.eq.2) then 
 			if(indepLBaxis.le.1) axisProcGrid_To_CordAxis = 2
 			if(indepLBaxis.eq.2) axisProcGrid_To_CordAxis = 1
 		end if
 	end function axisProcGrid_To_CordAxis
	 
	!---------------------------------------------------------------------------------------
	! Initial determination of the number of MPI tasks along indepLBaxis, different locations in the 2D ProcGrid  
	!---------------------------------------------------------------------------------------
	 
	 subroutine DetermineProcGridCount
		 integer :: i,j
		 integer :: nproc_indepLBaxis! number of proc. along each indepLBaxis
		 integer :: count
		 
		 
		 nproc_indepLBaxis = nproc / (isizeProcGrid*jsizeProcGrid)
		 if(nproc_indepLBaxis.eq.0) then 
		 	print*,'Error :: looks like there are not enough mpi tasks to distribute min. no. of subdomains'
			STOP 
		 end if 
		 
		 !by default, first evenly distribute proc. along every indepLB axis
		 count = 0
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				 ProcGrid(i,j)%count = nproc_indepLBaxis
				 count = count + nproc_indepLBaxis
			 end do 
		 end do 
		 
		 count = nproc - count ! this many proc. are still left to be placed
		 do while(count.gt.0)
			
			 do i=0,isizeProcGrid-1
				 do j=0,jsizeProcGrid-1
					 
					 if(count.gt.0) then 
						 ProcGrid(i,j)%count = ProcGrid(i,j)%count +1
						 count = count - 1
					 end if

				 end do 
			 end do 
			
		 end do
		  
	 end subroutine DetermineProcGridCount
	 
	 subroutine InitProcGrid
		 integer :: i,j,k, rank
		 
		 rank=0
		 do i=0,isizeProcGrid-1
			 do j=0,jsizeProcGrid-1
				 allocate(ProcGrid(i,j)%procs(0:ProcGrid(i,j)%count-1))
				 
				 do k=0,ProcGrid(i,j)%count-1
					 ProcGrid(i,j)%procs(k)=rank
					 if(rank.eq.proc) then 
						 iproc = i; jproc = j; kproc = k; 
					 end if 
					 rank = rank+1 
				 end do 
			 end do
		 end do
		 
		 !default initialisations of the borders, borders along the indepLBaxis is then updated later
		 call InitBorders1D(nSubDomainsX,xborders,box_bounds(1),box_bounds(2))
		 call InitBorders1D(nSubDomainsY,yborders,box_bounds(3),box_bounds(4))
		 call InitBorders1D(nSubDomainsZ,zborders,box_bounds(5),box_bounds(6))
		 
		 !set the borders along the indepLBaxis; overrides the borders defined above
		 call SetProcIndepBorders
	 
	 end subroutine InitProcGrid
	 
	 
	 subroutine SetProcIndepBorders
		 select case (indepLBaxis)
	     	case(0)
				procxind = kproc; procyind= iproc; proczind=jproc;
				nSubDomainsX=ProcGrid(iproc,jproc)%count
				call InitBorders1D(nSubDomainsX,xborders,box_bounds(1),box_bounds(2))
				call CopyBorders(nSubDomainsX,xborders,ProcGrid(iproc,jproc)%borders)
		 	case(1)
				procxind = iproc; procyind= kproc; proczind=jproc;
				nSubDomainsY=ProcGrid(iproc,jproc)%count
				call InitBorders1D(nSubDomainsY,yborders,box_bounds(3),box_bounds(4))
				call CopyBorders(nSubDomainsY,yborders,ProcGrid(iproc,jproc)%borders)
		 	case(2)
				procxind = iproc; procyind= jproc; proczind=kproc;
				nSubDomainsZ=ProcGrid(iproc,jproc)%count
				call InitBorders1D(nSubDomainsZ,zborders,box_bounds(5),box_bounds(6))
				call CopyBorders(nSubDomainsZ,zborders,ProcGrid(iproc,jproc)%borders)
		 end select
	 end subroutine SetProcIndepBorders
	 
	 
	 subroutine InitBorders1D(nSubDomains,borders,b1,b2)
		 integer, dimension(:), allocatable :: borders
		 integer :: nSubDomains, b1, b2
		 integer :: i 
		 if(allocated(borders)) deallocate(borders)
		 allocate(borders(0:nSubDomains))
		 
		 borders(0) = b1
		 do i=0,nSubDomains
			 borders(i)	= b1 + i*( (b2-b1)/nSubDomains )
		 end do 
		 borders(nSubDomains) = b2 
		 
	 end subroutine InitBorders1D
	 
	 subroutine CopyBorders(nSubDomains,source,dest)
		 integer :: nSubDomains
		 integer, dimension(:), allocatable :: source, dest
		 if(allocated(dest)) deallocate(dest)
		 allocate(dest(0:nSubDomains))
		 dest = source 
	 end subroutine CopyBorders
	 
  	!---------------------------------------------------------------------------------------
  	! Send/recv the updated indepLBaxis proc list to/from any proc. column in the 2D ProcGrid 
  	!---------------------------------------------------------------------------------------
	 
 	 subroutine SendIndepBordersTo(isend,jsend)
 		 integer :: isend, jsend !i,j indices of the column where mpi proc. will receive updated borders 
 		 integer :: n, count, ksend

 		 if(isend.eq.iproc.and.jsend.eq.jproc) return

		 count = ProcGrid(iproc,jproc)%count
		 
		 !first count how many sends will be posted
		 ksend = kproc  
		 num_send_indep_borders = 0 
		 do while ( ksend+1 .le. ProcGrid(isend,jsend)%count )
			if(count.ge.0) num_send_indep_borders = num_send_indep_borders +1
			ksend = ksend + ProcGrid(iproc,jproc)%count
		 end do

		 !allocate and initialise 'mpi_req' for sending the borders 
		 if(allocated(mpi_req_send_indep_borders)) deallocate(mpi_req_send_indep_borders)
		 allocate(mpi_req_send_indep_borders(num_send_indep_borders))
		 do n=1,num_send_indep_borders
			mpi_req_send_indep_borders(n) = MPI_REQUEST_NULL
		 end do 
		 
		 
		 !now send the borders
		 ksend = kproc  
		 n=0
		 do while ( ksend+1 .le. ProcGrid(isend,jsend)%count )
 			 
			 if(count.gt.0) then 
				n=n+1
				call MPI_ISEND(ProcGrid(iproc,jproc)%borders,count+1,MPI_INTEGER,ProcGrid(isend,jsend)%procs(ksend),1,MPI_COMM_WORLD, mpi_req_send_indep_borders(n))	
			 end if 		 
 			 ksend = ksend + ProcGrid(iproc,jproc)%count
 		 
		 end do 

 	 end subroutine SendIndepBordersTo
	 
 	 subroutine RecvIndepBordersFrom(irecv,jrecv)
 		 integer :: irecv, jrecv !i,j indices of the column where mpi proc. will receive updated borders 
 		 integer :: n,count, comm_proc

 		 if(irecv.eq.iproc.and.jrecv.eq.jproc) return
		 
 		 count = ProcGrid(irecv,jrecv)%count
		 if(allocated(ProcGrid(irecv,jrecv)%borders)) deallocate(ProcGrid(irecv,jrecv)%borders)
 		 allocate(ProcGrid(irecv,jrecv)%borders(0:count))

		 mpi_req_recv_indep_borders = MPI_REQUEST_NULL
 		 if(count.gt.0) then 
			comm_proc = ProcGrid(irecv,jrecv)%procs(mod(kproc,count))
			call MPI_IRECV(ProcGrid(irecv,jrecv)%borders,count+1,MPI_INTEGER,comm_proc,1,MPI_COMM_WORLD, mpi_req_recv_indep_borders)
		 end if 
	 
 	 end subroutine RecvIndepBordersFrom

	 subroutine WaitSendRecvIndepBorders
		integer :: n

		!wait for all sends to complete
		do n = 1, num_send_indep_borders
			call MPI_WAIT(mpi_req_send_indep_borders(n),MPI_STATUS_IGNORE)
		end do 

		!wait for the recv to complete
		call MPI_WAIT(mpi_req_recv_indep_borders,MPI_STATUS_IGNORE)
	 
	end subroutine WaitSendRecvIndepBorders
	
	 !---------------------------------------------------------------------------------------
	 ! Fill in the list of ngbrs for comm.
	 !---------------------------------------------------------------------------------------
	 subroutine SetNgbr
		 		 
		 if(allocated(ngbr_send)) then 
			 call DeallocateNgbrList(ngbr_send)
			 deallocate(ngbr_send)
		 end if 
		 if(allocated(ngbr_recv)) then 
			 call DeallocateNgbrList(ngbr_recv) 
			 deallocate(ngbr_recv)
		 end if 

#ifndef twoD 		 
		 allocate(ngbr_send(6),ngbr_recv(6)) ! 6 faces
		 		 
		 if(indepLBaxis.eq.0) call SetMultipleNgbr(iproc,jproc,(/ 3, 4, 5, 6 /))
		 if(indepLBaxis.eq.1) call SetMultipleNgbr(iproc,jproc,(/ 1, 2, 5, 6 /))
		 if(indepLBaxis.eq.2) call SetMultipleNgbr(iproc,jproc,(/ 1, 2, 3, 4 /))
		 
		 call SingleNgbr(ProcGrid(iproc,jproc)%count,ProcGrid(iproc,jproc)%borders,ProcGrid(iproc,jproc)%procs,indepLBaxis,ngbr_send,ngbr_recv)

#else

		 allocate(ngbr_send(4),ngbr_recv(4)) ! 4 faces
		 
		 if(indepLBaxis.eq.0) call SetMultipleNgbr(iproc,jproc,(/ 3, 4, 0 ,0 /))
		 if(indepLBaxis.eq.1) call SetMultipleNgbr(iproc,jproc,(/ 1, 2, 0 ,0 /))
		 if(indepLBaxis.eq.2) then 
			 print*, ' indepLBaxis = 2 ! should be 0 (x-axis) or 1 (y-axis) in 2D simulations '
		 end if 
		 
		 call SingleNgbr(ProcGrid(iproc,jproc)%count,ProcGrid(iproc,jproc)%borders,ProcGrid(iproc,jproc)%procs,indepLBaxis,ngbr_send,ngbr_recv)
#endif 

		 
		 !allocate prtl buffer meomory
		 call allocate_comm_prtl_ngbrs(ngbr_send)
		 call allocate_comm_prtl_ngbrs(ngbr_recv)
		 
	 end subroutine SetNgbr
	 
	 !---------------------------------------------------------------------------------------
	 ! mutiple ngbrs and corresponding edges not along indepLB axis
	 !---------------------------------------------------------------------------------------
	 
	 subroutine SetMultipleNgbr( i0,j0, face )
		 integer :: i0,j0
		 integer, dimension(4) :: face	 
		 integer :: i1,j1, i2, j2
		 
		 i1 = intPrdcDomain(i0-1,isizeProcGrid); j1 =j0;
		 i2 = intPrdcDomain(i0+1,isizeProcGrid); j2 =j0;	
		 
		 call SendIndepBordersTo(i1,j0)
		 call RecvIndepBordersFrom(i2,j0)
		 call WaitSendRecvIndepBorders	

		 call SendIndepBordersTo(i2,j0)
		 call RecvIndepBordersFrom(i1,j0)
		 call WaitSendRecvIndepBorders	

		 		 
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i1,j1)%count, ProcGrid(i1,j1)%borders, ProcGrid(i1,j1)%procs, ngbr_send(face(1)) )
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i1,j1)%count, ProcGrid(i1,j1)%borders, ProcGrid(i1,j1)%procs, ngbr_recv(face(1)) )
		 
	     call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i2,j2)%count, ProcGrid(i2,j2)%borders, ProcGrid(i2,j2)%procs, ngbr_send(face(2)) )
	     call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i2,j2)%count, ProcGrid(i2,j2)%borders, ProcGrid(i2,j2)%procs, ngbr_recv(face(2)) )
		 
#ifndef twoD
			 
		 i1 = i0; j1 =intPrdcDomain(j0-1,jsizeProcGrid);
		 i2 = i0; j2 =intPrdcDomain(j0+1,jsizeProcGrid);

		 call SendIndepBordersTo(i0,j1)
		 call RecvIndepBordersFrom(i0,j2)
		 call WaitSendRecvIndepBorders	

		 call SendIndepBordersTo(i0,j2)
		 call RecvIndepBordersFrom(i0,j1)
		 call WaitSendRecvIndepBorders
		 
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i1,j1)%count, ProcGrid(i1,j1)%borders, ProcGrid(i1,j1)%procs, ngbr_send(face(3)) )
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i1,j1)%count, ProcGrid(i1,j1)%borders, ProcGrid(i1,j1)%procs, ngbr_recv(face(3)) )
		 
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i2,j2)%count, ProcGrid(i2,j2)%borders, ProcGrid(i2,j2)%procs, ngbr_send(face(4)) )
		 call MultipleNgbr(ProcGrid(i0,j0)%count, ProcGrid(i0,j0)%borders, ProcGrid(i2,j2)%count, ProcGrid(i2,j2)%borders, ProcGrid(i2,j2)%procs, ngbr_recv(face(4)) )

#endif		 
		  	 
	 end subroutine SetMultipleNgbr
	 
	 !---------------------------------------------------------------------------------------
	 ! set mutiple ngbrs and corresponding edges from known overlapping "borders"
	 !---------------------------------------------------------------------------------------
	 	 
	 subroutine MultipleNgbr( size, borders, size_ngbr, borders_ngbr, proc_ngbr, list)
		 integer :: size, size_ngbr
		 integer, dimension(0:size) :: borders
	 	 integer, dimension(0:size_ngbr)   :: borders_ngbr
		 integer, dimension(0:size_ngbr-1) :: proc_ngbr 
		 type(ngbr_list) :: list
		 integer :: n , p1, p2, num_ngbrs, m
		 
		 num_ngbrs = 0
		 if(kproc+1.le.size) num_ngbrs = count_ngbrs(size_ngbr,borders_ngbr,borders(kproc),borders(kproc+1))

		 call allocate_ngbr_list(list, num_ngbrs)

		 if(num_ngbrs.eq.0) return

		 m=1
		 do n=0,size_ngbr-1
			 p1 = max(borders_ngbr(n),borders(kproc))
			 p2 = min(borders_ngbr(n+1),borders(kproc+1))
			 if(p2.gt.p1) then
				 
				 list%ngbr(m)%proc = proc_ngbr(n)
				 
				 list%edges(m) = p1 - borders(kproc) +3 ! local cordinate
				 list%edges(m+1) = p2 - borders(kproc) +3 
				 m=m+1
			 end if
		 end do
		 
		 ! include ghost cells as well, these edges are used only on sending particles and setting fld comm domain
		 call extend_edge_range(list%edges)

	 end subroutine MultipleNgbr
	 
	 !---------------------------------------------------------------------------------------
	 ! 2 ngbrs along the indepLB axis and edges
	 !---------------------------------------------------------------------------------------
	 
	 subroutine SingleNgbr(size, borders, procs, axis, send, recv)
		 type(ngbr_list), dimension(:) :: send, recv   
		 integer :: size, axis
		 integer, dimension(0:size)   :: borders
		 integer, dimension(0:size-1) :: procs
		 integer :: n
		 
		 call allocate_ngbr_list(send(1+2*axis),1) 
		 call allocate_ngbr_list(recv(1+2*axis),1) 
		 call allocate_ngbr_list(send(2+2*axis),1) 
		 call allocate_ngbr_list(recv(2+2*axis),1) 
		 
		 do n=0,size-1
			 if(proc.eq.procs(n)) then  
				 send(1+2*axis)%ngbr(1)%proc = procs(intPrdcDomain(n-1,size))
				 recv(1+2*axis)%ngbr(1)%proc = procs(intPrdcDomain(n-1,size))
				 send(2+2*axis)%ngbr(1)%proc = procs(intPrdcDomain(n+1,size))
				 recv(2+2*axis)%ngbr(1)%proc = procs(intPrdcDomain(n+1,size))  
				 
				 send(1+2*axis)%edges(1) = 1
				 send(1+2*axis)%edges(2) = 3
				 recv(1+2*axis)%edges(1) = 1
				 recv(1+2*axis)%edges(2) = 3
				 
				 send(2+2*axis)%edges(1) = borders(kproc+1)    - borders(kproc) +3 
				 send(2+2*axis)%edges(2) = borders(kproc+1) +2 - borders(kproc) +3 
				 recv(2+2*axis)%edges(1) = borders(kproc+1)    - borders(kproc) +3
				 recv(2+2*axis)%edges(2) = borders(kproc+1) +2 - borders(kproc) +3
				  
				 exit
			 end if
		 end do 
		 
	 end subroutine SingleNgbr 
	 
	 integer function intPrdcDomain(i,size)
		 integer :: i, size
		 intPrdcDomain = i 
		 if(i.lt.0) intPrdcDomain=i+size
		 if(i.gt.size-1) intPrdcDomain=i-size
	 end function intPrdcDomain
	 
	 subroutine DeallocateNgbrList(list)
		 type(ngbr_list), dimension(:), allocatable :: list
		 integer :: n,m,i
		 
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 
				 if(allocated(list(n)%ngbr(m)%p)) deallocate(list(n)%ngbr(m)%p)
				 
				 if(allocated(list(n)%ngbr(m)%Fldx)) deallocate(list(n)%ngbr(m)%Fldx)
				 if(allocated(list(n)%ngbr(m)%Fldy)) deallocate(list(n)%ngbr(m)%Fldy)
				 if(allocated(list(n)%ngbr(m)%Fldz)) deallocate(list(n)%ngbr(m)%Fldz)
	
			 end do
			 
			 if(allocated(list(n)%ngbr)) deallocate(list(n)%ngbr)
			 if(allocated(list(n)%edges)) deallocate(list(n)%edges)
		 
		 end do	 
	 
	 end subroutine DeallocateNgbrList
	
	 
	 !---------------------------------------------------------------------------------
	 !  BC is periodic by default; turn off the periodic BC if the BC is set to something else
	 !---------------------------------------------------------------------------------	
	    subroutine SetBoundaryNgbr
			if(bc_face(1)%type_fld.ne.'prdc'.or.bc_face(2)%type_fld.ne.'prdc') call PeriodicX_Off
			if(bc_face(3)%type_fld.ne.'prdc'.or.bc_face(4)%type_fld.ne.'prdc') call PeriodicY_Off
			if(bc_face(5)%type_fld.ne.'prdc'.or.bc_face(6)%type_fld.ne.'prdc') call PeriodicZ_Off
	    end subroutine SetBoundaryNgbr	
		
	 	subroutine PeriodicOff(side)
	 		integer :: side
	 		select case (side)
	 			case(1,2)
	 				call PeriodicX_Off
	 			case(3,4)
	 				call PeriodicY_Off
	 			case(5,6)
	 				call PeriodicZ_Off
	 	    end select
	 	end subroutine PeriodicOff
	 	
		subroutine PeriodicX_Off
		
	 		if(procxind.eq.0) call SetNgbrNull(ngbr_send(1))
	 		if(procxind.eq.0) call SetNgbrNull(ngbr_recv(1))
		
	 		if(procxind.eq.nSubDomainsX-1) call SetNgbrNull(ngbr_send(2))
	 		if(procxind.eq.nSubDomainsX-1) call SetNgbrNull(ngbr_recv(2))
		
	 	end subroutine PeriodicX_Off	
	 	
		subroutine PeriodicY_Off		
	 		if(procyind.eq.0) call SetNgbrNull(ngbr_send(3))
	 		if(procyind.eq.0) call SetNgbrNull(ngbr_recv(3))
		
	 		if(procyind.eq.nSubDomainsY-1) call SetNgbrNull(ngbr_send(4))
	 		if(procyind.eq.nSubDomainsY-1) call SetNgbrNull(ngbr_recv(4))
	 	end subroutine PeriodicY_Off
	 	
		subroutine PeriodicZ_Off
#ifdef twoD
			print*,'Warning :: setting BC along the z-axis in 2D simulation is invalid'
#endif					
	 		if(proczind.eq.0) call SetNgbrNull(ngbr_send(5))
	 		if(proczind.eq.0) call SetNgbrNull(ngbr_recv(5))
		
	 		if(proczind.eq.nSubDomainsZ-1) call SetNgbrNull(ngbr_send(6))
	 		if(proczind.eq.nSubDomainsZ-1) call SetNgbrNull(ngbr_recv(6))
	 	end subroutine PeriodicZ_Off
	
	 	subroutine SetNgbrNull(list)
	 	    type(ngbr_list) :: list
	 		integer :: m
	 		do m = 1, list%num_ngbrs
	 			list%ngbr(m)%proc = MPI_PROC_NULL
	 		end do
	 	end subroutine SetNgbrNull	
	 
	 !---------------------------------------------------------------------------------------------------------
	 ! if the edges are used in binning particle outliers then extend by 2 both sides to include the ghost cells 
	 !---------------------------------------------------------------------------------------------------------
	 subroutine extend_edge_range(edges)
		 integer :: arr_size
		 integer, dimension(:) :: edges
		 arr_size = size(edges)
		 edges(1)= edges(1) -2
		 edges(arr_size) = edges(arr_size) +2
	 end subroutine extend_edge_range
	 
	 integer function count_ngbrs(size_ngbr,borders_ngbr,b1,b2)
		 integer :: size_ngbr, b1, b2
	 	 integer, dimension(0:size_ngbr) :: borders_ngbr
		 integer :: n, p1, p2 
		 count_ngbrs =0
		 do n=0,size_ngbr-1
			 p1 = max(borders_ngbr(n),b1)
			 p2 = min(borders_ngbr(n+1),b2)
			 if(p2.gt.p1) count_ngbrs = count_ngbrs+1 
		 end do	
	 end function count_ngbrs
	 
	 subroutine allocate_ngbr_list(list, size)
		 integer :: size
		 type(ngbr_list) :: list
		 list%num_ngbrs = size	 
		 if(allocated(list%ngbr))  deallocate(list%ngbr)
		 if(allocated(list%edges)) deallocate(list%edges)
		 allocate(list%ngbr(size),list%edges(size+1)) 
	 end subroutine allocate_ngbr_list
	 
	 !---------------------------------------------------------------------------------------------------------
	 ! set indices that define range of boundary cells whose fld data are communicated
	 !---------------------------------------------------------------------------------------------------------
	 subroutine SetSendRecvFldDomain
		 integer :: n, off, nlayers, tag_offset
		 do n=1,size(ngbr_send)
			 off = 2; nlayers = 3; tag_offset = 10000*n;   
			 if(mod(n,2).eq.0) then 
				  off = -5; nlayers = 2; 
			 end if 
			 call SetCommFldDomain(ngbr_send(n)%ngbr,ngbr_send(n)%edges, n, off, nlayers, tag_offset)
			 if(proc.eq.1) print*,'edges', ngbr_send(n)%edges, off	
		 end do 
		 call AllocCommFld(ngbr_send)  
		 
		 do n=1,size(ngbr_recv)
			 off = 0; nlayers = 2; tag_offset = 10000*(n+1) 
			 if(mod(n,2).eq.0) then 
				  off = -3; nlayers = 3; tag_offset = 10000*(n-1)
			 end if 
			 call SetCommFldDomain(ngbr_recv(n)%ngbr,ngbr_recv(n)%edges, n, off, nlayers, tag_offset)			 
		 end do
		 call AllocCommFld(ngbr_recv)  
	 end subroutine SetSendRecvFldDomain
	 

	 
	 subroutine SetCommFldDomain(ngbr, edges, face, off, nlayers, tag_offset)
		 type(neighbor), dimension(:) :: ngbr
		 integer, dimension(:) :: edges	 
		 integer :: face, tag_offset
		 integer :: off, nlayers
		 integer :: n, e1 , e2
		 
		 !step 1 : set all ghost cells to be exchanged
		 do n=1,size(ngbr)
		
			 ngbr(n)%ind(1)=1; ngbr(n)%ind(2)=mx; ngbr(n)%ind(3)=1; ngbr(n)%ind(4)=my; ngbr(n)%ind(5)=1; ngbr(n)%ind(6)=mz;
			 
			 call set_layers( face2axis(face), off, nlayers, ngbr(n)%ind)
			 
			 ngbr(n)%mpi_tag_offset = tag_offset
			 if(proc.eq.1) print*,n,'ind are',ngbr(n)%ind

		 end do
		 
		 if(face2axis(face).eq.indepLBaxis) return ! no mutiple neighbors
		 
		 !step 2; slice for actual neghbrs: set layers for all neighbors along the independent load balancing axis
		 do n=1,size(ngbr)

! 			 e1 = edges(n)
! 			 e2 = edges(n+1)
!
! 			 if(n.ne.1)          e1 = e1 -2
! 			 if(n.ne.size(ngbr)) e2 = e2 +2
			 
			 e1 = edges(n)
			 e2 = edges(n+1) -1
			 ! edges are extended, so add/substract ghost cells from the ends
			 if(n.eq.1) e1 = e1 + 2
			 if(n.eq.size(ngbr)) e2 = e2 -2
			 
			  
			 ngbr(n)%ind (1+2*indepLBaxis) = e1 
			 ngbr(n)%ind (2+2*indepLBaxis) = e2
			 if(proc.eq.1) print*,n,'ind are',ngbr(n)%ind

		 end do 		 
		 
	 end subroutine SetCommFldDomain
	 

	 	 
	 
	 !set the range of indices corresposning to 2D layers in the subdomain
	 subroutine set_layers(axis, off, nlayers, ind)
		 integer :: axis, off, nlayers
		 integer, dimension(6) :: ind
		 integer :: l1,l2, off1
		 
		 l1 = 1 + 2*axis
		 l2 = 2 + 2*axis
		 
		 off1 = off
		 if(axis.eq.0 .and. off.lt.0) off1 = off+mx
		 if(axis.eq.1 .and. off.lt.0) off1 = off+my
		 if(axis.eq.2 .and. off.lt.0) off1 = off+mz
		 
		 ind(l1) = off1 + 1
		 ind(l2) = off1 + 1 + (nlayers-1)		 	  
	 end subroutine set_layers
	 
	 
	 subroutine AllocCommFld(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m
		 do n = 1,size(list)
			 do m = 1,list(n)%num_ngbrs
				 call allocate_buff_flds(list(n)%ngbr(m),list(n)%ngbr(m)%ind)
			 end do 
		 end do 	 
	 end subroutine AllocCommFld
	 
	 subroutine allocate_buff_flds(ngbr,ind)
		 type(neighbor) :: ngbr
		 integer, dimension(:) :: ind
		 integer :: sx,sy,sz	 
		 if(allocated(ngbr%Fldx)) deallocate(ngbr%Fldx,ngbr%Fldy,ngbr%Fldz)
		 sx=ind(2)-ind(1) +1 
		 sy=ind(4)-ind(3) +1
		 sz=ind(6)-ind(5) +1 
		 allocate(ngbr%Fldx(sx,sy,sz))	
		 allocate(ngbr%Fldy(sx,sy,sz))	
		 allocate(ngbr%Fldz(sx,sy,sz))
		 ngbr%mpi_req_fldx = MPI_REQUEST_NULL
		 ngbr%mpi_req_fldy = MPI_REQUEST_NULL
		 ngbr%mpi_req_fldz = MPI_REQUEST_NULL			 
	 end subroutine allocate_buff_flds
	 
	 
	 !---------------------------------------------------------------------------------------------------------
	 ! Prtl cordinates need to be shifted before/after send/recv  
	 !---------------------------------------------------------------------------------------------------------
	 subroutine SetSendRecvPrtl
		 
		 call SetSendRecvPrtlOffsets(ngbr_send)
		 call SetSendRecvPrtlOffsets(ngbr_recv)
		 
	 end subroutine SetSendRecvPrtl
	 
	 subroutine SetSendRecvPrtlOffsets(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: m,n
			 
		 do n=1,size(list)
			 do m=1,list(n)%num_ngbrs
				 list(n)%ngbr(m)%xshift = 0 
				 list(n)%ngbr(m)%yshift = 0 
				 list(n)%ngbr(m)%zshift = 0
				 if(n.eq.2)  list(n)%ngbr(m)%xshift = xlen
				 if(n.eq.4)  list(n)%ngbr(m)%yshift = ylen
				 if(n.eq.6)  list(n)%ngbr(m)%zshift = zlen
			 
				 if(indepLBaxis.eq.0 .and. list(n)%edges(m).gt.xmin) list(n)%ngbr(m)%xshift = list(n)%edges(m) - xmin
				 if(indepLBaxis.eq.1 .and. list(n)%edges(m).gt.ymin) list(n)%ngbr(m)%yshift = list(n)%edges(m) - ymin
				 if(indepLBaxis.eq.2 .and. list(n)%edges(m).gt.zmin) list(n)%ngbr(m)%zshift = list(n)%edges(m) - zmin
				 
			 end do 
		 end do
		 
	 end subroutine SetSendRecvPrtlOffsets
	 
	 !---------------------------------------------------------------------------------------------------------
	 ! face(1-6) to axis(1-3) : e.g., 1(left),2(right)- along x-axis=1; 
	 !---------------------------------------------------------------------------------------------------------
	 
	 integer function face2axis(face)
		 integer :: face
		 face2axis=0
		 select case (face)
	 	 	case(1,2)
				face2axis = 0
			case(3,4)
				face2axis = 1 
			case(5,6)
				face2axis = 2
		 end select
	 end function face2axis
	 
    !--------------------------------------------------------------------------
 	! Create groups for limited MPI comm. along ProcGrid axis 
 	!--------------------------------------------------------------------------
 	subroutine SetProcGridCommGroup
		integer :: color
 		call MPI_Comm_split(MPI_COMM_WORLD, (jsizeProcGrid-1)*jproc+iproc, kproc, comm_indepLBaxis)
		
		color = kproc
		if(kproc.ne.0) color = MPI_UNDEFINED
		call MPI_Comm_split(MPI_COMM_WORLD, color , (jsizeProcGrid-1)*jproc+iproc, comm_kproc0) ! used to communicate among kproc =0 proc. only
 	end subroutine SetProcGridCommGroup
	
	subroutine FreeProcGridCommGroup
		call MPI_comm_free(comm_indepLBaxis)
		if(comm_kproc0.ne.MPI_COMM_NULL) call MPI_comm_free(comm_kproc0)
	end subroutine FreeProcGridCommGroup
		 
end module subdomains