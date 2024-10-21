module loadbalance_procgrid
	use parameters 
	use vars
	use communication
	use subdomains
	use loadbalance_indep
	implicit none
	integer, dimension(4) :: positive_proc_grid_ind = (/ -huge(1), huge(1), -huge(1), huge(1) /)
	
contains 
	
	!The implementation below assumes that indepLBaxis = 0  (x-axis)

    subroutine BalanceLoadProcGrid
		logical :: balance_load

		balance_load = .false.
		
		!periodic load balancing
		if(modulo(t,load_balance_yz_period).eq.0) balance_load = .true.

		if(expand_positive_proc_grid_left(3, positive_proc_grid_ind(1), nSubDomainsY, yborders)) balance_load = .true. 
		if(expand_positive_proc_grid_right(4, positive_proc_grid_ind(2), nSubDomainsY, yborders)) balance_load = .true. 

		if(expand_positive_proc_grid_left(5, positive_proc_grid_ind(3), nSubDomainsZ, zborders)) balance_load = .true. 
		if(expand_positive_proc_grid_right(6, positive_proc_grid_ind(4), nSubDomainsZ, zborders)) balance_load = .true. 

		
		if(balance_load) call AdjustProcGrid

		if(t.lt.100) then
			bc_face(3)%speed = 0.1
			bc_face(4)%speed = 0.0
		end if 

	end subroutine BalanceLoadProcGrid 
	
	!adjust the number of proc. across different locations in ProcGrid that are placed along indepLBaxis
	subroutine AdjustProcGrid
		real(dbpsn) :: load_this, load_tot ! load along indepLBaxis of this proc.
		real(dbpsn), dimension(:,:), allocatable :: load_indepLBaxis
		integer, dimension(:,:), allocatable :: count_new
		real,    dimension(:,:), allocatable :: count_ideal
		integer :: i,j, m,n,k, diff, ntransfer, naccept, ind_proc
		integer, dimension(:), allocatable :: borders_new, proc_new 
		integer :: count
		integer, dimension(:), allocatable :: list_free_proc
		integer :: nfree_proc, free_proc, x1, x2

		!limit of the entire domain, must be same (by design) on every proc.  
		x1 = xborders(0)
		x2 = xborders(nSubDomainsX)
		
		!allocate local variables
		allocate(load_indepLBaxis(0:isizeProcGrid-1,0:jsizeProcGrid-1),count_new(0:isizeProcGrid-1,0:jsizeProcGrid-1), count_ideal(0:isizeProcGrid-1,0:jsizeProcGrid-1)) 
		
		!STEP-0: calc. load along different columns 
		!total load on individual group of proc. working along each indpepLBaxis 
		
		load_indepLBaxis = 0
		
		call load_this_proc(load_this)
		load_indepLBaxis(iproc,jproc) = load_this
		
		call MPI_ALLREDUCE(MPI_IN_PLACE,load_indepLBaxis(0:isizeProcGrid-1,0:jsizeProcGrid-1), isizeProcGrid*jsizeProcGrid , MPI_DOUBLE_PRECISION, mpi_sum, MPI_COMM_WORLD)
		
		!in the case of moving boundaries
		call load_boundaries(load_indepLBaxis)

		!calculate total load and also make a copy proc. count for further update  
		load_tot = 0
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				count_new(i,j)=ProcGrid(i,j)%count
				load_tot = load_tot + load_indepLBaxis(i,j)
			end do 
		end do
		
		!STEP-1 : calculate ideal number of proc. for each column given the current load
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				count_ideal(i,j)= nproc * load_indepLBaxis(i,j) / load_tot
			end do 
		end do
		print*,'count _old', count_new

		        
		!STEP-2: even out the proc. count
		allocate(list_free_proc(nproc))
		nfree_proc = 0 
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1

				diff = int( count_new(i,j) - count_ideal(i,j) )
				ind_proc =  ProcGrid(i,j)%count-1

				if(diff.gt.0) then 

					do m=0,isizeProcGrid-1
						do n=0,jsizeProcGrid-1
                            
							naccept = ceiling( count_ideal(m,n) - count_new(m,n) )
							
							if( naccept .gt.0 ) then 	
								
								ntransfer = min(diff, naccept)
								count_new(i,j) = count_new(i,j) - ntransfer
								count_new(m,n) = count_new(m,n) + ntransfer
								diff = diff - ntransfer
								
								!append the free proc list 
								do k=1,ntransfer
									nfree_proc = nfree_proc+1
									list_free_proc(nfree_proc) =  ProcGrid(i,j)%procs(ind_proc)
									ind_proc = ind_proc -1 
								end do 						
							
							end if 

							if(diff.eq.0) exit 
			
						end do 
					end do 

				end if 

				if(diff.gt.0) count_new(i,j) = count_new(i,j) + diff
			
			end do 
		end do 

		!print*,'count ideal',count_ideal , t, proc
		!print*,'count_new ', count_new, t, proc



		!Now check if a thresold (<2) number of proc. need to be replaced. If not, then return
		count = 0
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				count = max(count, abs(count_new(i,j) - ProcGrid(i,j)%count) ) 
			end do 
		end do 
		if (count.lt.2) then 
			deallocate(load_indepLBaxis,count_new,count_ideal)
			deallocate(list_free_proc)
			return
		end if

		!STEP-3 : reduce proc. count on significantly underloaded ProcGrid locations
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				
				if(count_new(i,j).lt.ProcGrid(i,j)%count) then 					
	                
					allocate(proc_new(0:count_new(i,j)-1))
					
					!prepare the list of new proc. that will work in a group along indepLBaxis 
					if(count_new(i,j).gt.0) proc_new(0:count_new(i,j)-1) = ProcGrid(i,j)%procs(0:count_new(i,j)-1)
					
					!set new borders with reduced number of proc.
					if(i.eq.iproc .and. j.eq.jproc) then
				
						allocate(borders_new(0:count_new(i,j)))
						
						!new borders along indepLBaxis
						call determine_new_borders(ProcGrid(i,j)%count,ProcGrid(i,j)%borders,count_new(i,j),borders_new,x1,x2)
						
						!send-recv data to set new borders
						call SetNewBorders(ProcGrid(i,j)%count,ProcGrid(i,j)%borders,ProcGrid(i,j)%procs,count_new(i,j),borders_new,proc_new)
						
						deallocate(ProcGrid(i,j)%borders)
						allocate(ProcGrid(i,j)%borders(0:count_new(i,j)))
						ProcGrid(i,j)%borders(0:count_new(i,j)) = borders_new(0:count_new(i,j))
						deallocate(borders_new)
						
						!update borders and proc. count along the x-axis
						nSubDomainsX = count_new(i,j)
						deallocate(xborders)
						allocate(xborders(0:nSubDomainsX))
						xborders = ProcGrid(iproc,jproc)%borders
						
					end if
					
					!update proc list and count
					deallocate(ProcGrid(i,j)%procs)
					allocate(ProcGrid(i,j)%procs(0:count_new(i,j)-1))

					
					if(count_new(i,j).gt.0) ProcGrid(i,j)%procs(0:count_new(i,j)-1) = proc_new(0:count_new(i,j)-1)
					
					ProcGrid(i,j)%count = count_new(i,j)
					deallocate(proc_new)
					
				end if
			
			end do 
		end do   
	
	    call MPI_BARRIER(MPI_COMM_WORLD)
		
		
		!STEP-4 : place freed proc. in new locations in ProcGrid
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				
				if(count_new(i,j).gt.ProcGrid(i,j)%count) then 
						
					allocate(proc_new(0:count_new(i,j)-1))
					!prepare the list of new proc. that will work in a group along indepLBaxis 
					if(ProcGrid(i,j)%count-1 .gt.0) proc_new(0:ProcGrid(i,j)%count-1) = ProcGrid(i,j)%procs(0:ProcGrid(i,j)%count-1)
					
					do n = 1,  count_new(i,j) - ProcGrid(i,j)%count
						free_proc =   list_free_proc(nfree_proc)
						nfree_proc = nfree_proc -1 

						proc_new(ProcGrid(i,j)%count-1 +n) = free_proc
						
						!if this proc. is being replaced, set ProcGrid cord. for this proc 
						if(proc.eq.free_proc) then
							deallocate(ProcGrid(iproc,jproc)%borders) !deallocate borders from older proc grid location
							!update position in proc grid
							iproc = i
							jproc = j
							kproc = ProcGrid(i,j)%count-1 +n
							procxind = kproc
							procyind = iproc
							proczind = jproc
							
							!get the existing borders at the new proc grid location
							if(allocated(ProcGrid(i,j)%borders)) deallocate(ProcGrid(i,j)%borders) 
							allocate(ProcGrid(i,j)%borders(0:ProcGrid(i,j)%count))
							if(ProcGrid(i,j)%count.gt.0) then
								call MPI_RECV(ProcGrid(i,j)%borders,ProcGrid(i,j)%count+1,MPI_INTEGER,ProcGrid(i,j)%procs(0),1,MPI_COMM_WORLD, MPI_STATUS_IGNORE)
							end if
								
						end if 

						!send the existing border to new proc. joining this column 
						if(iproc.eq.i.and.jproc.eq.j.and.kproc.eq.0) then
							call MPI_SEND(ProcGrid(i,j)%borders,ProcGrid(i,j)%count+1,MPI_INTEGER,free_proc,1,MPI_COMM_WORLD) 
						endif

					end do
					
					!set new borders with reduced number of proc.
					if(i.eq.iproc .and. j.eq.jproc) then
				
						allocate(borders_new(0:count_new(i,j)))
						
						!new borders along indepLBaxis
						call determine_new_borders(ProcGrid(i,j)%count,ProcGrid(i,j)%borders,count_new(i,j),borders_new,x1,x2)
						
						!send-recv data to set new borders
						print*,'setting new borders started', proc , 'iproc',iproc

						call SetNewBorders(ProcGrid(i,j)%count,ProcGrid(i,j)%borders,ProcGrid(i,j)%procs,count_new(i,j),borders_new,proc_new)
						print*,'setting new borders done', proc

						
						deallocate(ProcGrid(i,j)%borders)
						allocate(ProcGrid(i,j)%borders(0:count_new(i,j)))
						ProcGrid(i,j)%borders(0:count_new(i,j)) = borders_new(0:count_new(i,j))
						deallocate(borders_new)
						
						!update borders and proc. count along the x-axis
						nSubDomainsX = count_new(i,j)
						deallocate(xborders)
						allocate(xborders(0:nSubDomainsX))
						xborders = ProcGrid(iproc,jproc)%borders
						
					
					end if
					
						
					!update proc list and count
					deallocate(ProcGrid(i,j)%procs)
					allocate(ProcGrid(i,j)%procs(0:count_new(i,j)-1))
					
					if(count_new(i,j).gt.0) ProcGrid(i,j)%procs(0:count_new(i,j)-1) = proc_new(0:count_new(i,j)-1)
					ProcGrid(i,j)%count = count_new(i,j)
					
					deallocate(proc_new)
					
				
				end if
			
			end do 
		end do
		if(proc.eq.1) print*,'proc 1 iproc is',iproc

		print*,'I am here',proc
		
		call SetSubDomainGridSize
		call SetPrtlBoundaries
		call InitAuxFld
	
	    call MPI_BARRIER(MPI_COMM_WORLD)

		call ResetCommSubdomainsLB

		call update_positive_proc_grid_ind
				
		deallocate(load_indepLBaxis,count_new,count_ideal)
		deallocate(list_free_proc)

		print*,'I am finished here',proc

	
	end subroutine AdjustProcGrid
	
	subroutine load_this_proc(tot)
		real(dbpsn) :: tot
		real(psn) :: my1, mz1
		
		my1 = min(real(yborders(procyind+1),dbpsn), bc_face(4)%pos_fld + bc_face(4)%attn_thickness) - max(real(yborders(procyind),dbpsn), bc_face(3)%pos_fld - bc_face(3)%attn_thickness ) 
		my1 = max(0.0_dbpsn,my1)

		mz1 = min(real(zborders(proczind+1),dbpsn), bc_face(6)%pos_fld + bc_face(6)%attn_thickness) - max(real(zborders(proczind),dbpsn), bc_face(5)%pos_fld - bc_face(5)%attn_thickness ) 
		mz1 = max(0.0_dbpsn,mz1)
		
		tot = np + (mx-5)*my1*mz1
		 
	end subroutine load_this_proc

	subroutine load_boundaries(load)
		real(dbpsn), dimension(0:isizeProcGrid-1,0:jsizeProcGrid-1) :: load
		real(dbpsn) :: tot

		tot = my*mz*(xborders(nSubDomainsX) - xborders(0))

		if(expand_positive_proc_grid_left(3, positive_proc_grid_ind(1), nSubDomainsY, yborders)) load(positive_proc_grid_ind(1)-1,:) = tot 
		if(expand_positive_proc_grid_right(4, positive_proc_grid_ind(2), nSubDomainsY, yborders)) load(positive_proc_grid_ind(2)+1,:) = tot

		if(expand_positive_proc_grid_left(5, positive_proc_grid_ind(3), nSubDomainsZ, zborders)) load(:,positive_proc_grid_ind(3)-1) = tot
		if(expand_positive_proc_grid_right(6, positive_proc_grid_ind(4), nSubDomainsZ, zborders)) load(:,positive_proc_grid_ind(4)-1) = tot

	end subroutine load_boundaries

	subroutine update_positive_proc_grid_ind
		integer :: i,j
		positive_proc_grid_ind(1) = isizeProcGrid -1
		positive_proc_grid_ind(2) = 0
		positive_proc_grid_ind(3) = jsizeProcGrid -1 
		positive_proc_grid_ind(4) = 0
		do i=0,isizeProcGrid-1
			do j=0,jsizeProcGrid-1
				if(ProcGrid(i,j)%count.gt.0) then 
					positive_proc_grid_ind(1) = min(positive_proc_grid_ind(1),i)
					positive_proc_grid_ind(2) = max(positive_proc_grid_ind(2),i)
					positive_proc_grid_ind(3) = min(positive_proc_grid_ind(3),j)
					positive_proc_grid_ind(4) = max(positive_proc_grid_ind(4),j)
				end if 
			end do 
		end do 
	end subroutine update_positive_proc_grid_ind

	!check if the moving boundaries may go into regions where no proc. are placed
	logical function expand_positive_proc_grid_left(face, proc_grid_ind, nSubDomains, borders)
		integer :: face, proc_grid_ind, nSubDomains
		integer, dimension(0:nSubDomains) :: borders

		expand_positive_proc_grid_left = .false.
		if( bc_face(face)%speed .lt.0 .and. proc_grid_ind .gt. 0) then
			 if( bc_face(face)%pos_fld - bc_face(face)%attn_thickness .lt. borders(proc_grid_ind) +1 ) expand_positive_proc_grid_left = .true.
		end if 
	end function expand_positive_proc_grid_left

	logical function expand_positive_proc_grid_right(face, proc_grid_ind, nSubDomains, borders)
		integer :: face, proc_grid_ind, nSubDomains
		integer, dimension(0:nSubDomains) :: borders

		expand_positive_proc_grid_right = .false.
		if( bc_face(face)%speed .gt.0 .and. proc_grid_ind .lt. nSubDomains-1) then
			if( bc_face(face)%pos_fld + bc_face(face)%attn_thickness .gt. borders(proc_grid_ind+1) -1 ) expand_positive_proc_grid_right = .true.
	    end if 
    end function expand_positive_proc_grid_right



	subroutine determine_new_borders(size,borders,size_new,borders_new,x1,x2)
		integer :: size, size_new
		integer, dimension(0:size)     :: borders
		integer, dimension(0:size_new) :: borders_new
		integer :: n, ncells_tot, ncells, x1, x2 
		
		ncells_tot = x2-x1  
		ncells = ncells_tot/size_new 
		
		borders_new(0) = x1
		do n=1,size_new
			borders_new(n) = borders_new(n-1)+ncells
            borders_new(n) = min(borders_new(n), x2-4*(size_new-n) ) !make sure that there are enough cells (4 each at least) for the rest
		end do
		borders_new(size_new) = x2
		
	end subroutine determine_new_borders

	
end module loadbalance_procgrid