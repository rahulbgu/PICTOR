module communication 
     use parameters
     use vars
	 use mpi_f08 
     implicit none 
     !include 'mpif.h'
	 
	 
	integer :: ierr ! used globally
	type(MPI_DATATYPE) :: mpi_psn
	type(MPI_DATATYPE) :: mpi_prtltype
     
contains 
     subroutine StartMPI
          integer :: extent,nblocks
          type(MPI_DATATYPE), dimension(0:7) :: oldtypes
		  integer , dimension(0:7) :: blockcounts, offsets
          call mpi_init(ierr)
          call mpi_comm_rank(MPI_COMM_WORLD,proc,ierr)
          call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
          
          !  a) set the precision for mpi data transfer b)Create data type for particle transfer 
          ! **Warning**: The following must be modified if the structure of particle datatype is changed  
          nblocks=0
          select case(psn)
          case(kind(1.0d0))
             mpi_psn=MPI_DOUBLE_PRECISION !mpi precision 
			
             !reals
             offsets(nblocks)=0
             oldtypes(nblocks)=MPI_DOUBLE_PRECISION
             blockcounts(nblocks)=10
             call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,extent,ierr)
             nblocks=nblocks+1
			 
             !integers
             offsets(nblocks)=offsets(nblocks-1)+blockcounts(nblocks-1)*extent
             oldtypes(nblocks)=MPI_INTEGER 
             blockcounts(nblocks)=4
             call MPI_TYPE_EXTENT(MPI_INTEGER,extent,ierr)
             nblocks=nblocks+1 

          case(kind(1.0e0))
             mpi_psn=MPI_REAL !mpi precision
           
             !reals
             offsets(nblocks)=0
             oldtypes(nblocks)=MPI_REAL 
             blockcounts(nblocks)=10
             call MPI_TYPE_EXTENT(MPI_REAL,extent,ierr)
             nblocks=nblocks+1       
			 
             !integers
             offsets(nblocks)=offsets(nblocks-1)+blockcounts(nblocks-1)*extent
             oldtypes(nblocks)=MPI_INTEGER 
             blockcounts(nblocks)=4
             call MPI_TYPE_EXTENT(MPI_INTEGER,extent,ierr)
             nblocks=nblocks+1   

          end select           
          
          call MPI_TYPE_STRUCT(nblocks,blockcounts,offsets,oldtypes,mpi_prtltype,ierr)
          call MPI_TYPE_COMMIT(mpi_prtltype, ierr)
		  
		call InitTimer
     
     end subroutine StartMPI

     function get_current_time()
          real(dbpsn) :: get_current_time
           
#ifdef OPEN_MP
          get_current_time=omp_get_wtime()
#else		  		  
          get_current_time=MPI_Wtime()
#endif	
     end function get_current_time 

     subroutine StartTimer(tind)
          integer :: tind
          exec_time(tind)=get_current_time()	  
     end subroutine StartTimer

     subroutine StopTimer(tind)
          integer :: tind
          exec_time(tind)=get_current_time()	- exec_time(tind)		  
     end subroutine StopTimer
 
	subroutine InitTimer
		 integer :: n
		 do n=1,size(exec_time)
         	     exec_time(n)=get_current_time()				 
		 end do
	 end subroutine InitTimer
     
	 subroutine Abort(err_code)
		 integer :: err_code
		 if(err_code.eq.12) then 
			 if(proc.eq.0) print*,'Out of memory! Too many particles on GPU!'
		 end if
		 call MPI_Abort(MPI_COMM_WORLD, err_code, ierr)
	 end subroutine Abort
     
end module communication 