module communication 
     use parameters
     use vars
     implicit none 
     include 'mpif.h'

    integer :: mpi_psn,mpi_prtltype
     
contains 
     subroutine StartMPI
          integer :: extent,nblocks
          integer,dimension(0:7) :: oldtypes, blockcounts, offsets
          call mpi_init(ierr)
          call mpi_comm_rank(MPI_COMM_WORLD,proc,ierr)
          call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
          
          !  a) set the precision for mpi data transfer b)Create data type for particle transfer 
          ! **Warning**: The following must be modified if the structure of particle datatype is changed  
          nblocks=0
          select case(psn)
          case(kind(1.0d0))
             mpi_psn=MPI_DOUBLE_PRECISION !mpi precision 
			 
             offsets(nblocks)=0
             oldtypes(nblocks)=MPI_REAL8
             blockcounts(nblocks)=8
             call MPI_TYPE_EXTENT(MPI_REAL8,extent,ierr)
             nblocks=nblocks+1
			 
             offsets(nblocks)=offsets(nblocks-1)+blockcounts(nblocks-1)*extent
             oldtypes(nblocks)=MPI_INTEGER 
             blockcounts(nblocks)=2
             call MPI_TYPE_EXTENT(MPI_INTEGER,extent,ierr)
             nblocks=nblocks+1  

          case(kind(1.0e0))
             mpi_psn=MPI_REAL !mpi precision
           
             offsets(nblocks)=0
             oldtypes(nblocks)=MPI_REAL 
             blockcounts(nblocks)=8
             call MPI_TYPE_EXTENT(MPI_REAL,extent,ierr)
             nblocks=nblocks+1       
			 
             offsets(nblocks)=offsets(nblocks-1)+blockcounts(nblocks-1)*extent
             oldtypes(nblocks)=MPI_INTEGER 
             blockcounts(nblocks)=2
             call MPI_TYPE_EXTENT(MPI_INTEGER,extent,ierr)
             nblocks=nblocks+1         
          end select           
          
          call MPI_TYPE_STRUCT(nblocks,blockcounts,offsets,oldtypes,mpi_prtltype,ierr)
          call MPI_TYPE_COMMIT(mpi_prtltype, ierr)
     
     end subroutine StartMPI
     subroutine StartTimer(tind)
          integer :: tind
#ifdef OPEN_MP
          exec_time(tind)=omp_get_wtime()
#else		  		  
          exec_time(tind)=MPI_Wtime()
#endif		  
     end subroutine StartTimer
     subroutine StopTimer(tind)
          integer :: tind
#ifdef OPEN_MP		  
          exec_time(tind)=omp_get_wtime()-exec_time(tind)
#else		  
          exec_time(tind)=MPI_Wtime()-exec_time(tind)
#endif		  
     end subroutine StopTimer
     
	 subroutine Abort(err_code)
		 integer :: err_code
		 if(err_code.eq.12) then 
			 if(proc.eq.0) print*,'Low memory! Too many particles on GPU!'
		 end if
		 call MPI_Abort(MPI_COMM_WORLD, err_code, ierr)
	 end subroutine Abort
     
end module communication 