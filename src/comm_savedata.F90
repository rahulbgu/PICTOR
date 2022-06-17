module comm_savedata
     use parameters
     use vars
     use communication
     use savedata_routines
     implicit none 
contains 
     
     subroutine GatherEnergy
          real(psn), dimension(4) :: energy_this
          energy=0
          call CalcEnergy(energy_this)
          call MPI_ALLREDUCE(energy_this,energy,4,mpi_psn,mpi_sum,MPI_COMM_WORLD,ierr)
          !if(proc.eq.0) print*,'Total Energy',energy(1)+energy(2)+energy(3)+energy(4)
     end subroutine GatherEnergy
     !---------------------------------------------------------------------------------------------
     !The following subroutines are to facilitate saving particle data 
     !---------------------------------------------------------------------------------------------

     subroutine ReduceToSavePrtlSize
           call MPI_ALLGATHER(tosave_prtl_arr_size,1,MPI_INTEGER,prtl_arr_size_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
     end subroutine ReduceToSavePrtlSize
     

! !----The following subroutine are for corsscheck purpose, should not be part of the final version ---------------
        subroutine GetTotalNP
               integer :: np_total,np_temp,i,np0
               type(MPI_Status) :: stat  
               save np0
               if(proc.eq.0) then
                    np_total=np
                    do i=1,nproc-1
                       call MPI_RECV(np_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
                       np_total=np_total+np_temp
                   end do
                   print*,'Total Number of particle in the Simulations at step no.',t,' is ::',np_total
                    if(t.eq.0) np0=np_total
                    !if(np_total.ne.np0) STOP 'Some Particles have gone missing'
               else
                    call MPI_SEND(np,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
               end if
        end subroutine GetTotalNP
!         subroutine GetTotalPositiveA
!                integer :: a_total,a_temp,i,j,a0
!                type(MPI_Status) :: stat  
!                save a0
!                a_total=0
!                do j=1,prtl_arr_size
!                     if(p(j)%a.ne.0) a_total=a_total+1
!                end do
!
!                if(proc.eq.0) then
!                     do i=1,nproc-1
!                        call MPI_RECV(a_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
!                        a_total=a_total+a_temp
!                    end do
!                    print*,'Total Number of ACTIVE particles in the Simulations at step no.',t,' is ::',a_total
!                     if(t.eq.1) a0=a_total
!                     if(a_total.ne.a0) STOP 'No. of active particles have gone missing'
!                else
!                     call MPI_SEND(a_total,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
!                end if
!         end subroutine GetTotalPositiveA
!           subroutine GetTotalIon
!                integer :: ion_total,ion_temp,i,j,ion0
!                type(MPI_Status) :: stat  
!                save ion0
!                ion_total=0
!                do j=1,prtl_arr_size
!                     if(p(j)%flv.eq.1) ion_total=ion_total+1
!                end do
!
!                if(proc.eq.0) then
!                     do i=1,nproc-1
!                        call MPI_RECV(ion_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
!                        ion_total=ion_total+ion_temp
!                    end do
!                    print*,'Total Number of Ions in the Simulations at step no.',t,' is ::',ion_total
!                     if(t.eq.0) ion0=ion_total
!                     if(ion_total.ne.ion0) STOP 'No. of ions have gone missing'
!                else
!                     call MPI_SEND(ion_total,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
!                end if
!           end subroutine GetTotalIon
          subroutine GetPositiveQ
               integer :: q_total,q_temp,i,j,q0
               type(MPI_Status) :: stat  
               save q0
               q_total=0
               do j=1,used_prtl_arr_size
                    if(qp(j).gt.0) q_total=q_total+1
               end do

               if(proc.eq.0) then
                    do i=1,nproc-1
                       call MPI_RECV(q_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
                       q_total=q_total+q_temp
                   end do
                   print*,'Total Number of positive Q in the Simulations at step no.',t,' is ::',q_total
                    if(t.eq.0) q0=q_total
                    if(q_total.ne.q0) STOP 'No. of postive Q has changed'
               else
                    call MPI_SEND(q_total,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
               end if

          end subroutine GetPositiveQ
          subroutine GetNegativeQ
               integer :: q_total,q_temp,i,j,q0
               type(MPI_Status) :: stat  
               save q0
               q_total=0
               do j=1,used_prtl_arr_size
                    if(qp(j).eq.-1.0) q_total=q_total+1
               end do

               if(proc.eq.0) then
                    do i=1,nproc-1
                       call MPI_RECV(q_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
                       q_total=q_total+q_temp
                   end do
                   print*,'Total Number of Negative Q in the Simulations at step no.',t,' is ::',q_total
                    if(t.eq.0) q0=q_total
                    !if(q_total.ne.q0) STOP 'No. of postive Q has changed'
               else
                    call MPI_SEND(q_total,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
               end if

          end subroutine GetNegativeQ

end module comm_savedata