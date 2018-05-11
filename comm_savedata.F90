module comm_savedata
     use parameters
     use vars
     use communication
     use savedata_routines
     implicit none 
contains 
     subroutine SendFldSizeToMaster
          integer,dimension(3) :: fld_size_this
          fld_size_this(1)=fdatax
          fld_size_this(2)=fdatay
          fld_size_this(3)=fdataz
          call MPI_SEND(fld_size_this,3,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)     
     end subroutine SendFldSizeToMaster
     subroutine RecvFldSizeAtMaster(proc_no)
          integer :: proc_no
          integer, dimension(MPI_STATUS_SIZE) :: stat
          call MPI_RECV(fld_size_all(proc_no,:),3,MPI_INTEGER,proc_no,0,MPI_COMM_WORLD,stat,ierr)
     end subroutine RecvFldSizeAtMaster
     subroutine SendEMFldToMaster
          integer :: fdsize 
          fdsize=fdatax*fdatay*fdataz
          call CollectFld(Ex,1)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,1,MPI_COMM_WORLD,ierr)
          call CollectFld(Ey,3)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,2,MPI_COMM_WORLD,ierr)
          call CollectFld(Ez,5)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,3,MPI_COMM_WORLD,ierr)
          call CollectFld(Bx,7)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,4,MPI_COMM_WORLD,ierr)
          call CollectFld(By,9)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,5,MPI_COMM_WORLD,ierr)
          call CollectFld(Bz,11)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,6,MPI_COMM_WORLD,ierr)
     end subroutine SendEMFldToMaster
     subroutine SendCurrToMaster
          integer :: fdsize 
          fdsize=fdatax*fdatay*fdataz
          !currents are interpolated to grid points in the same way as E 
          call CollectFld(Jx,1)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,7,MPI_COMM_WORLD,ierr)
          call CollectFld(Jy,3)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,8,MPI_COMM_WORLD,ierr)
          call CollectFld(Jz,5)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,9,MPI_COMM_WORLD,ierr)
          
     end subroutine SendCurrToMaster
     subroutine SendDensityToMaster
          integer :: fdsize 
          fdsize=fdatax*fdatay*fdataz
          call CollectFld(Jx,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,10,MPI_COMM_WORLD,ierr)
          call CollectFld(Jy,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,11,MPI_COMM_WORLD,ierr)
               
     end subroutine SendDensityToMaster
     subroutine SendVelFldToMaster_ion
          integer :: fdsize 
          fdsize=fdatax*fdatay*fdataz
          call CollectFld(Jx,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,12,MPI_COMM_WORLD,ierr)
          call CollectFld(Jy,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,13,MPI_COMM_WORLD,ierr)
          call CollectFld(Jz,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,14,MPI_COMM_WORLD,ierr)
          
     end subroutine SendVelFldToMaster_ion
     subroutine SendVelFldToMaster_elec
          integer :: fdsize 
          fdsize=fdatax*fdatay*fdataz
          call CollectFld(Jx,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,15,MPI_COMM_WORLD,ierr)
          call CollectFld(Jy,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,16,MPI_COMM_WORLD,ierr)
          call CollectFld(Jz,13)
          call MPI_SEND(fdata,fdsize,mpi_psn,0,17,MPI_COMM_WORLD,ierr)
          
     end subroutine SendVelFldToMaster_elec
     subroutine senddivEToMaster
          call CalcDivE
          call CollectFld(Jx,13)
          call MPI_SEND(fdata,fdatax*fdatay*fdataz,mpi_psn,0,18,MPI_COMM_WORLD,ierr)
     end subroutine senddivEToMaster
     
     subroutine RecvFldAtMaster(proc_no,tag)
          integer :: fs1,fs2,fs3,proc_no,tag
          integer, dimension(MPI_STATUS_SIZE) :: stat
          fs1=fld_size_all(proc_no,1)
          fs2=fld_size_all(proc_no,2)
          fs3=fld_size_all(proc_no,3)
        allocate(fdata(fs1,fs2,fs3))
          call MPI_RECV(fdata,fs1*fs2*fs3,mpi_psn,proc_no,tag,MPI_COMM_WORLD,stat,ierr)
     end subroutine RecvFldAtMaster

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
          subroutine SendPrtlToMaster
                call CollectPrtl(1)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,1,MPI_COMM_WORLD,ierr) !send x data first
               call CollectPrtl(2)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,2,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(3)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,3,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(4)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,4,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(5)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,5,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(6)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,6,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(7)
               call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,7,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(8)
               call MPI_SEND(pdata_int,tosave_prtl_arr_size,MPI_INTEGER,0,8,MPI_COMM_WORLD,ierr) 
#ifdef mulflvr               
                call CollectPrtl(9)
               call MPI_SEND(pdata_int,tosave_prtl_arr_size,MPI_INTEGER,0,9,MPI_COMM_WORLD,ierr) 
                call CollectPrtl(10)
               call MPI_SEND(pdata_int,tosave_prtl_arr_size,MPI_INTEGER,0,10,MPI_COMM_WORLD,ierr) 
#endif
               if(save_prtl_local_fld) then
                     call CollectPrtl(11)
                    call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,11,MPI_COMM_WORLD,ierr) 
                    call CollectPrtl(12)
                   call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,12,MPI_COMM_WORLD,ierr) 
                    call CollectPrtl(13)
                     call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,13,MPI_COMM_WORLD,ierr) 
                     call CollectPrtl(14)
                    call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,14,MPI_COMM_WORLD,ierr) 
                    call CollectPrtl(15)
                   call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,15,MPI_COMM_WORLD,ierr) 
                    call CollectPrtl(16)
                     call MPI_SEND(pdata_real,tosave_prtl_arr_size,mpi_psn,0,16,MPI_COMM_WORLD,ierr) 
               end if
          
          end subroutine SendPrtlToMaster
          subroutine SendPrtlSizeToMaster
               call MPI_SEND(tosave_prtl_arr_size,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)          
          end subroutine SendPrtlSizeToMaster
          subroutine RecvPrtlSizeAtMaster(proc_no)
               integer :: proc_no
               integer, dimension(MPI_STATUS_SIZE) :: stat
               call MPI_RECV(prtl_arr_size_all(proc_no),1,MPI_INTEGER,proc_no,0,MPI_COMM_WORLD,stat,ierr)
          end subroutine RecvPrtlSizeAtMaster
          subroutine RecvPrtlAtMaster(proc_no,tag,dtype)
               integer :: proc_no,tag,dtype
               integer, dimension(MPI_STATUS_SIZE)::stat
               select case (dtype)
                case(1)
                   allocate(pdata_real(prtl_arr_size_all(proc_no)))
                 call MPI_RECV(pdata_real,prtl_arr_size_all(proc_no),mpi_psn,proc_no,tag,MPI_COMM_WORLD,stat,ierr)               
                 case(2)
                      allocate(pdata_int(prtl_arr_size_all(proc_no)))
                    call MPI_RECV(pdata_int,prtl_arr_size_all(proc_no),MPI_INTEGER,proc_no,tag,MPI_COMM_WORLD,stat,ierr)
               end select                
          end subroutine RecvPrtlAtMaster

          subroutine ReduceToSavePrtlSize
               call MPI_ALLGATHER(tosave_prtl_arr_size,1,MPI_INTEGER,prtl_arr_size_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
          end subroutine ReduceToSavePrtlSize
     

! !----The following subroutine are for corsscheck purpose, should not be part of the final version ---------------
        subroutine GetTotalNP
               integer :: np_total,np_temp,i,np0
               integer, dimension(MPI_STATUS_SIZE) :: stat
               save np0
               if(proc.eq.0) then
                    np_total=np
                    do i=1,nproc-1
                       call MPI_RECV(np_temp,1,MPI_INTEGER,i,i,MPI_COMM_WORLD,stat,ierr)
                       np_total=np_total+np_temp
                   end do
                   print*,'Total Number of particle in the Simulations at step no.',t,' is ::',np_total
                    if(t.eq.0) np0=np_total
                    if(np_total.ne.np0) STOP 'Some Particles have gone missing'
               else
                    call MPI_SEND(np,1,MPI_INTEGER,0,proc,MPI_COMM_WORLD,ierr)
               end if
        end subroutine GetTotalNP
!         subroutine GetTotalPositiveA
!                integer :: a_total,a_temp,i,j,a0
!                integer, dimension(MPI_STATUS_SIZE) :: stat
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
!                integer, dimension(MPI_STATUS_SIZE) :: stat
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
               integer, dimension(MPI_STATUS_SIZE) :: stat
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
               integer, dimension(MPI_STATUS_SIZE) :: stat
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