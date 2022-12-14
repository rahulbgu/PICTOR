module DerivedQnty
     use parameters
     use vars 
     use communication
contains 
     
!--------------------------------------------------------------------------------------------------     
! Mean Quantities: Averagin is done over the entire simualtion box   
!--------------------------------------------------------------------------------------------------          
     subroutine CalcSumVecFldSQ(Fldx,Fldy,Fldz,sum)
          real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz
          integer :: i,j,k
          integer :: mpi_err
          real(psn):: sum,sum_this
          sum_this=0.0_psn
#ifndef twoD 
          do k=3,mz-3
#else
        do k=1,1
#endif
          do j=3,my-3
               do i=3,mx-3
                    sum_this=sum_this+(Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2)
               end do
          end do
         end do
          call MPI_ALLREDUCE(sum_this,sum,1,mpi_psn,mpi_sum,MPI_COMM_WORLD,mpi_err)     
     end subroutine CalcSumVecFldSQ
         
     
end module DerivedQnty