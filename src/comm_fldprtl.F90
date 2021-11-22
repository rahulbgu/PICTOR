module comm_fldprtl
     use communication
#ifdef cyl
     use cyl_comm_fldprtl
#endif	 
     implicit none
                                        
contains 
	
!----------------------------------------------------------------------------------------------
! The following subroutines are used for tranferring particles between MPI proc. 
!----------------------------------------------------------------------------------------------
     
     subroutine SendRecvPrtlSize
          integer :: stat(MPI_STATUS_SIZE)     
          integer :: mpi_err          
          call MPI_SENDRECV(lpcross,1,MPI_INTEGER,lproc,0,rinp_count,1,MPI_INTEGER,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(rpcross,1,MPI_INTEGER,rproc,1,linp_count,1,MPI_INTEGER,lproc,1,MPI_COMM_WORLD,stat,mpi_err)          
          call MPI_SENDRECV(tpcross,1,MPI_INTEGER,tproc,2,binp_count,1,MPI_INTEGER,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(bpcross,1,MPI_INTEGER,bproc,3,tinp_count,1,MPI_INTEGER,tproc,3,MPI_COMM_WORLD,stat,mpi_err)
#ifndef twoD
        call MPI_SENDRECV(upcross,1,MPI_INTEGER,uproc,4,dinp_count,1,MPI_INTEGER,dproc,4,MPI_COMM_WORLD,stat,mpi_err)
        call MPI_SENDRECV(dpcross,1,MPI_INTEGER,dproc,5,uinp_count,1,MPI_INTEGER,uproc,5,MPI_COMM_WORLD,stat,mpi_err)
#endif           
!size of test particles 
          call MPI_SENDRECV(ltpcross,1,MPI_INTEGER,lproc,0,rintp_count,1,MPI_INTEGER,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(rtpcross,1,MPI_INTEGER,rproc,1,lintp_count,1,MPI_INTEGER,lproc,1,MPI_COMM_WORLD,stat,mpi_err)          
          call MPI_SENDRECV(ttpcross,1,MPI_INTEGER,tproc,2,bintp_count,1,MPI_INTEGER,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(btpcross,1,MPI_INTEGER,bproc,3,tintp_count,1,MPI_INTEGER,tproc,3,MPI_COMM_WORLD,stat,mpi_err)
#ifndef twoD
        call MPI_SENDRECV(utpcross,1,MPI_INTEGER,uproc,4,dintp_count,1,MPI_INTEGER,dproc,4,MPI_COMM_WORLD,stat,mpi_err)
        call MPI_SENDRECV(dtpcross,1,MPI_INTEGER,dproc,5,uintp_count,1,MPI_INTEGER,uproc,5,MPI_COMM_WORLD,stat,mpi_err)
#endif    
     end subroutine SendRecvPrtlSize
     subroutine SendRecvPrtl
          integer :: stat(MPI_STATUS_SIZE)     
          integer :: mpi_err     

          call MPI_SENDRECV(loutp(1:lcross),lcross,mpi_prtltype,lproc,0,rinp(1:rinp_count+rintp_count),rinp_count+rintp_count,mpi_prtltype,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(routp(1:rcross),rcross,mpi_prtltype,rproc,1,linp(1:linp_count+lintp_count),linp_count+lintp_count,mpi_prtltype,lproc,1,MPI_COMM_WORLD,stat,mpi_err)          
          call MPI_SENDRECV(toutp(1:tcross),tcross,mpi_prtltype,tproc,2,binp(1:binp_count+bintp_count),binp_count+bintp_count,mpi_prtltype,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(boutp(1:bcross),bcross,mpi_prtltype,bproc,3,tinp(1:tinp_count+tintp_count),tinp_count+tintp_count,mpi_prtltype,tproc,3,MPI_COMM_WORLD,stat,mpi_err)                         
#ifndef twoD
          call MPI_SENDRECV(uoutp(1:ucross),ucross,mpi_prtltype,uproc,4,dinp(1:dinp_count+dintp_count),dinp_count+dintp_count,mpi_prtltype,dproc,4,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(doutp(1:dcross),dcross,mpi_prtltype,dproc,5,uinp(1:uinp_count+uintp_count),uinp_count+uintp_count,mpi_prtltype,uproc,5,MPI_COMM_WORLD,stat,mpi_err)                         
#endif                    
     end subroutine SendRecvPrtl
	 
	 
	 
!----------------------------------------------------------------------------------------------
! The following subroutines are used to field communications	 
!----------------------------------------------------------------------------------------------
	 
	 
     
     subroutine ExchangeYZEdgeCurrent
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          dcount1=3*my*mz
          call MPI_SENDRECV(Jx(1:3,:,:),dcount1,mpi_psn,lproc,0,buff_rJx,dcount1,mpi_psn,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jx(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,buff_lJx,dcount1,mpi_psn,lproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jy(1:3,:,:),dcount1,mpi_psn,lproc,0,buff_rJy,dcount1,mpi_psn,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jy(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,buff_lJy,dcount1,mpi_psn,lproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jz(1:3,:,:),dcount1,mpi_psn,lproc,0,buff_rJz,dcount1,mpi_psn,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jz(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,buff_lJz,dcount1,mpi_psn,lproc,1,MPI_COMM_WORLD,stat,mpi_err)
#ifdef cyl 
          call  ExchangeYZEdgeCurrent_Axis
#endif 		     
     end subroutine ExchangeYZEdgeCurrent
     subroutine ExchangeZXEdgeCurrent
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          dcount1=3*mx*mz
          call MPI_SENDRECV(Jx(:,1:3,:),dcount1,mpi_psn,bproc,0,buff_tJx,dcount1,mpi_psn,tproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jx(:,my-2:my,:),dcount1,mpi_psn,tproc,1,buff_bJx,dcount1,mpi_psn,bproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jy(:,1:3,:),dcount1,mpi_psn,bproc,0,buff_tJy,dcount1,mpi_psn,tproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jy(:,my-2:my,:),dcount1,mpi_psn,tproc,1,buff_bJy,dcount1,mpi_psn,bproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jz(:,1:3,:),dcount1,mpi_psn,bproc,0,buff_tJz,dcount1,mpi_psn,tproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jz(:,my-2:my,:),dcount1,mpi_psn,tproc,1,buff_bJz,dcount1,mpi_psn,bproc,1,MPI_COMM_WORLD,stat,mpi_err)     
     end subroutine ExchangeZXEdgeCurrent
     subroutine ExchangeXYEdgeCurrent
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          dcount1=3*mx*my
          call MPI_SENDRECV(Jx(:,:,1:3),dcount1,mpi_psn,dproc,0,buff_uJx,dcount1,mpi_psn,uproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jx(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,buff_dJx,dcount1,mpi_psn,dproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jy(:,:,1:3),dcount1,mpi_psn,dproc,0,buff_uJy,dcount1,mpi_psn,uproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jy(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,buff_dJy,dcount1,mpi_psn,dproc,1,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Jz(:,:,1:3),dcount1,mpi_psn,dproc,0,buff_uJz,dcount1,mpi_psn,uproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Jz(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,buff_dJz,dcount1,mpi_psn,dproc,1,MPI_COMM_WORLD,stat,mpi_err)     
     end subroutine ExchangeXYEdgeCurrent
     
     
     subroutine ExchangeYZEdgeField(Fldx,Fldy,Fldz)
          real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz 
          integer :: dcount1,dcount2,mpi_err
          integer :: stat(MPI_STATUS_SIZE)
		  dcount1=3*my*mz
          dcount2=2*my*mz
          call MPI_SENDRECV(Fldx(3:5,:,:),dcount1,mpi_psn,lproc,1,Fldx(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldx(mx-4:mx-3,:,:),dcount2,mpi_psn,rproc,2,Fldx(1:2,:,:),dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldy(3:5,:,:),dcount1,mpi_psn,lproc,1,Fldy(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldy(mx-4:mx-3,:,:),dcount2,mpi_psn,rproc,2,Fldy(1:2,:,:),dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldz(3:5,:,:),dcount1,mpi_psn,lproc,1,Fldz(mx-2:mx,:,:),dcount1,mpi_psn,rproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldz(mx-4:mx-3,:,:),dcount2,mpi_psn,rproc,2,Fldz(1:2,:,:),dcount2,mpi_psn,lproc,2,MPI_COMM_WORLD,stat,mpi_err)
#ifdef cyl
          call ExchangeYZEdgeField_Axis(Fldx,Fldy,Fldz) 
#endif 		  
                    
     end subroutine ExchangeYZEdgeField
     subroutine ExchangeZXEdgeField(Fldx,Fldy,Fldz)
          real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz 
          integer :: dcount1,dcount2,mpi_err
          integer :: stat(MPI_STATUS_SIZE)
          dcount1=3*mx*mz
          dcount2=2*mx*mz
          call MPI_SENDRECV(Fldx(:,3:5,:),dcount1,mpi_psn,bproc,1,Fldx(:,my-2:my,:),dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldx(:,my-4:my-3,:),dcount2,mpi_psn,tproc,2,Fldx(:,1:2,:),dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldy(:,3:5,:),dcount1,mpi_psn,bproc,1,Fldy(:,my-2:my,:),dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldy(:,my-4:my-3,:),dcount2,mpi_psn,tproc,2,Fldy(:,1:2,:),dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldz(:,3:5,:),dcount1,mpi_psn,bproc,1,Fldz(:,my-2:my,:),dcount1,mpi_psn,tproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldz(:,my-4:my-3,:),dcount2,mpi_psn,tproc,2,Fldz(:,1:2,:),dcount2,mpi_psn,bproc,2,MPI_COMM_WORLD,stat,mpi_err)               
     end subroutine ExchangeZXEdgeField
     subroutine ExchangeXYEdgeField(Fldx,Fldy,Fldz)
          real(psn), dimension(mx,my,mz) :: Fldx,Fldy,Fldz 
          integer :: dcount1,dcount2,mpi_err
          integer :: stat(MPI_STATUS_SIZE)
          dcount1=3*mx*my
          dcount2=2*mx*my
          call MPI_SENDRECV(Fldx(:,:,3:5),dcount1,mpi_psn,dproc,1,Fldx(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldx(:,:,mz-4:mz-3),dcount2,mpi_psn,uproc,2,Fldx(:,:,1:2),dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldy(:,:,3:5),dcount1,mpi_psn,dproc,1,Fldy(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldy(:,:,mz-4:mz-3),dcount2,mpi_psn,uproc,2,Fldy(:,:,1:2),dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
          
          call MPI_SENDRECV(Fldz(:,:,3:5),dcount1,mpi_psn,dproc,1,Fldz(:,:,mz-2:mz),dcount1,mpi_psn,uproc,1,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(Fldz(:,:,mz-4:mz-3),dcount2,mpi_psn,uproc,2,Fldz(:,:,1:2),dcount2,mpi_psn,dproc,2,MPI_COMM_WORLD,stat,mpi_err)
     end subroutine ExchangeXYEdgeField

    !The following subroutines are used to comminucate only the immediate outer current layer 
     subroutine ExchangeYZEdgeCurrent1(J1)
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          real(psn),dimension(mx,my,mz) :: J1
          dcount1=my*mz
          call MPI_SENDRECV(J1(3:3,:,:),dcount1,mpi_psn,lproc,0,J1(mx-2:mx-2,:,:),dcount1,mpi_psn,rproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(J1(mx-3:mx-3,:,:),dcount1,mpi_psn,rproc,1,J1(2:2,:,:),dcount1,mpi_psn,lproc,1,MPI_COMM_WORLD,stat,mpi_err)

     end subroutine ExchangeYZEdgeCurrent1
     
     subroutine ExchangeZXEdgeCurrent1(J1)
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          real(psn),dimension(mx,my,mz) :: J1
          dcount1=mx*mz
          call MPI_SENDRECV(J1(:,3:3,:),dcount1,mpi_psn,bproc,0,J1(:,my-2:my-2,:),dcount1,mpi_psn,tproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(J1(:,my-3:my-3,:),dcount1,mpi_psn,tproc,1,J1(:,2:2,:),dcount1,mpi_psn,bproc,1,MPI_COMM_WORLD,stat,mpi_err)          
     end subroutine ExchangeZXEdgeCurrent1
     
     subroutine ExchangeXYEdgeCurrent1(J1)
          integer :: stat(MPI_STATUS_SIZE)
          integer :: dcount1,mpi_err
          real(psn),dimension(mx,my,mz) :: J1
          dcount1=mx*my
          call MPI_SENDRECV(J1(:,:,3:3),dcount1,mpi_psn,dproc,0,J1(:,:,mz-2:mz-2),dcount1,mpi_psn,uproc,0,MPI_COMM_WORLD,stat,mpi_err)
          call MPI_SENDRECV(J1(:,:,mz-3:mz-3),dcount1,mpi_psn,uproc,1,J1(:,:,2:2),dcount1,mpi_psn,dproc,1,MPI_COMM_WORLD,stat,mpi_err)          
     end subroutine ExchangeXYEdgeCurrent1

     
     !the following subroutines are useful in case of smoothening the EM field 
     subroutine ExchangeYZEdgeEMFld1(Fld) !Not needed 
          integer :: stat(MPI_STATUS_SIZE)
          real(psn),dimension(mx,my,mz) :: Fld
          call MPI_SENDRECV(Fld(3:3,:,:),my*mz,mpi_psn,lproc,0,Fld(mx-2:mx-2,:,:),my*mz,mpi_psn,rproc,0,MPI_COMM_WORLD,stat,ierr)
          call MPI_SENDRECV(Fld(mx-3:mx-3,:,:),my*mz,mpi_psn,rproc,1,Fld(2:2,:,:),my*mz,mpi_psn,lproc,1,MPI_COMM_WORLD,stat,ierr)
     end subroutine ExchangeYZEdgeEMFld1
     
     subroutine ExchangeZXEdgeEMFld1(Fld) !not needed 
          integer :: stat(MPI_STATUS_SIZE)
          real(psn),dimension(mx,my,mz) :: Fld
          call MPI_SENDRECV(Fld(:,3:3,:),mx*mz,mpi_psn,bproc,0,Fld(:,my-2:my-2,:),mx*mz,mpi_psn,tproc,0,MPI_COMM_WORLD,stat,ierr)
          call MPI_SENDRECV(Fld(:,my-3:my-3,:),mx*mz,mpi_psn,tproc,1,Fld(:,2:2,:),mx*mz,mpi_psn,bproc,1,MPI_COMM_WORLD,stat,ierr)     
     end subroutine ExchangeZXEdgeEMFld1
     
end module comm_fldprtl