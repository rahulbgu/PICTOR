!-----------------------------------------------------
!subroutines used to generate vel and gamma spectrum
!-----------------------------------------------------
module prtl_stats
	use parameters
	use vars
	use communication
    use interpolation, only : InterpEMfield	
	implicit none 
	
    real(dbpsn),dimension(:,:), allocatable:: spec_gamma 
    real,dimension(:),   allocatable:: Gamma_spec_bin
    integer :: Gamma_spec_binlen
    real(psn) :: gmax_allflv,gmin_allflv,local_gmax_allflv,local_gmin_allflv 
    
	!varaibles used to genrate prtl velocity spectrum 
    real(dbpsn),dimension(:,:), allocatable:: spec_speed 
    real,dimension(:),   allocatable:: Speed_spec_bin
    integer :: Speed_spec_binlen
	
    !varaibles used to generate mean of various prtl averaged quantities 
    real(dbpsn),dimension(:), allocatable :: SumNprtl,SumQ, SumQxKE, SumQxExVx,SumQxEyVy,SumQxEzVz,SumQxPx2,SumQxPy2,SumQxPz2
    real(dbpsn),dimension(:,:),allocatable :: Speed_SumQxExVx,Speed_SumQxEyVy,Speed_SumQxEzVz 
    !real(dbpsn),dimension(:), allocatable  :: SumQxLR
contains 


subroutine InitPrtlStatVars
     call InitSpeedSpecBin !defined in savedata_routines
     allocate(SumNprtl(Nflvr),SumQ(Nflvr),SumQxKE(Nflvr))  !used for saving prtl averaged info  
     allocate(SumQxExVx(Nflvr),SumQxEyVy(Nflvr),SumQxEzVz(Nflvr))
     allocate(SumQxPx2(Nflvr),SumQxPy2(Nflvr),SumQxPz2(Nflvr))
     allocate(Speed_SumQxExVx(Nflvr,Speed_spec_binlen),Speed_SumQxEyVy(Nflvr,Speed_spec_binlen),Speed_SumQxEzVz(Nflvr,Speed_spec_binlen)) 
    !allocate(SumQxLR(Nflvr))
     allocate(spec_gamma(Nflvr,0),Gamma_spec_bin(0))
end subroutine InitPrtlStatVars
		
!-------------------------------------------------------------------------------------------------
!Allocate arrays according to the current bin size and range to be covered 
!-------------------------------------------------------------------------------------------------
subroutine CreateGammaSpecBin     
  integer :: i
  Gamma_spec_binlen=ceiling((log(gmax_allflv)-log(gmin_allflv))/Gamma_spec_binwidth)+1
  if(size(Gamma_spec_bin).lt.Gamma_spec_binlen) then 
      deallocate(spec_gamma,Gamma_spec_bin)
	  allocate(spec_gamma(Nflvr,Gamma_spec_binlen),Gamma_spec_bin(Gamma_spec_binlen)) 
	  do i=1,Gamma_spec_binlen
	      Gamma_spec_bin(i)=gmin_allflv*exp(real(i-1)*Gamma_spec_binwidth)
	  end do
   end if
end subroutine CreateGammaSpecBin

subroutine InitSpeedSpecBin     
  integer :: i
  Speed_spec_binlen=ceiling(1.0_psn/Speed_spec_binwidth)+1
  allocate(spec_speed(Nflvr,Speed_spec_binlen),Speed_spec_bin(Speed_spec_binlen)) 
  do i=1,Speed_spec_binlen
      Speed_spec_bin(i)=Speed_spec_binwidth*((i-1)**2)
  end do
end subroutine InitSpeedSpecBin


!-------------------------------------------------------------------------------------------------
!The following subroutines are used to calculate spectrum of all particles in the simulation box
!-------------------------------------------------------------------------------------------------
subroutine CalcGmaxGminLocal_all_prtl
  integer   :: n
  real(psn) :: gamma 
  local_gmin_allflv=1.0_psn
  local_gmax_allflv=1.0_psn
  do n=1,used_prtl_arr_size          
       if(qp(n).eq.0) cycle
       gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
       if(gamma.gt.local_gmax_allflv) local_gmax_allflv=gamma
       if(gamma.lt.local_gmin_allflv) local_gmin_allflv=gamma !optmisation possible      
  end do       
  do n=1,used_test_prtl_arr_size          
       if(qtp(n).eq.0) cycle
       gamma=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
       if(gamma.gt.local_gmax_allflv) local_gmax_allflv=gamma
       if(gamma.lt.local_gmin_allflv) local_gmin_allflv=gamma !optmisation possible      
  end do   
end subroutine CalcGmaxGminLocal_all_prtl

subroutine CalcGammaSpectrum_all_prtl
  integer :: n,ch,bin_ind
  real(psn) ::gamma
  spec_gamma=0.0_psn
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle
       ch=flvp(n)
       gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
       bin_ind=floor((log(gamma)-log(gmin_allflv))/Gamma_spec_binwidth)+1
       spec_gamma(ch,bin_ind)=spec_gamma(ch,bin_ind)+ abs(qp(n))
  end do
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle
       ch=flvtp(n)
       gamma=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
       bin_ind=floor((log(gamma)-log(gmin_allflv))/Gamma_spec_binwidth)+1
       spec_gamma(ch,bin_ind)=spec_gamma(ch,bin_ind)+ abs(qtp(n))
  end do
end subroutine CalcGammaSpectrum_all_prtl

subroutine CalcSpeedSpectrum_all_prtl
  integer :: n,ch,bin_ind
  real(psn) ::speed
  spec_speed=0
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle
       ch=flvp(n)
       speed=1.0_psn - 1.0_psn/(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n))               
       bin_ind=floor(speed/Speed_spec_binwidth) +1
       spec_speed(ch,bin_ind)=spec_speed(ch,bin_ind)+abs(qp(n))
  end do
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle
       ch=flvtp(n)
       speed=1.0_psn - 1.0_psn/(1.0_psn+utp(n)*utp(n)+vtp(n)*vtp(n)+wtp(n)*wtp(n))               
       bin_ind=floor(speed/Speed_spec_binwidth) +1
       spec_speed(ch,bin_ind)=spec_speed(ch,bin_ind)+abs(qtp(n))
  end do
end subroutine CalcSpeedSpectrum_all_prtl

!-------------------------------------------------------------------------------------------------
!spectrum of particles in a subdomain, x1,..  are in global cordinate 
!-------------------------------------------------------------------------------------------------
subroutine CalcGmaxGminLocalInSubDomain(xi,xf,yi,yf,zi,zf)
  integer   :: xi,xf,yi,yf,zi,zf
  integer   :: x1,x2,y1,y2,z1,z2
  integer   :: n
  real(psn) :: gamma
  
  call DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
  
  local_gmin_allflv=1.0_psn
  local_gmax_allflv=1.0_psn
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle
	   if((xp(n).lt.x1.or.xp(n).gt.x2).or.(yp(n).lt.y1.or.yp(n).gt.y2).or.(zp(n).lt.z1.or.zp(n).gt.z2)) cycle
       gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
       if(gamma.gt.local_gmax_allflv) local_gmax_allflv=gamma
       if(gamma.lt.local_gmin_allflv) local_gmin_allflv=gamma !optmisation possible
  end do
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle
	   if((xtp(n).lt.x1.or.xtp(n).gt.x2).or.(ytp(n).lt.y1.or.ytp(n).gt.y2).or.(ztp(n).lt.z1.or.ztp(n).gt.z2)) cycle
       gamma=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
       if(gamma.gt.local_gmax_allflv) local_gmax_allflv=gamma
       if(gamma.lt.local_gmin_allflv) local_gmin_allflv=gamma !optmisation possible
  end do
end subroutine CalcGmaxGminLocalInSubDomain

subroutine CalcGammaSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
  integer   :: xi,xf,yi,yf,zi,zf
  integer   :: x1,x2,y1,y2,z1,z2
  integer :: n,ch,bin_ind
  real(psn) ::gamma
  
  call DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
  
  spec_gamma=0.0_psn
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle
	   if((xp(n).lt.x1.or.xp(n).gt.x2).or.(yp(n).lt.y1.or.yp(n).gt.y2).or.(zp(n).lt.z1.or.zp(n).gt.z2)) cycle
       ch=flvp(n)
       gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
       bin_ind=floor((log(gamma)-log(gmin_allflv))/Gamma_spec_binwidth)+1
       spec_gamma(ch,bin_ind)=spec_gamma(ch,bin_ind)+ abs(qp(n))
  end do
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle
	   if((xtp(n).lt.x1.or.xtp(n).gt.x2).or.(ytp(n).lt.y1.or.ytp(n).gt.y2).or.(ztp(n).lt.z1.or.ztp(n).gt.z2)) cycle
       ch=flvtp(n)
       gamma=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
       bin_ind=floor((log(gamma)-log(gmin_allflv))/Gamma_spec_binwidth)+1
       spec_gamma(ch,bin_ind)=spec_gamma(ch,bin_ind)+ abs(qtp(n))
  end do
end subroutine CalcGammaSpectrumInSubDomain

subroutine CalcSpeedSpectrumInSubDomain(xi,xf,yi,yf,zi,zf)
  integer   :: xi,xf,yi,yf,zi,zf
  integer   :: x1,x2,y1,y2,z1,z2
  integer :: n,ch,bin_ind
  real(psn) ::speed
  
  call DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
  
  spec_speed=0
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle
       ch=flvp(n)
       speed=1.0_psn - 1.0_psn/(1.0_psn+up(n)*up(n)+vp(n)*vp(n)+wp(n)*wp(n))
       bin_ind=floor(speed/Speed_spec_binwidth) +1
       spec_speed(ch,bin_ind)=spec_speed(ch,bin_ind)+abs(qp(n))
  end do
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle
       ch=flvtp(n)
       speed=1.0_psn - 1.0_psn/(1.0_psn+utp(n)*utp(n)+vtp(n)*vtp(n)+wtp(n)*wtp(n))
       bin_ind=floor(speed/Speed_spec_binwidth) +1
       spec_speed(ch,bin_ind)=spec_speed(ch,bin_ind)+abs(qtp(n))
  end do
end subroutine CalcSpeedSpectrumInSubDomain

subroutine DomainBoundary_GlobalToLocalCord(xi,xf,yi,yf,zi,zf,x1,x2,y1,y2,z1,z2)
    integer, intent(IN)  :: xi,xf,yi,yf,zi,zf
    integer, intent(INOUT)  :: x1,x2,y1,y2,z1,z2
    x1=xi-xborders(procxind(proc))+3
    x2=xf-xborders(procxind(proc))+3
    y1=yi-yborders(procyind(proc))+3
    y2=yf-yborders(procyind(proc))+3 
#ifdef twoD
    z1=1+zi !range of local cordinate is [1,2]
	z2=1+zf
#else 	 	 
    z1=zi-zborders(proczind(proc))+3
    z2=zf-zborders(proczind(proc))+3
#endif 
end subroutine DomainBoundary_GlobalToLocalCord
	 
!-------------------------------------------------------------------------------------------------------
! The following subroutines are used to computed mean values of several particle related quantities, such as mean KE
!-------------------------------------------------------------------------------------------------------

subroutine SumPrtlQntyMain
  integer :: n,bin_ind
  real(psn) :: gamma,invg,speed
  real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz

  SumNprtl=0
  SumQ=0
  SumQxKE=0
  SumQxExVx=0
  SumQxEyVy=0
  SumQxEzVz=0
  SumQxPx2=0
  SumQxPy2=0
  SumQxPz2=0
  Speed_SumQxExVx=0
  Speed_SumQxEyVy=0
  Speed_SumQxEzVz=0
  do n=1,used_prtl_arr_size
       if(qp(n).eq.0) cycle

       SumNprtl(flvp(n))=SumNprtl(flvp(n))+1
       SumQ(flvp(n))=SumQ(flvp(n))+qp(n)
       gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
       invg=1.0_psn/gamma
       SumQxKE(flvp(n))=SumQxKE(flvp(n))+qp(n)*(gamma-1) ! kinetic energy is in units of m_0c^2 (m_0=rest mass)
  
       call InterpEMfield(xp(n),yp(n),zp(n),pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)
       SumQxExVx(flvp(n))=SumQxExVx(flvp(n))+qp(n)*up(n)*invg*pEx
       SumQxEyVy(flvp(n))=SumQxEyVy(flvp(n))+qp(n)*vp(n)*invg*pEy
       SumQxEzVz(flvp(n))=SumQxEzVz(flvp(n))+qp(n)*wp(n)*invg*pEz
       SumQxPx2(flvp(n)) =SumQxPx2(flvp(n)) +qp(n)*up(n)**2
       SumQxPy2(flvp(n)) =SumQxPy2(flvp(n)) +qp(n)*vp(n)**2
       SumQxPz2(flvp(n)) =SumQxPz2(flvp(n)) +qp(n)*wp(n)**2
  
       speed=sqrt((gamma-1.0_psn)*(gamma+1.0_psn))*invg
       bin_ind= floor(speed/Speed_spec_binwidth)+1
       Speed_SumQxExVx(flvp(n),bin_ind)=Speed_SumQxExVx(flvp(n),bin_ind)+qp(n)*up(n)*invg*pEx
       Speed_SumQxEyVy(flvp(n),bin_ind)=Speed_SumQxEyVy(flvp(n),bin_ind)+qp(n)*vp(n)*invg*pEy
       Speed_SumQxEzVz(flvp(n),bin_ind)=Speed_SumQxEzVz(flvp(n),bin_ind)+qp(n)*wp(n)*invg*pEz          
  end do     
  
  do n=1,used_test_prtl_arr_size
       if(qtp(n).eq.0) cycle

       SumNprtl(flvtp(n))=SumNprtl(flvtp(n))+1
       SumQ(flvtp(n))=SumQ(flvtp(n))+qtp(n)
       gamma=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
       invg=1.0_psn/gamma
       SumQxKE(flvtp(n))=SumQxKE(flvtp(n))+qtp(n)*(gamma-1) ! kinetic energy is in units of m_0c^2 (m_0=rest mass)
  
       call InterpEMfield(xtp(n),ytp(n),ztp(n),pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)
       SumQxExVx(flvtp(n))=SumQxExVx(flvtp(n))+qtp(n)*utp(n)*invg*pEx
       SumQxEyVy(flvtp(n))=SumQxEyVy(flvtp(n))+qtp(n)*vtp(n)*invg*pEy
       SumQxEzVz(flvtp(n))=SumQxEzVz(flvtp(n))+qtp(n)*wtp(n)*invg*pEz
       SumQxPx2(flvtp(n)) =SumQxPx2(flvtp(n)) +qtp(n)*utp(n)**2
       SumQxPy2(flvtp(n)) =SumQxPy2(flvtp(n)) +qtp(n)*vtp(n)**2
       SumQxPz2(flvtp(n)) =SumQxPz2(flvtp(n)) +qtp(n)*wtp(n)**2
  
       speed=sqrt((gamma-1.0_psn)*(gamma+1.0_psn))*invg
       bin_ind= floor(speed/Speed_spec_binwidth)+1
       Speed_SumQxExVx(flvtp(n),bin_ind)=Speed_SumQxExVx(flvtp(n),bin_ind)+qtp(n)*utp(n)*invg*pEx
       Speed_SumQxEyVy(flvtp(n),bin_ind)=Speed_SumQxEyVy(flvtp(n),bin_ind)+qtp(n)*vtp(n)*invg*pEy
       Speed_SumQxEzVz(flvtp(n),bin_ind)=Speed_SumQxEzVz(flvtp(n),bin_ind)+qtp(n)*wtp(n)*invg*pEz          
  end do  
             
end subroutine SumPrtlQntyMain 
subroutine CalcSpeedBinMeanEV(ArrSum,ArrMean)
  real(dbpsn), dimension(Nflvr,Speed_spec_binlen), intent(in):: ArrSum
  real, dimension(Nflvr,Speed_spec_binlen),intent(inout) :: ArrMean 
  integer :: i,j
  do j=1,Speed_spec_binlen
     do i=1,Nflvr
           if(spec_speed(i,j).eq.0) then 
                 ArrMean(i,j)=0
            else 
               ArrMean(i,j)=real(ArrSum(i,j)/spec_speed(i,j))                    
            end if 
     end do
end do   
end subroutine CalcSpeedBinMeanEV


!---------------------------------------------------------------------------------------------------------------------------------
!   subroutines to communicate data related to spectrum
!---------------------------------------------------------------------------------------------------------------------------------
       subroutine GetGmaxGminGlobal
		 integer :: mpi_err  
            call MPI_ALLREDUCE(local_gmax_allflv,gmax_allflv,1,mpi_psn,mpi_max,MPI_COMM_WORLD,mpi_err)
            call MPI_ALLREDUCE(local_gmin_allflv,gmin_allflv,1,mpi_psn,mpi_min,MPI_COMM_WORLD,mpi_err)
       end subroutine GetGmaxGminGlobal
  
       subroutine ReduceGammaSpectrum
		   integer :: mpi_err
            if(proc.eq.0) then 
              call MPI_REDUCE(MPI_IN_PLACE,spec_gamma,Nflvr*Gamma_spec_binlen,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
            else 
                call MPI_REDUCE(spec_gamma,  spec_gamma,Nflvr*Gamma_spec_binlen,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
            end if
       end subroutine ReduceGammaSpectrum
	   
       subroutine AllReduceGammaSpectrum
		   integer :: mpi_err
           call MPI_AllREDUCE(MPI_IN_PLACE,spec_gamma,Nflvr*Gamma_spec_binlen,MPI_REAL8,mpi_sum,MPI_COMM_WORLD,mpi_err)
       end subroutine AllReduceGammaSpectrum
            
       subroutine ReduceSpeedSpectrum
		   integer :: mpi_err
            if(proc.eq.0) then 
              call MPI_REDUCE(MPI_IN_PLACE,spec_speed,Nflvr*Speed_spec_binlen,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
            else 
                call MPI_REDUCE(spec_speed,  spec_speed,Nflvr*Speed_spec_binlen,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
            end if
       end subroutine ReduceSpeedSpectrum
       subroutine AllReduceSpeedSpectrum
		   integer :: mpi_err
           call MPI_AllREDUCE(MPI_IN_PLACE,spec_speed,Nflvr*Speed_spec_binlen,MPI_REAL8,mpi_sum,MPI_COMM_WORLD,mpi_err)
       end subroutine AllReduceSpeedSpectrum
	   
!----------------------------------------------------------------------------------------------------------------
! Following subroutines are used to communicate data related to "prtl_average" 
!----------------------------------------------------------------------------------------------------------------
subroutine ReduceSumPrtlQntyAll
    call ReduceReal8Arr1(sumNprtl,Nflvr)
    call ReduceReal8Arr1(SumQ,Nflvr)
    call ReduceReal8Arr1(SumQxKE,Nflvr)
    call ReduceReal8Arr1(SumQxExVx,Nflvr)
    call ReduceReal8Arr1(SumQxEyVy,Nflvr)
    call ReduceReal8Arr1(SumQxEzVz,Nflvr)
    call ReduceReal8Arr1(SumQxPx2,Nflvr)
    call ReduceReal8Arr1(SumQxPy2,Nflvr)
    call ReduceReal8Arr1(SumQxPz2,Nflvr)
    call ReduceReal8Arr2(Speed_SumQxExVx,Nflvr,Speed_spec_binlen)
    call ReduceReal8Arr2(Speed_SumQxEyVy,Nflvr,Speed_spec_binlen)
    call ReduceReal8Arr2(Speed_SumQxEzVz,Nflvr,Speed_spec_binlen)
    !call ReduceReal8Arr1(SumQxLR,Nflvr)

end subroutine ReduceSumPrtlQntyAll


!--------------------------------------------------------------------------------------------------------------
! Following subroutines are useful for reducing arrays on master proc. 
!--------------------------------------------------------------------------------------------------------------
subroutine ReduceReal8Arr1(Arr,sizex)
    integer :: sizex
    real(dbpsn), dimension(sizex) :: Arr
    integer :: mpi_err
    if(proc.eq.0) then 
       call MPI_REDUCE(MPI_IN_PLACE,Arr,sizex,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
    else 
        call MPI_REDUCE(Arr,  Arr,sizex,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
    end if
end subroutine ReduceReal8Arr1
subroutine ReduceReal8Arr2(Arr,sizex,sizey)
    integer :: sizex,sizey
    real(dbpsn), dimension(sizex,sizey) :: Arr
    integer :: mpi_err
    if(proc.eq.0) then 
       call MPI_REDUCE(MPI_IN_PLACE,Arr,sizex*sizey,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
    else 
        call MPI_REDUCE(Arr,  Arr,sizex*sizey,MPI_REAL8,mpi_sum,0,MPI_COMM_WORLD,mpi_err)
    end if
end subroutine ReduceReal8Arr2 
	   
	
end module prtl_stats