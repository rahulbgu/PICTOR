module savedata_gpu
	use parameters
	use vars
	use var_gpu 
#ifdef cyl
     use cyl_savedata_routines  
#endif 			
	use cudafor
    integer, device :: fdataxi_gpu,fdatayi_gpu,fdatazi_gpu,procx
	integer, device :: xb_min_gpu,yb_min_gpu,zb_min_gpu
	integer, device :: fsave_ratio_gpu
	integer, device :: tosave_prtl_size_gpu
	integer         :: tosave_prtl_size_host
	real(psn), dimension(:,:,:), allocatable, device :: fdata_gpu
	
contains 
	!inititalis  fdata nd borderes 
	 subroutine InitSaveDataGPU(xi,yi,zi,xb1,yb1,zb1) 
		 integer :: xi,yi,zi,xb1,yb1,zb1
		 fdataxi_gpu=xi 
		 fdatayi_gpu=yi
		 fdatazi_gpu=zi
		 xb_min_gpu=xb1
		 yb_min_gpu=yb1
		 zb_min_gpu=zb1 
		 fsave_ratio_gpu=fsave_ratio
		 procx=procxind(proc)
	 end subroutine InitSaveDataGPU 
	 
	 
	 subroutine CalcPrtlEnergyGPU(elc_energy,ion_energy)
		 real(psn) :: elc_energy, ion_energy
		 integer :: size = 1024
		 real , dimension(size)         :: elc_energy_wide , ion_energy_wide
		 real , dimension(size), device :: elc_energy_wide_gpu, ion_energy_wide_gpu
		 integer :: n, kc, indi,indf
		 
		  elc_energy_wide=0 ; ion_energy_wide =0;
		  elc_energy_wide_gpu = elc_energy_wide; ion_energy_wide_gpu = ion_energy_wide;
		  
	      do kc=1,Nchunk_prtl_gpu
	   		   indi=(kc-1)*chunk_size_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
			   call  CalcPrtlEnergyKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(up_gpu,vp_gpu,wp_gpu,qp_gpu,flvp_gpu,elc_energy_wide_gpu,ion_energy_wide_gpu,size,indi,indf)
		  end do 
		 

		  elc_energy_wide = elc_energy_wide_gpu; ion_energy_wide = ion_energy_wide_gpu;
		  do n=1,size
			  elc_energy = elc_energy + elc_energy_wide(n)
			  ion_energy = ion_energy + ion_energy_wide(n)
		  end do
		  
	 end subroutine CalcPrtlEnergyGPU
	 
	 attributes(global) subroutine CalcPrtlEnergyKernel(u,v,w,q,flv,elc_energy,ion_energy,size,indi,indf)
		 integer, value :: size,indi,indf
		 real, dimension(size) :: elc_energy, ion_energy
		 real, dimension(:) :: q,u,v,w
		 integer, dimension(:) :: flv
		 real    :: eng
		 integer :: n,m,stat
		 
	     n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
	     if(n.gt.indf) return 
	     
		 m=mod(n-1,size)+1
		 
         eng = (sqrt(1.0+u(n)**2+v(n)**2+w(n)**2) - 1.0) * abs(q(n))     
         if(flv(n).eq.1)  stat=atomicAdd( ion_energy(m) , eng ) 
         if(flv(n).eq.2)  stat=atomicAdd( elc_energy(m) , eng ) 
	 end subroutine CalcPrtlEnergyKernel
	 
     subroutine CalcPrtlDensityGPU(ch)
		 integer :: ch
		 integer :: kc,indi,indf
		 
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jx_gpu,mx,my,mz)			
	      do kc=1,Nchunk_prtl_gpu
	   		   indi=(kc-1)*chunk_size_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
               call CalcPrtlDensityGPUKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(ch,Jx_gpu,mx,my,mz,&
			   qp_gpu,xp_gpu,yp_gpu,zp_gpu,flvp_gpu,indi,indf,xb_min_gpu,yb_min_gpu,zb_min_gpu,fdataxi_gpu,fdatayi_gpu,fdatazi_gpu,fsave_ratio_gpu,procx)
          end do

		  call NormaliseFldDensity1_GPU

     end subroutine CalcPrtlDensityGPU
	 
     subroutine CalcTestPrtlDensityGPU(ch)
		 integer :: ch
		 integer :: kc,indi,indf
	 

		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jx_gpu,mx,my,mz)			
	      do kc=1,Nchunk_test_prtl_gpu
	   		   indi=(kc-1)*chunk_size_test_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
               call CalcPrtlDensityGPUKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(ch,Jx_gpu,mx,my,mz,&
			   qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,flvtp_gpu,indi,indf,xb_min_gpu,yb_min_gpu,zb_min_gpu,fdataxi_gpu,fdatayi_gpu,fdatazi_gpu,fsave_ratio_gpu,procx)
          end do

		  call NormaliseFldDensity1_GPU

     end subroutine CalcTestPrtlDensityGPU
	 
     subroutine CalcPrtlChargeFluxGPU(ch)
		 integer :: ch
		 integer :: kc,indi,indf
	 
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jx_gpu,mx,my,mz)	
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jy_gpu,mx,my,mz)			
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jz_gpu,mx,my,mz)			
		  		
	      do kc=1,Nchunk_prtl_gpu
	   		   indi=(kc-1)*chunk_size_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
               call CalcPrtlChargeFluxGPUKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(ch,Jx_gpu,Jy_gpu,Jz_gpu,mx,my,mz,&
			   qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,flvp_gpu,indi,indf,xb_min_gpu,yb_min_gpu,zb_min_gpu,fdataxi_gpu,fdatayi_gpu,fdatazi_gpu,fsave_ratio_gpu,procx)
          end do

		  call NormaliseFldDensity3_GPU		  
     end subroutine CalcPrtlChargeFluxGPU
	 
     subroutine CalcTestPrtlChargeFluxGPU(ch)
		 integer :: ch
		 integer :: kc,indi,indf	 

		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jx_gpu,mx,my,mz)	
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jy_gpu,mx,my,mz)			
		  call ResetMatrixGPUKernelSaveData<<<grid,tBlock>>>(Jz_gpu,mx,my,mz)			
		  		
	      do kc=1,Nchunk_test_prtl_gpu
	   		   indi=(kc-1)*chunk_size_test_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
               call CalcPrtlChargeFluxGPUKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(ch,Jx_gpu,Jy_gpu,Jz_gpu,mx,my,mz,&
			   qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,flvtp_gpu,indi,indf,xb_min_gpu,yb_min_gpu,zb_min_gpu,fdataxi_gpu,fdatayi_gpu,fdatazi_gpu,fsave_ratio_gpu,procx)
          end do

		  call NormaliseFldDensity3_GPU
     end subroutine CalcTestPrtlChargeFluxGPU
	 
	 !-------------------------------------------------------------------------------------------------------------------------
	 !      Normalisation
	 !-------------------------------------------------------------------------------------------------------------------------
	 subroutine NormaliseFldDensity1_GPU
		 real :: cell_vol

#ifdef twoD		 
		 cell_vol=(fsave_ratio**2)
#else
         cell_vol=(fsave_ratio**3)
#endif

#ifndef cyl		  
		  call NormaliseFldDensity1GPUKernel<<<grid,tBlock>>>(Jx_gpu,mx,my,mz,cell_vol)
		  Jx=Jx_gpu
#endif

#ifdef cyl		  
          Jx=Jx_gpu
		  call NormaliseFldDensity1
#endif		 
		 
	 end subroutine NormaliseFldDensity1_GPU
	 
	 
	 subroutine NormaliseFldDensity3_GPU
		 real :: cell_vol

#ifdef twoD		 
		 cell_vol=(fsave_ratio**2)
#else
         cell_vol=(fsave_ratio**3)
#endif

#ifndef cyl		 
		  call NormaliseFldDensity1GPUKernel<<<grid,tBlock>>>(Jx_gpu,mx,my,mz,cell_vol)
		  call NormaliseFldDensity1GPUKernel<<<grid,tBlock>>>(Jy_gpu,mx,my,mz,cell_vol)
		  call NormaliseFldDensity1GPUKernel<<<grid,tBlock>>>(Jz_gpu,mx,my,mz,cell_vol)
          Jx=Jx_gpu
		  Jy=Jy_gpu
		  Jz=Jz_gpu
#endif		  	  

#ifdef cyl
          Jx=Jx_gpu
		  Jy=Jy_gpu
		  Jz=Jz_gpu
		  call NormaliseFldDensity3
#endif		 

	 end subroutine NormaliseFldDensity3_GPU
	 
	 attributes(global) subroutine CalcPrtlDensityGPUKernel(ch,Fld,mx,my,mz,q,x,y,z,flv,indi,indf,xb,yb,zb,fdataxi,fdatayi,fdatazi,fsave_ratio,procx)
	       integer, value :: ch,mx,my,mz
		   integer :: procx	
           real  , dimension(mx,my,mz) :: Fld 		           
 	       real, dimension(:) :: q,x,y,z
 	       integer, dimension(:)  :: flv
		   integer, value :: indi,indf
 		   
		   integer  :: xb,yb,zb
 		   integer  :: fdataxi,fdatayi,fdatazi 
		   integer  :: fsave_ratio
		   
   		   integer :: n 
   		   integer :: stat
		   integer :: i,j,k,ip,jp,kp
           real    :: Wx,Wy,Wz,Wxp,Wyp,Wzp
#ifdef twoD
           real,dimension(4) :: Den
#else 
           real,dimension(8) :: Den
#endif 	
		   
		   
		   n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		   if(n.gt.indf) return 
		   
               if(flv(n).eq.ch) then
                    call  DownsampleGridIndexGPU(x(n),y(n),z(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,fsave_ratio,xb,yb,zb,fdataxi,fdatayi,fdatazi,mx,my,mz,procx)
                    Den(1)= Wx *Wy *Wz*q(n)
                    Den(2)= Wxp*Wy *Wz*q(n)
                    Den(3)= Wx *Wyp*Wz*q(n)
                    Den(4)= Wxp*Wyp*Wz*q(n)
#ifndef twoD
                    Den(5)= Wx *Wy *Wzp*q(n)
                    Den(6)= Wxp*Wy *Wzp*q(n)
                    Den(7)= Wx *Wyp*Wzp*q(n)
                    Den(8)= Wxp*Wyp*Wzp*q(n)
#endif


                    stat=atomicAdd(Fld(i  ,j  ,k  ), Den(1))
                    stat=atomicAdd(Fld(ip ,j  ,k  ), Den(2))
                    stat=atomicAdd(Fld(i  ,jp ,k  ), Den(3))
                    stat=atomicAdd(Fld(ip ,jp ,k  ), Den(4))
#ifndef twoD
                    stat=atomicAdd(Fld(i  ,j  ,kp ), Den(5))
                    stat=atomicAdd(Fld(ip ,j  ,kp ), Den(6))
                    stat=atomicAdd(Fld(i  ,jp ,kp ), Den(7))
                    stat=atomicAdd(Fld(ip ,jp ,kp ), Den(8))
#endif

             end if 
	 end subroutine CalcPrtlDensityGPUKernel 
	 
	 
	 attributes(global) subroutine CalcPrtlChargeFluxGPUKernel(ch,Fldx,Fldy,Fldz,mx,my,mz,q,x,y,z,u,v,w,flv,indi,indf,xb,yb,zb,fdataxi,fdatayi,fdatazi,fsave_ratio,procx)
	       integer, value :: ch, mx,my,mz
		   integer :: x1,x2,y1,y2,z1,z2,procx
           real  , dimension(mx,my,mz) :: Fldx,Fldy,Fldz 	
		   real, dimension(:) :: q,x,y,z,u,v,w
		   integer, dimension(:)  :: flv
		   integer, value :: indi,indf	           
		   integer  :: xb,yb,zb
 		   integer  :: fdataxi,fdatayi,fdatazi 
		   integer  :: fsave_ratio
		   
   		   integer :: n 
   		   integer :: stat
		   integer :: i,j,k,ip,jp,kp
           real    :: Wx,Wy,Wz,Wxp,Wyp,Wzp
           real    ::vx,vy,vz,invg
#ifdef twoD
           real,dimension(4) :: Velp 
#else 
           real,dimension(8) :: Velp
#endif 	

		   n = blockDim%x * (blockIdx%x - 1) + threadIdx%x +(indi-1)
		   if(n.gt.indf) return 
		   
           if(flv(n).eq.ch) then
                call  DownsampleGridIndexGPU(x(n),y(n),z(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,fsave_ratio,xb,yb,zb,fdataxi,fdatayi,fdatazi,mx,my,mz,procx)
           
                    invg=1/sqrt(1+u(n)**2+v(n)**2+w(n)**2)
                    vx=u(n)*invg
                    vy=v(n)*invg
                    vz=w(n)*invg


                    Velp(1)= Wx *Wy *Wz*vx*q(n)
                    Velp(2)= Wxp*Wy *Wz*vx*q(n)
                    Velp(3)= Wx *Wyp*Wz*vx*q(n)
                    Velp(4)= Wxp*Wyp*Wz*vx*q(n)
                    stat=atomicAdd(Fldx(i  ,j  ,k  ), Velp(1))
                    stat=atomicAdd(Fldx(ip ,j  ,k  ), Velp(2))  
                    stat=atomicAdd(Fldx(i  ,jp ,k  ), Velp(3))
                    stat=atomicAdd(Fldx(ip ,jp ,k  ), Velp(4))
#ifndef twoD
					Velp(5)= Wx *Wy *Wzp*vx*q(n)
					Velp(6)= Wxp*Wy *Wzp*vx*q(n)
					Velp(7)= Wx *Wyp*Wzp*vx*q(n)
					Velp(8)= Wxp*Wyp*Wzp*vx*q(n)
                    stat=atomicAdd(Fldx(i  ,j  ,kp ), Velp(5))
                    stat=atomicAdd(Fldx(ip ,j  ,kp ), Velp(6))
                    stat=atomicAdd(Fldx(i  ,jp ,kp ), Velp(7))
                    stat=atomicAdd(Fldx(ip ,jp ,kp ), Velp(8))
#endif

					Velp(1)= Wx *Wy *Wz*vy*q(n)
					Velp(2)= Wxp*Wy *Wz*vy*q(n)
					Velp(3)= Wx *Wyp*Wz*vy*q(n)
					Velp(4)= Wxp*Wyp*Wz*vy*q(n)
					stat=atomicAdd(Fldy(i  ,j  ,k  ), Velp(1))
					stat=atomicAdd(Fldy(ip ,j  ,k  ), Velp(2))  
					stat=atomicAdd(Fldy(i  ,jp ,k  ), Velp(3))
					stat=atomicAdd(Fldy(ip ,jp ,k  ), Velp(4))
#ifndef twoD
					Velp(5)= Wx *Wy *Wzp*vy*q(n)
					Velp(6)= Wxp*Wy *Wzp*vy*q(n)
					Velp(7)= Wx *Wyp*Wzp*vy*q(n)
					Velp(8)= Wxp*Wyp*Wzp*vy*q(n)
					stat=atomicAdd(Fldy(i  ,j  ,kp ), Velp(5))
					stat=atomicAdd(Fldy(ip ,j  ,kp ), Velp(6))
					stat=atomicAdd(Fldy(i  ,jp ,kp ), Velp(7))
					stat=atomicAdd(Fldy(ip ,jp ,kp ), Velp(8))
#endif
					Velp(1)= Wx *Wy *Wz*vz*q(n)
					Velp(2)= Wxp*Wy *Wz*vz*q(n)
					Velp(3)= Wx *Wyp*Wz*vz*q(n)
					Velp(4)= Wxp*Wyp*Wz*vz*q(n)
					stat=atomicAdd(Fldz(i  ,j  ,k  ), Velp(1))
					stat=atomicAdd(Fldz(ip ,j  ,k  ), Velp(2))  
					stat=atomicAdd(Fldz(i  ,jp ,k  ), Velp(3))
					stat=atomicAdd(Fldz(ip ,jp ,k  ), Velp(4))
#ifndef twoD
					Velp(5)= Wx *Wy *Wzp*vz*q(n)
					Velp(6)= Wxp*Wy *Wzp*vz*q(n)
					Velp(7)= Wx *Wyp*Wzp*vz*q(n)
					Velp(8)= Wxp*Wyp*Wzp*vz*q(n)
					stat=atomicAdd(Fldz(i  ,j  ,kp ), Velp(5))
					stat=atomicAdd(Fldz(ip ,j  ,kp ), Velp(6))
					stat=atomicAdd(Fldz(i  ,jp ,kp ), Velp(7))
					stat=atomicAdd(Fldz(ip ,jp ,kp ), Velp(8))
#endif
		   end if
	 end subroutine CalcPrtlChargeFluxGPUKernel
	 
	 
     attributes(device) subroutine DownsampleGridIndexGPU(x,y,z,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,fsave_ratio,xb,yb,zb,fdataxi,fdatayi,fdatazi,mx,my,mz,procx)
	      integer :: procx 
		  real :: x,y,z
          integer :: i,j,k,ip,jp,kp
          real :: Wx,Wy,Wz,Wxp,Wyp,Wzp
          real :: fp
		  integer  :: fsave_ratio
		  integer  :: xb,yb,zb
		  integer  :: fdataxi,fdatayi,fdatazi 
		  integer  :: mx,my,mz
		  real :: binlen
		  
		       binlen=real(fsave_ratio)/2 
               i=floor2realGPU(real(x-3+xb-fdataxi),fsave_ratio)+fdataxi-xb+3
               j=floor2realGPU(real(y-3+yb-fdatayi),fsave_ratio)+fdatayi-yb+3
               
               ip=i+fsave_ratio
               jp=j+fsave_ratio
               fp=x-i-binlen+0.5
               Wxp=max(min(fp,1.0),0.0)
               Wx =max(min(1-fp,1.0),0.0)
               fp =y-j-binlen+0.5
               Wyp=max(min(fp,1.0),0.0)
               Wy =max(min(1-fp,1.0),0.0)
               
               i=max(1,i)
               j=max(1,j)
               ip=min(mx,ip)
               jp=min(my,jp)

#ifndef twoD
               k=floor2realGPU(real(z-3+zb-fdatazi),fsave_ratio)+fdatazi-zb+3
               kp=k+fsave_ratio
               fp =z-k-binlen+0.5
               Wzp=max(min(fp,1.0),0.0)
               Wz =max(min(1-fp,1.0),0.0)     
               k=max(1,k)
               kp=min(mz,kp)
#else
               Wz=1
               k=1
#endif
#ifdef cyl
               call AxisPrtlWtGPU(x,Wx,Wxp,procx,fsave_ratio)	 
#endif
     end subroutine DownsampleGridIndexGPU
	 
     integer attributes(device) function floor2realGPU(n1,n2)
          real    :: n1
          integer :: n2
          if(n1.ge.0) then 
               floor2realGPU=n2*int(n1/n2)
          else 
               floor2realGPU=n2*int(n1/n2)-n2
          end if
     end function floor2realGPU
	 
  	attributes(global) subroutine NormaliseFldDensity1GPUKernel(Fld,mx,my,mz,norm)
  		integer, value :: mx,my,mz		
  		real(psn), dimension(mx,my,mz) :: Fld
        real, value :: norm 
  		integer :: i,j,k

		
  		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
  		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
 #ifndef twoD 		
  		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
 #else
        k=1
 #endif  
          if((i.le.mx).and.(j.le.my).and.(k.le.mz)) Fld(i,j,k)=Fld(i,j,k)/norm
  	end subroutine NormaliseFldDensity1GPUKernel
	
 	attributes(global) subroutine ResetMatrixGPUKernelSaveData(Fld,mx,my,mz)
 		integer, value :: mx,my,mz
 		real(psn), dimension(mx,my,mz) :: Fld
  		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif  
         if((i.le.mx).and.(j.le.my).and.(k.le.mz)) Fld(i,j,k)=0.0
 	end subroutine ResetMatrixGPUKernelSaveData
	
!-----------------------------------------------------------------------	
! Particles subroutines
!-----------------------------------------------------------------------
subroutine CollectPrtlGPU(vid)
     integer :: vid
	 !Note that prtl send recv buffer arrays are used for saving prtl data, and therefore only a limited fraction of prtl data can be saved     
                 select case (vid)
                 case(1) 
	 			   pdata_real(1:tosave_prtl_size_host)=xp_send_gpu(1:tosave_prtl_size_host)					   
                 case(2)
 			       pdata_real(1:tosave_prtl_size_host)=yp_send_gpu(1:tosave_prtl_size_host)					   
                 case(3)
		           pdata_real(1:tosave_prtl_size_host)=zp_send_gpu(1:tosave_prtl_size_host)					   
                 case(4)
	               pdata_real(1:tosave_prtl_size_host)=up_send_gpu(1:tosave_prtl_size_host)					   
                 case(5)
                   pdata_real(1:tosave_prtl_size_host)=vp_send_gpu(1:tosave_prtl_size_host)					   
                 case(6)
                   pdata_real(1:tosave_prtl_size_host)=wp_send_gpu(1:tosave_prtl_size_host)					   
                 case(7)
                   pdata_real(1:tosave_prtl_size_host)=qp_send_gpu(1:tosave_prtl_size_host)					   
                 case(8)    
                   pdata_int(1:tosave_prtl_size_host)=tagp_send_gpu(1:tosave_prtl_size_host)					   
                 case(9)
                   pdata_int(1:tosave_prtl_size_host)=flvp_send_gpu(1:tosave_prtl_size_host)					   
			     case(10)
                   pdata_real(1:tosave_prtl_size_host)=var1p_send_gpu(1:tosave_prtl_size_host)					   
                 case(11)
 			       pdata_real(1:tosave_prtl_size_host)=xp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(12)
		           pdata_real(1:tosave_prtl_size_host)=yp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(13)
	               pdata_real(1:tosave_prtl_size_host)=zp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(14)
                   pdata_real(1:tosave_prtl_size_host)=up_recv_gpu(1:tosave_prtl_size_host)					   
                 case(15)
                   pdata_real(1:tosave_prtl_size_host)=vp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(16)
                   pdata_real(1:tosave_prtl_size_host)=wp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(17)
                   pdata_real(1:tosave_prtl_size_host)=up_recv_gpu(1:tosave_prtl_size_host)					   
                 case(18)
                   pdata_real(1:tosave_prtl_size_host)=vp_recv_gpu(1:tosave_prtl_size_host)					   
                 case(19)
                   pdata_real(1:tosave_prtl_size_host)=wp_recv_gpu(1:tosave_prtl_size_host)					   
                 end select      
end subroutine CollectPrtlGPU
	
     subroutine CalcPrtlGPU
		 integer :: kc,indi,indf
		 tosave_prtl_size_host=0
         tosave_prtl_size_gpu=tosave_prtl_size_host
	      do kc=1,Nchunk_prtl_gpu
	   		   indi=(kc-1)*chunk_size_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_prtl_gpu+used_prtl_chunk(kc)
               call CollectPrtlGPUKernel<<<ceiling(real(used_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qp_gpu,xp_gpu,yp_gpu,zp_gpu,up_gpu,vp_gpu,wp_gpu,var1p_gpu,tagp_gpu,flvp_gpu,indi,indf,&
               qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu,tagp_send_gpu,flvp_send_gpu,tosave_prtl_size_gpu)
          end do
		  !test particles
	      do kc=1,Nchunk_test_prtl_gpu
	   		   indi=(kc-1)*chunk_size_test_prtl_gpu+1
	   		   indf=(kc-1)*chunk_size_test_prtl_gpu+used_test_prtl_chunk(kc)
               call CollectPrtlGPUKernel<<<ceiling(real(used_test_prtl_chunk(kc))/NthreadsGPU), NthreadsGPU>>>(qtp_gpu,xtp_gpu,ytp_gpu,ztp_gpu,utp_gpu,vtp_gpu,wtp_gpu,var1tp_gpu,tagtp_gpu,flvtp_gpu,indi,indf,&
               qp_send_gpu,xp_send_gpu,yp_send_gpu,zp_send_gpu,up_send_gpu,vp_send_gpu,wp_send_gpu,var1p_send_gpu,tagp_send_gpu,flvp_send_gpu,tosave_prtl_size_gpu)
          end do
          tosave_prtl_size_host=tosave_prtl_size_gpu		  	
  		  tosave_prtl_arr_size=tosave_prtl_size_host  
     end subroutine CalcPrtlGPU
	 
	 subroutine CalcPrtlLocalFieldGPU
	          !Note that on gpu the data is stored in recv buffer 
              call CollectPrtlLocalFldGPUKernel<<<ceiling(real(tosave_prtl_size_host)/NthreadsGPU), NthreadsGPU>>>(TexEx_gpu,TexEy_gpu,TexEz_gpu,TexBx_gpu,TexBy_gpu,TexBz_gpu,mx,my,mz,&
			  xp_send_gpu,yp_send_gpu,zp_send_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,up_recv_gpu,vp_recv_gpu,wp_recv_gpu,tosave_prtl_size_gpu)
	 end subroutine CalcPrtlLocalFieldGPU
	 
	 subroutine CalcPrtlLocalCurrGPU
              call CollectPrtlLocalCurrGPUKernel<<<ceiling(real(tosave_prtl_size_host)/NthreadsGPU), NthreadsGPU>>>(Jx_gpu,Jy_gpu,Jz_gpu,mx,my,mz,&
			  xp_send_gpu,yp_send_gpu,zp_send_gpu,xp_recv_gpu,yp_recv_gpu,zp_recv_gpu,tosave_prtl_size_gpu)
	 end subroutine CalcPrtlLocalCurrGPU
	 
	 
	 attributes(global) subroutine CollectPrtlGPUKernel(q1,x1,y1,z1,u1,v1,w1,var1,tag1,flv1,i1,i2,q2,x2,y2,z2,u2,v2,w2,var2,tag2,flv2,index)
	     real, dimension(:) ::q1,x1,y1,z1,u1,v1,w1,var1
		 integer, dimension(:) :: flv1,tag1  
	     real, dimension(:) ::q2,x2,y2,z2,u2,v2,w2,var2
		 integer, dimension(:) :: flv2,tag2  
		 integer, value :: i1,i2
		 integer :: index
		 integer :: InsertAt
		 integer :: n
		  
		 n = blockDim%x * (blockIdx%x - 1) + threadIdx%x + (i1-1) 
		 if(n.gt.i2) return 
		 if(tag1(n).ne.0) then 
			 InsertAt=atomicinc(index,1000000000)
		 
			 InsertAt=InsertAt+1
								  
			 q2(InsertAt)=q1(n)		 
		     x2(InsertAt)=x1(n)
			 y2(InsertAt)=y1(n)
			 z2(InsertAt)=z1(n)
			 u2(InsertAt)=u1(n)
			 v2(InsertAt)=v1(n)
			 w2(InsertAt)=w1(n)
			 var2(InsertAt)=var1(n)
			 flv2(InsertAt)=flv1(n)
			 tag2(InsertAt)=tag1(n)
		 end if 
	 end subroutine CollectPrtlGPUKernel
	 
	 attributes(global) subroutine CollectPrtlLocalFldGPUKernel(Ex,Ey,Ez,Bx,By,Bz,mx,my,mz,x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,Nprtl)
	    integer, value :: mx,my,mz
	    real, dimension(mx,my,mz)  :: Ex,Ey,Ez,Bx,By,Bz
        real, dimension(:) ::x,y,z
        real, dimension(:) ::pEx,pEy,pEz,pBx,pBy,pBz
        
		integer :: Nprtl 
        integer :: n

		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: wt 
#else 
        real,dimension(8) :: wt
#endif 	
        n = blockDim%x * (blockIdx%x - 1) + threadIdx%x  
		if(n.le.Nprtl) then 
		  i=x(n)
		  j=y(n)
		  dx=x(n)-i
		  dy=y(n)-j
#ifndef twoD
          k=z(n)
          dz=z(n)-k
	      wt(1)=(1.0-dx)*(1.0-dy)*(1.0-dz)
	      wt(2)=dx          *(1.0-dy)*(1.0-dz)
	      wt(3)=(1.0-dx) *dy          *(1.0-dz)
	      wt(4)=dx          *dy          *(1.0-dz)
	      wt(5)=(1.0-dx)*(1.0-dy)*dz
	      wt(6)=dx          *(1.0-dy)*dz
	      wt(7)=(1.0-dx)*dy          *dz
	      wt(8)=dx          *dy          *dz
#else
          k=1
		  wt(1)=(1.0-dx)*(1.0-dy)
		  wt(2)=dx          *(1.0-dy)
		  wt(3)=(1.0-dx) *dy
		  wt(4)=dx          *dy
#endif


#ifdef twoD
          pEx(n)=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)
		  pEy(n)=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)
		  pEz(n)=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)
          pBx(n)=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)
		  pBy(n)=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)
		  pBz(n)=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)
#else
          pEx(n)=wt(1)*Ex(i,j,k)+wt(2)*Ex(i+1,j,k)+wt(3)*Ex(i,j+1,k)+wt(4)*Ex(i+1,j+1,k)+wt(5)*Ex(i,j,k+1)+wt(6)*Ex(i+1,j,k+1)+wt(7)*Ex(i,j+1,k+1)+wt(8)*Ex(i+1,j+1,k+1)
          pEy(n)=wt(1)*Ey(i,j,k)+wt(2)*Ey(i+1,j,k)+wt(3)*Ey(i,j+1,k)+wt(4)*Ey(i+1,j+1,k)+wt(5)*Ey(i,j,k+1)+wt(6)*Ey(i+1,j,k+1)+wt(7)*Ey(i,j+1,k+1)+wt(8)*Ey(i+1,j+1,k+1)
          pEz(n)=wt(1)*Ez(i,j,k)+wt(2)*Ez(i+1,j,k)+wt(3)*Ez(i,j+1,k)+wt(4)*Ez(i+1,j+1,k)+wt(5)*Ez(i,j,k+1)+wt(6)*Ez(i+1,j,k+1)+wt(7)*Ez(i,j+1,k+1)+wt(8)*Ez(i+1,j+1,k+1)
          pBx(n)=wt(1)*Bx(i,j,k)+wt(2)*Bx(i+1,j,k)+wt(3)*Bx(i,j+1,k)+wt(4)*Bx(i+1,j+1,k)+wt(5)*Bx(i,j,k+1)+wt(6)*Bx(i+1,j,k+1)+wt(7)*Bx(i,j+1,k+1)+wt(8)*Bx(i+1,j+1,k+1)
          pBy(n)=wt(1)*By(i,j,k)+wt(2)*By(i+1,j,k)+wt(3)*By(i,j+1,k)+wt(4)*By(i+1,j+1,k)+wt(5)*By(i,j,k+1)+wt(6)*By(i+1,j,k+1)+wt(7)*By(i,j+1,k+1)+wt(8)*By(i+1,j+1,k+1)
          pBz(n)=wt(1)*Bz(i,j,k)+wt(2)*Bz(i+1,j,k)+wt(3)*Bz(i,j+1,k)+wt(4)*Bz(i+1,j+1,k)+wt(5)*Bz(i,j,k+1)+wt(6)*Bz(i+1,j,k+1)+wt(7)*Bz(i,j+1,k+1)+wt(8)*Bz(i+1,j+1,k+1)	  
#endif
         end if 
     end subroutine CollectPrtlLocalFldGPUKernel
	 
	 attributes(global) subroutine CollectPrtlLocalCurrGPUKernel(Jx,Jy,Jz,mx,my,mz,x,y,z,pJx,pJy,pJz,Nprtl)
	    integer, value :: mx,my,mz
	    real, dimension(mx,my,mz)  :: Jx,Jy,Jz
        real, dimension(:) ::x,y,z
		real, dimension(:) ::pJx,pJy,pJz
        
		integer :: Nprtl 
        integer :: n

		integer :: i,j,k
		real :: dx,dy,dz
#ifdef twoD
        real,dimension(4) :: wt,FldVec 
#else 
        real,dimension(8) :: wt,FldVec
#endif 	
        n = blockDim%x * (blockIdx%x - 1) + threadIdx%x  
		if(n.le.Nprtl) then
						
  		  i=x(n)
  		  j=y(n)
  		  dx=x(n)-i
  		  dy=y(n)-j
#ifndef twoD
          k=z(n)
          dz=z(n)-k
  	      wt(1)=(1.0-dx)*(1.0-dy)*(1.0-dz)
  	      wt(2)=dx          *(1.0-dy)*(1.0-dz)
  	      wt(3)=(1.0-dx) *dy          *(1.0-dz)
  	      wt(4)=dx          *dy          *(1.0-dz)
  	      wt(5)=(1.0-dx)*(1.0-dy)*dz
  	      wt(6)=dx          *(1.0-dy)*dz
  	      wt(7)=(1.0-dx)*dy          *dz
  	      wt(8)=dx          *dy          *dz
#else
          k=1
  		  wt(1)=(1.0-dx)*(1.0-dy)
  		  wt(2)=dx          *(1.0-dy)
  		  wt(3)=(1.0-dx) *dy
  		  wt(4)=dx          *dy
#endif

		  FldVec(1)=Jx(i-1,j,k)+Jx(i,j,k)
		  FldVec(2)=Jx(i,j,k)+Jx(i+1,j,k)
		  FldVec(3)=Jx(i-1,j+1,k)+Jx(i,j+1,k)
		  FldVec(4)=Jx(i,j+1,k)+Jx(i+1,j+1,k)
#ifndef twoD
		  FldVec(5)=Jx(i-1,j,k+1)+Jx(i,j,k+1)
		  FldVec(6)=Jx(i,j,k+1)+Jx(i+1,j,k+1)
		  FldVec(7)=Jx(i-1,j+1,k+1)+Jx(i,j+1,k+1)
		  FldVec(8)=Jx(i,j+1,k+1)+Jx(i+1,j+1,k+1)
#endif

#ifdef twoD
          pJx(n)=0.5*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#else
          pJx(n)=0.5*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#endif

		  FldVec(1)=Jy(i,j-1,k)+Jy(i,j,k)
		  FldVec(2)=Jy(i+1,j-1,k)+Jy(i+1,j,k)
		  FldVec(3)=Jy(i,j,k)+Jy(i,j+1,k)
		  FldVec(4)=Jy(i+1,j,k)+Jy(i+1,j+1,k)
#ifndef twoD
		  FldVec(5)=Jy(i,j-1,k+1)+Jy(i,j,k+1)
		  FldVec(6)=Jy(i+1,j-1,k+1)+Jy(i+1,j,k+1)
		  FldVec(7)=Jy(i,j,k+1)+Jy(i,j+1,k+1)
		  FldVec(8)=Jy(i+1,j,k+1)+Jy(i+1,j+1,k+1)
#endif

#ifdef twoD
          pJy(n)=0.5*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#else
          pJy(n)=0.5*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#endif

#ifndef twoD
		  FldVec(1)=Jz(i,j,k-1)+Jz(i,j,k)
		  FldVec(2)=Jz(i+1,j,k-1)+Jz(i+1,j,k)
		  FldVec(3)=Jz(i,j+1,k-1)+Jz(i,j+1,k)
		  FldVec(4)=Jz(i+1,j+1,k-1)+Jz(i+1,j+1,k)
		  FldVec(5)=Jz(i,j,k)+Jz(i,j,k+1)
		  FldVec(6)=Jz(i+1,j,k)+Jz(i+1,j,k+1)
		  FldVec(7)=Jz(i,j+1,k)+Jz(i,j+1,k+1)
		  FldVec(8)=Jz(i+1,j+1,k)+Jz(i+1,j+1,k+1)
	      pJz(n)=0.5*(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4)+wt(5)*FldVec(5)+wt(6)*FldVec(6)+wt(7)*FldVec(7)+wt(8)*FldVec(8))
#else
		  FldVec(1)=Jz(i,j,k)
		  FldVec(2)=Jz(i+1,j,k)
		  FldVec(3)=Jz(i,j+1,k)
		  FldVec(4)=Jz(i+1,j+1,k)
          pJz(n)=(wt(1)*FldVec(1)+wt(2)*FldVec(2)+wt(3)*FldVec(3)+wt(4)*FldVec(4))
#endif
		end if
	
	end subroutine CollectPrtlLocalCurrGPUKernel
		
	
	
end module savedata_gpu