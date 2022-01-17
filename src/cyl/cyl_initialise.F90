module cyl_initialise
	use parameters
	use vars 
	use cyl_vars 
	use communication
	use memory
	implicit none
contains 

	subroutine InitAll_cyl
		call InitFldVars_cyl
	end subroutine InitAll_cyl
	
	subroutine InitParam_cyl
		dtheta=2.0_psn*pi/ny
		ax_perm_area=4.0_psn/ny
		inv_dtheta=1.0_psn/dtheta
		Nelc_uniform_cyl=0
	end subroutine InitParam_cyl
#ifdef gpu
    subroutine InitAll_cyl_gpu 
		allocate(Bz_ax_gpu(mz),Ey_ax_gpu(my,mz),Ey_ax_host(my,mz))
        tBlock_ax = dim3 (1 ,64 ,4)
        grid_ax = dim3(1, ceiling(real(my)/tBlock_ax%y), ceiling(real(mz)/tBlock_ax%z)) 
	end subroutine InitAll_cyl_gpu 
#endif 	
	subroutine InitPrtlArrSize_cyl
		 real(psn) :: rmin,rmax
		 
	     BC_Rmin_Prtl=-1; BC_Rmax_Prtl=-1;
	     BC_Rmin_Fld=-1; BC_Rmax_Fld=-1;
         
		 rmin=rborders(procxind(proc))
		 rmax=rborders(procxind(proc)+1)
		 if(procxind(proc).eq.0) rmin=0

		 
		 Nelc_uniform_cyl=0.5_psn*epc*dtheta*(rmax**2-rmin**2)*(ymax-ymin)*(zmax-zmin)
		 prtl_arr_size=prtl_arr_size0!min(int(3.0*Nelc_uniform_cyl),prtl_arr_size0)
		 
#ifdef gpu
         prtl_arr_size=gpu_prtl_arr_size 
#endif 
		
	end subroutine InitPrtlArrSize_cyl

	subroutine InitPrtlTransferInOutArr_cyl 
		if(inc_axis) then
			if(procxind(proc).eq.0) then 
			    allocate(rinp_count_axis(0:nSubDomainsY-1),rintp_count_axis(0:nSubDomainsY-1),rpcross_axis(0:nSubDomainsY-1),rtpcross_axis(0:nSubDomainsY-1))
			    rinp_count_axis=0;rintp_count_axis=0;rpcross_axis=0;rtpcross_axis=0;
			    allocate(rinp_axis(nSubDomainsY*rinp_size,0:nSubDomainsY-1),routp_axis(nSubDomainsY*routp_size,0:nSubDomainsY-1))
		    end if 
		end if 	
	end subroutine InitPrtlTransferInOutArr_cyl
    
	subroutine InitFldVars_cyl
		if(procxind(proc).eq.0) then 
			allocate(Bz_ax(mz))
		    Bz_ax=0
		end if 
	end subroutine InitFldVars_cyl
	
	subroutine Initboundaries_cyl
		integer :: i

		if(inc_axis) then
		    xborders(0)=0
			xborders(1)=4
			do i=2,nSubDomainsX
				xborders(i)=xborders(i-1)+fsave_ratio*( ((nx-4)/(nSubDomainsX-1))/fsave_ratio)
			end do
			xborders(nSubDomainsX)=nx
		end if
		
		call SetRborders

	end subroutine Initboundaries_cyl
	subroutine SetRborders
		rborders=xborders-0.5_psn
	end subroutine SetRborders
	
	subroutine UpdateDomainSkelton_cyl				 
		 if(inc_axis.eqv..false.) then 
			call CheckNum_MPI_Task_cyl 
	 		if(procxind(proc).eq.0) lproc=MPI_PROC_NULL
	 		if(procxind(proc).eq.nSubDomainsX-1) rproc=MPI_PROC_NULL
		 end if 
		 if(inc_axis.eqv..true.) then 
			call CheckNum_MPI_Task_IncAxis
			call IncAxisProcGrid
			call IncAxisUpdateNeighbor
		 end if   		
	end subroutine UpdateDomainSkelton_cyl
	
	subroutine IncAxisProcGrid
        integer :: i,j,k, m
#ifdef twoD
        do k=0,0
#else                 
        do k=0,nSubDomainsZ-1
#endif               
              do j=0,nSubDomainsY-1
                  do i=1,nSubDomainsX-1
                        proc_grid(i,j,k)=(i-1)+j*(nSubDomainsX-1)+k*(nSubDomainsX-1)*nSubDomainsY
                  end do
              end do
         end do	
#ifdef twoD
        m= (nSubDomainsX-1)*nSubDomainsY 
        do k=0,0
#else                 
        m= (nSubDomainsX-1)*nSubDomainsY*nSubDomainsZ
		do k=0,nSubDomainsZ-1
#endif	 
             proc_grid(0,:,k)=k + m 
		end do 
	end subroutine IncAxisProcGrid
	
	subroutine IncAxisDomainSize
		if((inc_axis.eqv..true.).and.procxind(proc).eq.0) then
            my=ny+5
			ymax=my-2
			ylen=my-5.0_psn
	    end if 
	end subroutine IncAxisDomainSize
	
	subroutine IncAxisUpdateNeighbor
		integer :: i,j,k
		integer :: ip,im,jp,jm,kp,km
#ifdef twoD
        do k=0,0
#else                 
        do k=0,nSubDomainsZ-1
#endif               
              do j=0,nSubDomainsY-1
                  do i=0,nSubDomainsX-1
					  procxind(proc_grid(i,j,k))=i
					  procyind(proc_grid(i,j,k))=j
					  proczind(proc_grid(i,j,k))=k
					  if(proc.eq.proc_grid(i,j,k)) then 
						  ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km= k-1;
						  if(ip.gt.nSubDomainsX-1) ip=0
						  if(jp.gt.nSubDomainsY-1) jp=0
						  if(im.lt.0) im=nSubDomainsX-1
						  if(jm.lt.0) jm=nSubDomainsY-1						
						  lproc=proc_grid(im,j,k)
						  rproc=proc_grid(ip,j,k)
						  tproc=proc_grid(i,jp,k)
						  bproc=proc_grid(i,jm,k)
#ifndef twoD
                          if(kp.gt.nSubDomainsZ-1) kp=0 
						  if(km.lt.0) km=nSubDomainsZ-1						  
						  uproc=proc_grid(i,j,kp)
						  dproc=proc_grid(i,j,km)
#endif 
					  end if 
                  end do
              end do
         end do
		 
#ifdef twoD
			 do k=0,0
#else                 
			 do k=0,nSubDomainsZ-1
#endif 			 
			      procyind(proc_grid(0,0,k))=0 
		     end do
		 
		 if(procxind(proc).eq.1) lproc=MPI_PROC_NULL
 		 if(procxind(proc).eq.nSubDomainsX-1) rproc=MPI_PROC_NULL
   		 if(procxind(proc).eq.0) then 
 			lproc=MPI_PROC_NULL
 			rproc=MPI_PROC_NULL
 		 end if 		 			
	end subroutine IncAxisUpdateNeighbor
	
			
    subroutine CheckNum_MPI_Task_cyl 
#ifdef twoD          
          if(nproc.ne.nSubDomainsX*nSubDomainsY) then
                    STOP 'Error :: Number of MPI tasks does not match (nSubDomainsX*nSubDomainsY) the number declared in parameter.F90 file'
          end if
#else 
          if(nproc.ne.nSubDomainsX*nSubDomainsY*nSubDomainsZ) then
                    STOP 'Error :: Number of MPI tasks does not match (nSubDomainsX*nSubDomainsY*nSubDomainsZ) the number declared in parameter.F90 file'
          end if
#endif     
     end subroutine CheckNum_MPI_Task_cyl	
	 
     subroutine CheckNum_MPI_Task_IncAxis 
#ifdef twoD          
           if(nproc.ne.((nSubDomainsX-1)*nSubDomainsY+1)) then
                     STOP 'Error :: the number of MPI tasks must be (nSubDomainsX-1)*nSubDomainsY + 1, because the axis lies entirly within one subdomain'
           end if
#else 
           if(nproc.ne.((nSubDomainsX-1)*nSubDomainsY*nSubDomainsZ+nSubDomainsZ)) then
                      STOP 'Error :: the number of MPI tasks must be (nSubDomainsX-1)*nSubDomainsY*nSubDomainsZ + nSubDomainsZ, because the axis lies entirly within one subdomain'
           end if
#endif     
      end subroutine CheckNum_MPI_Task_IncAxis		
	  
	  
	  subroutine ComplyBC_cyl
		  integer :: n 
		  real(psn) :: rmax_global,rmin_global,zmin_global,zmax_global
		  real(psn) :: rmax_local,rmin_local,zmin_local,zmax_local
		  
		  rmax_global=BC_Rmax_Prtl
		  rmin_global=BC_Rmin_Prtl
		  zmax_global=BC_Zmax_Prtl
		  zmin_global=BC_Zmin_Prtl
		  
		  if(BC_Rmax_Prtl.lt.0) rmax_global=rborders(nSubDomainsX)+1
		  if(BC_Rmin_Prtl.lt.0) rmin_global=rborders(0)-1
		  if(inc_axis) rmin_global=-1 
		  rmax_local=rmax_global-rborders(procxind(proc))+3
		  rmin_local=rmin_global-rborders(procxind(proc))+3
		  		  
#ifdef twoD		  
          zmax_local=2
          zmin_local=-1
#else 	  
		  if(BC_Zmax_Prtl.lt.0) zmax_global=zborders(nSubDomainsZ)+1
		  if(BC_Zmin_Prtl.lt.0) zmin_global=zborders(0)-1
		  zmax_local=zmax_global-zborders(proczind(proc))+3
		  zmin_local=zmin_global-zborders(proczind(proc))+3
#endif 		  

		  


          do n=1,prtl_arr_size
               if(qp(n).ne.0) then
				   if((xp(n).lt.rmin_local).or.(xp(n).gt.rmax_local)) call DeletePrtl(n)
				   if((zp(n).lt.zmin_local).or.(zp(n).gt.zmax_local)) call DeletePrtl(n)
			   end if
		  end do 
		  
          do n=1,test_prtl_arr_size
               if(qtp(n).ne.0) then
				   if((xtp(n).lt.rmin_local).or.(xtp(n).gt.rmax_local)) call DeleteTestPrtl(n)
				   if((ztp(n).lt.zmin_local).or.(ztp(n).gt.zmax_local)) call DeleteTestPrtl(n)
			   end if
		  end do 
		  
	  end subroutine ComplyBC_cyl
	

end module cyl_initialise