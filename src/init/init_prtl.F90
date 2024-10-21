module init_prtl
	use parameters
	use vars
	use mem_prtl
	use prtl_tag
	use prob_dist
#ifdef cyl
    use cyl_common
#endif	
contains 
	
	 !-------------------------------------------------------------------------------------------------
	 ! Total number of particles in the subodmain, assuming uniform distribution and "epc" particles per cell
	 !-------------------------------------------------------------------------------------------------	 
	 
	 integer function Nelc_uniform()
		  Nelc_uniform=epc*(mx-5)*(my-5)      
#ifndef twoD
          Nelc_uniform=Nelc_uniform*(mz-5) ! Initial Number of electrons  
#endif	 
#ifdef cyl
          Nelc_uniform = Nelc_uniform_cyl(real(xborders(procxind),dbpsn),real(xborders(procxind+1),dbpsn),real(yborders(procyind),dbpsn),real(yborders(procyind+1),dbpsn),real(zborders(proczind),dbpsn),real(zborders(proczind+1),dbpsn)) 
#endif
	 end function Nelc_uniform
	 

integer function EstimatePrtlCount(Den,Nuniform)
	integer :: n, count , Nuniform
	procedure(scalar_global) :: Den
	real(dbpsn) :: r1,r2,r3,rnd_acpt
	real(dbpsn) :: xglobal,yglobal,zglobal

	count = 0

	do n=1,10000
		call random_number(r1)
	    call random_number(r2)
		call random_number(r3)
		call random_number(rnd_acpt)
	
		xglobal= xborders(procxind) + r1*(xborders(procxind+1)-xborders(procxind))
		yglobal= yborders(procyind) + r2*(yborders(procyind+1)-yborders(procyind))
		zglobal= zborders(proczind) + r3*(zborders(proczind+1)-zborders(proczind))
		
#ifdef cyl
		call random_xglobal_cyl(r1,xglobal,real(xborders(procxind),dbpsn),real(xborders(procxind+1),dbpsn)) 
#endif
	
		if(rnd_acpt.le.Den(xglobal,yglobal,zglobal)) count = count + 1 
	
	end do 

    EstimatePrtlCount = int(1.1*Nuniform*(count/10000.0)) + 1000000
	
end function EstimatePrtlCount

integer function prtl_count_each_placement(pid)
	integer, dimension(:) :: pid
	integer :: n 
	prtl_count_each_placement = 0
	do n=1,size(pid)
		prtl_count_each_placement = prtl_count_each_placement + PSP_list(pid(n))%multiplicity
	end do

end function prtl_count_each_placement


subroutine SetQbyM(ind,value) ! Update this soubroutine such that the Flvrs can be added in any order, no incremently 
     real(psn), dimension(:), allocatable :: flvrqmTemp,FlvrChargeTemp
     integer, dimension(:), allocatable   :: FlvrSaveFldDataTemp, FlvrTypeTemp,FlvrSaveRatioTemp,CurrentTagIDTemp, CurrentTagProcIDTemp
	 type(FlvrProperty), dimension(:), allocatable :: flvr_prpt_temp
	 integer :: ind, n
     real(psn) :: value 
     if(ind.gt.Nflvr) then 
          allocate(flvrqmTemp(ind),FlvrChargeTemp(ind),FlvrSaveFldDataTemp(ind),FlvrTypeTemp(ind),FlvrSaveRatioTemp(ind),CurrentTagIDTemp(ind), CurrentTagProcIDTemp(ind))
          flvrqmTemp(1:Nflvr)=flvrqm(1:Nflvr)
	      FlvrChargeTemp(1:NFlvr)=FlvrCharge(1:Nflvr)
          FlvrSaveFldDataTemp(1:Nflvr)=FlvrSaveFldData(1:Nflvr)     
          FlvrTypeTemp(1:Nflvr)=FlvrType(1:Nflvr)
          FlvrSaveRatioTemp(1:Nflvr)=FlvrSaveRatio(1:Nflvr)
          CurrentTagIDTemp(1:Nflvr)=CurrentTagID(1:Nflvr)
		  CurrentTagProcIDTemp(1:Nflvr)=CurrentTagProcID(1:Nflvr)

          deallocate(flvrqm,FlvrCharge,FlvrSaveFldData,FlvrType,FlvrSaveRatio,CurrentTagID, CurrentTagProcID)
          call move_alloc(flvrqmTemp,flvrqm)
	      call move_alloc(FlvrChargeTemp,FlvrCharge)
          call move_alloc(FlvrSaveFldDataTemp,FlvrSaveFldData)
          call move_alloc(FlvrTypeTemp,FlvrType)
          call move_alloc(FlvrSaveRatioTemp,FlvrSaveRatio)
          call move_alloc(CurrentTagIDTemp,CurrentTagID)
		  call move_alloc(CurrentTagProcIDTemp,CurrentTagProcID)

		  allocate(flvr_prpt_temp(ind))
		  do n=1,Nflvr
			flvr_prpt_temp(n) = flvr_prpt(n)
		  end do
		  call move_alloc(flvr_prpt_temp,flvr_prpt) 


          Nflvr=Nflvr+1
      end if
          flvrqm(ind)=value
end subroutine SetQbyM


subroutine SetPhaseSpaceProperty(psp,Flvr,Density,Temperature,SpeedDist,Vmax,DriftVelocity,Multiplicity, Fraction)
	type(PhaseSpaceProperty) :: psp
	integer :: Flvr
	procedure(scalar_global), optional :: Density
	procedure(scalar_global), optional :: Temperature
	procedure(func1D),        optional :: SpeedDist
	procedure(vector_global), optional :: DriftVelocity  
	integer                 , optional :: Multiplicity
	real(psn), optional :: Vmax !the maximum particle speed in the plasma frame, default is c 
	real(psn), optional :: Fraction

	psp%Flvr = Flvr
	if(present(Density))       psp%Density => Density 
	if(present(Temperature))   psp%Temperature => Temperature 
	if(present(DriftVelocity)) psp%DriftVelocity => DriftVelocity  

	psp%Vmax = 1.0_psn
	psp%multiplicity =1 
	if(present(Vmax)) psp%Vmax = Vmax  
	if(present(Multiplicity)) psp%multiplicity = Multiplicity
	if(present(Fraction)) psp%fraction = Fraction


	if(present(SpeedDist)) then 
		psp%SpeedDist => SpeedDist
		allocate( psp%Table(PDF_TableSize), psp%PDF_Table(PDF_TableSize) )
		call InitPDFTable(PDF_TableSize, psp%Table, psp%PDF_Table, psp%SpeedDist, psp%Vmax)
	end if 

end subroutine SetPhaseSpaceProperty 


subroutine InsertPrtl_PSP(pid,xglobal,yglobal,zglobal,xlocal,ylocal,zlocal, weight)
	integer, dimension(:) :: pid
	real(dbpsn) :: xglobal, yglobal, zglobal
	real(psn) :: xlocal, ylocal, zlocal, ugamma, vgamma, wgamma , vdx, vdy, vdz, weight
	real(psn) :: Temp
	real(psn) :: charge
	integer :: tag_id, proc_id, i, n
 
	do i=1,size(pid)
		do n = 1, PSP_list(pid(i))%multiplicity
			 if( associated( PSP_list(pid(i))%Temperature)) then 
				 Temp = PSP_list(pid(i))%Temperature(xglobal,yglobal,zglobal)
				 call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
			 else if( associated( PSP_list(pid(i))%SpeedDist)) then 
				 call GetIsoVelGammaTable(PDF_TableSize, PSP_list(pid(i))%Table, PSP_list(pid(i))%PDF_Table,ugamma,vgamma,wgamma) 
		     end if

			 if( associated( PSP_list(pid(i))%DriftVelocity)) then 
				 call PSP_list(pid(i))%DriftVelocity(xglobal,yglobal,zglobal,vdx,vdy,vdz)
				 call AddDriftVel(ugamma,vgamma,wgamma,vdx,vdy,vdz)
			 end if   

			 call GetTag( PSP_list(pid(i))%Flvr, tag_id, proc_id)
			 charge = weight * FlvrCharge( PSP_list(pid(i))%Flvr) * PSP_list(pid(i))%fraction
			 call InsertParticleAt(used_prtl_arr_size+1,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma, charge, PSP_list(pid(i))%Flvr, 0.0_psn , tag_id, proc_id, flvrqm(PSP_list(pid(i))%Flvr) ,weight)
			 used_prtl_arr_size=used_prtl_arr_size+1
			 np=np+1
	 	
		end do 
    end do 

end subroutine InsertPrtl_PSP



!-----------------------------------------------------------------------------------------------------------------------
! Determine position of the particles and then call routines to initialise momenta and insert them into the main array
!-----------------------------------------------------------------------------------------------------------------------

	subroutine InitPrtlFromPSP_RandomPosition(pid, Density,x1,x2,y1,y2,z1,z2)
		integer, dimension(:) :: pid
		procedure(scalar_global) :: Density
	    real(dbpsn) :: x1,x2,y1,y2,z1,z2 !range of domain where particle are placed; global cordinate
		real(dbpsn) :: r1,r2,r3,rnd_acpt
		real(dbpsn) :: xglobal, yglobal, zglobal
		real(psn) :: xlocal, ylocal, zlocal
		integer   :: n, est_np, Nprtl
	
		!make sure that the particle array is large enough 
		est_np = prtl_count_each_placement(pid)*EstimatePrtlCount(Density, Nelc_uniform())
		if(used_prtl_arr_size + est_np .gt. prtl_arr_size) call ReshapePrtlArr( new_size=int(prtl_arr_size + est_np) )
	
	    Nprtl = epc*max(0.0_psn,x2-x1)*max(0.0_psn,y2-y1)*max(0.0_psn,z2-z1)
#ifdef cyl
        Nprtl = Nelc_uniform_cyl(x1,x2,y1,y2,z1,z2)
#endif		
		
		do n=1,Nprtl
			call random_number(r1)
		    call random_number(r2)
			call random_number(r3)
			call random_number(rnd_acpt)
	
			xglobal= x1 + r1*(x2-x1)
			yglobal= y1 + r2*(y2-y1)
			zglobal= z1 + r3*(z2-z1)

#ifdef cyl
			call random_xglobal_cyl(r1,xglobal,x1,x2) 
#endif			
				
			if(rnd_acpt.le.Density(xglobal,yglobal,zglobal)) then
			
			    xlocal= xglobal - xborders(procxind) + xmin
			    ylocal= yglobal - yborders(procyind) + ymin
			    zlocal= zglobal - zborders(proczind) + zmin

				call InsertPrtl_PSP(pid,xglobal,yglobal,zglobal, xlocal,ylocal,zlocal, 1.0_psn ) 
	 
			end if 
		end do 
	
	
	end subroutine InitPrtlFromPSP_RandomPosition


	subroutine InitPrtlFromPSP_UniformSeparation(pid, Density, x1_g,x2_g,y1_g,y2_g,z1_g,z2_g, max_displacement)
		integer, dimension(:) :: pid
		procedure(scalar_global) :: Density
		integer :: i,j,k, est_np
		real(dbpsn) :: x1_g,x2_g,y1_g,y2_g,z1_g,z2_g !range of domain where particle are placed; global cordinate
		real(psn) :: x1,x2,y1,y2,z1,z2 !range of domain where particle are placed; local cordinate
		integer   :: i1,i2,j1,j2,k1,k2, si,sj,sk, si_max, sj_max, sk_max 
		real(dbpsn) :: xglobal, yglobal, zglobal
		real(psn)   :: xlocal, ylocal, zlocal
		real(psn)   :: r1, r2, r3, r11
		real(psn) :: dpos_x, dpos_y, dpos_z, wt
		real(psn) :: max_displacement
	
		x1 = x1_g - xborders(procxind)+3; x2 = x2_g -xborders(procxind)+3;
		y1 = y1_g - yborders(procyind)+3; y2 = y2_g -yborders(procyind)+3;
		z1 = z1_g - zborders(proczind)+3; z2 = z2_g -zborders(proczind)+3;
#ifdef twoD
        z1 = zmin; z2 = zmax;
#endif		
		
	
	    
#ifdef twoD		
		if(epc_x*epc_y.ne.epc) then 
 			 STOP 'Err : change epc in the parameter file such that epc = epc_x*epc_y '
		end if
	
#else		
		if(epc_x*epc_y*epc_z.ne.epc) then 
			STOP 'Err: change epc in the parameter file such that epc = epc_x*epc_y*epc_z'
		end if 
#endif		
	
        dpos_x = 1.0_dbpsn/real(epc_x,dbpsn)
		dpos_y = 1.0_dbpsn/real(epc_y,dbpsn)
		dpos_z = 1.0_dbpsn/real(epc_z,dbpsn)

	
		i1 = floor(x1)
		i2 = ceiling(x2)
		j1 = floor(y1)
		j2 = ceiling(y2)
		k1 = floor(z1)
		k2 = ceiling(z2)

		!make sure that the particle array is large enough 
		est_np = epc*(i2-i1+1)*(j2-j1+1)*(k2-k1+1)*prtl_count_each_placement(pid)
		if(used_prtl_arr_size + est_np .gt. prtl_arr_size) call ReshapePrtlArr( new_size=int(prtl_arr_size + est_np) )
		
		
		si_max = epc_x; sj_max=epc_y; sk_max=epc_z;
#ifdef cyl
		if(RZtwoD) then
			j2 = j1 
			sj_max = 1
		end if 
#endif

		
	
#ifdef twoD
        do k=1,1
#else
        do k=k1,k2
#endif 		
			do j=j1,j2
				do i=i1,i2
					
					
					!placement within the cell at i,j,k
					do si=1,si_max
						do sj=1,sj_max
#ifndef twoD 							
							do sk=1,sk_max 
#else	
                            do sk=1,1 							
#endif								
								xlocal = i+ (si-1)*dpos_x 
								ylocal = j+ (sj-1)*dpos_y
								zlocal = k+ (sk-1)*dpos_z 
								
#ifdef cyl
								if(RZtwoD) then
									call random_number(r11)
									ylocal = r11*ny +3.0_psn
								end if
#endif									
								
								xglobal= xlocal + xborders(procxind) -3
								yglobal= ylocal + yborders(procyind) -3
								zglobal= zlocal + zborders(proczind) -3
#ifdef twoD
								zlocal = 1.0_psn
								zglobal = 0.0_psn
#endif					
			                    if( xlocal.ge.x1 .and. xlocal.lt.x2 .and. ylocal.ge.y1 .and. ylocal.lt.y2 .and. zlocal.ge.z1 .and. zlocal.lt.z2 ) then 
					
									wt = Density(xglobal,yglobal,zglobal)
#ifdef cyl
			                        call azimuthal_weight_prtl(wt,xglobal) !mutiply wt by r\dtheta
#endif						
					
									if(max_displacement.gt.0) then ! add a random displacement vector to particles position
										call random_number(r1)
										call random_number(r2)
										xlocal = xlocal + (2.0_psn*r1-1.0_psn)*max_displacement*dpos_x
										ylocal = ylocal + (2.0_psn*r2-1.0_psn)*max_displacement*dpos_y
#ifndef twoD 							
										call random_number(r3)
										zlocal = zlocal + (2.0_psn*r3-1.0_psn)*max_displacement*dpos_z
#endif							
					
									end if 
					
									if(wt.gt.0) call InsertPrtl_PSP(pid,xglobal,yglobal,zglobal, xlocal,ylocal,zlocal, wt ) 
				
								end if 
							
							end do 							
						end do 
					end do 
					
				
				end do 
			end do 
		end do
			
	
	end subroutine InitPrtlFromPSP_UniformSeparation

end module init_prtl