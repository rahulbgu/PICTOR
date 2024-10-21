module cyl_help_setup
	use parameters 
	use vars 
	use cyl_vars
	use prob_dist
	use help_setup
    use communication
    use mem_prtl
	use em_update
	implicit none
	
	real(psn) :: ri_local,rf_local,yi_local,yf_local,zi_local,zf_local
	real(psn) :: rmin_global2, rmax_global2
contains 
	
!---------------------------------------------------------------------------------
!  subroutines to initialise particles 
!---------------------------------------------------------------------------------		
	
    subroutine InitPrtl_cyl(Flvr,Density,Temperature,DriftVelocity,MaxDen)
		integer :: Flvr
		real(psn), external :: Density,Temperature
		external :: DriftVelocity
		real(psn), optional :: MaxDen
		real(psn)           :: max_den
		
		max_den=1.0_psn
		if(present(MaxDen)) max_den=MaxDen
		call InitPrtlPos_cyl(FlvrCharge(Flvr),Flvr,Density,max_den)
		call InitPrtlMom_cyl(Flvr,Temperature,DriftVelocity)
		call GetUsedPrtlIndex(initialised_prtl_ind)
	end subroutine InitPrtl_cyl
	
	subroutine InitPrtlPair_cyl(Flvr1,Flvr2,Density,Temperature1,Temperature2,DriftVelocity1,DriftVelocity2,MaxDen)
		integer :: Flvr1,Flvr2
		real(psn), external :: Density,Temperature1,Temperature2
		external :: DriftVelocity1,DriftVelocity2
		real(psn), optional :: MaxDen
		real(psn)           :: max_den
		max_den=1.0_psn
		if(present(MaxDen)) max_den=MaxDen
		call InitPrtlPos_cyl(FlvrCharge(Flvr1),Flvr1,Density,max_den)
		call CopyPrtlPos(Flvr2,Flvr1)
		call InitPrtlMom_cyl(Flvr1,Temperature1,DriftVelocity1)
		call InitPrtlMom_cyl(Flvr2,Temperature2,DriftVelocity2)
		call GetUsedPrtlIndex(initialised_prtl_ind)
	end subroutine InitPrtlPair_cyl 
	
!---------------------------------------------------------------------------------
!  subroutines to initialise Electric and Magnetic field
!---------------------------------------------------------------------------------		
    subroutine InitElectricField_cyl(Efld)
		external :: Efld
		integer :: i,j,k, i1
		real(psn) :: r,theta,z
		real(psn) :: e_x,e_y,e_z
	    i1=1
	    if(procxind.eq.0) i1=4 
#ifdef twoD
        do k=1,1
#else 	
	    do k=1,mz
#endif 
        do j=1,my 
		    do i=i1,mx
				r=i -3.0_psn + rborders(procxind)
				theta= 2*pi*(j -3.0_psn + yborders(procyind))/ny
#ifdef twoD    
                z=0.0_psn
#else				 				
				z= k - 3.0 + zborders(proczind)
#endif 				
				call Efld(r+0.5_psn,theta,z,e_x,e_y,e_z)
				Ex(i,j,k)=e_x
			    call Efld(r,theta+dtheta*0.5_psn,z,e_x,e_y,e_z)
				Ey(i,j,k)=e_y
#ifdef twoD				
				call Efld(r,theta,0.0_psn,e_x,e_y,e_z)
#else				
				call Efld(r,theta,z+0.5_psn,e_x,e_y,e_z)
#endif				
				
				Ez(i,j,k)=e_z
		    end do 
        end do  
        end do
		
		if(procxind.eq.0) call UpdateEFldAxis 
! 		if(procxind.eq.0) then
! 			Ex(1:3,:,:)=0.0_psn
! 			Ey(1:3,:,:)=0.0_psn
! 			Ez(1:3,:,:)=0.0_psn
! 		end if
	end subroutine InitElectricField_cyl
    subroutine InitMagneticField_cyl(Bfld)
		external :: Bfld
		integer :: i,j,k,i1
		real(psn) :: r,theta,z,zp_half
		real(psn) :: b_x,b_y,b_z,b_x1,b_x2
	    i1=1
	    if(procxind.eq.0) i1=4 
#ifdef twoD
        do k=1,1
#else 	
	    do k=1,mz
#endif 
        do j=1,my 
			
			theta= 2*pi*(j -3.0_psn + yborders(procyind))/ny
#ifdef twoD    
            z=0.0_psn
		    zp_half=0.0_psn
#else				 				
			z= k - 3.0 + zborders(proczind)
			zp_half=z+0.5_psn
#endif 	
		    do i=i1,mx
				r=i -3.0_psn + rborders(procxind)
			
				call Bfld(r,theta+0.5_psn*dtheta,zp_half,b_x,b_y,b_z)
				Bx(i,j,k)=b_x
			    call Bfld(r+0.5_psn,theta,zp_half,b_x,b_y,b_z)
				By(i,j,k)=b_y				
				call Bfld(r+0.5_psn,theta+dtheta*0.5_psn,z,b_x,b_y,b_z)	
				Bz(i,j,k)=b_z
		    end do 
			
			!By and Bz at the axis
			if(procxind.eq.0) then 
				r=0
				call Bfld(r+0.0_psn,theta+(pi/2.0_psn),z,b_x1,b_y,b_z)	
				call Bfld(r+0.0_psn,theta-(pi/2.0_psn),z,b_x2,b_y,b_z)	
				By(3,j,k) = 0.5_psn*(b_x1 - b_x2)
				!if nan then compute the limit r->0 
				call Bfld(r,theta+dtheta*0.5_psn,z,b_x,b_y,b_z)	
				Bz(3,j,k) =  b_z
			end if

			
        end do  
        end do 
		
		if(procxind.eq.0) call UpdateBFldAxis 


	end subroutine InitMagneticField_cyl
		
	
	subroutine InitPrtlMom_cyl(Flvr,Temperature,DriftVelocity)
		integer :: Flvr
		real(psn), external :: Temperature
		external :: DriftVelocity
		real(psn) :: TempThis,r,theta,z
		real(psn) :: ugamma,vgamma,wgamma
		real(psn)  :: vx,vy,vz
		integer :: i
		
		
		do i=initialised_prtl_ind+1,prtl_arr_size
			if((qp(i).ne.0).and.flvp(i).eq.Flvr) then  
				r=xp(i)-xmin+rborders(procxind)
				theta=(yp(i)-ymin+yborders(procyind))*dtheta
				z=zp(i)-zmin+zborders(proczind)
				TempThis=Temperature(r,theta,z)
				call GetVelGamma_MaxBolt(TempThis,ugamma,vgamma,wgamma)
				call DriftVelocity(r,theta,z,vx,vy,vz)
				call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
				up(i)=ugamma
				vp(i)=vgamma
				wp(i)=wgamma   
				
		    end if 
		end do 
	end subroutine InitPrtlMom_cyl


	!The method of estimating prtl arr size is now changed, change it to the method used in the rectangular version
	!the new version first estiamtes prtl arr size and then allocates the prtl arr and then intialise the particles 
 	subroutine InitPrtlPos_cyl(ch,flv,Den,max_den)
		integer :: flv
		real(psn) :: ch,max_den
		real(psn), external :: Den
		real(psn) :: r_global2,rmin_global2,rmax_global2,xlocal,ylocal,zlocal,yglobal,zglobal,theta
		integer   :: yglobal_max
		real(psn) :: rnd1,rnd2,rnd3,rnd_acpt
		integer   :: i,j,off,off_prev,init_block_this,estmt_size,Nelc_uniform_cyl_this
		integer   :: init_block=1000000
	    real(psn) :: rmin,rmax
		
	    rmin_global2=rborders(procxind)**2
	    rmax_global2=rborders(procxind+1)**2
		if(procxind.eq.0) rmin_global2=0
		
		yglobal_max = yborders(procyind+1)
		if(inc_axis) then
			if(procxind.eq.0) yglobal_max=ny
		end if
		
	    rmin=rborders(procxind)
	    rmax=rborders(procxind+1)
	    if(procxind.eq.0) rmin=0
		Nelc_uniform_cyl_this=0.5_psn*epc*dtheta*(rmax**2-rmin**2)*(ymax-ymin)*(zmax-zmin)*max_den
		
		if(Nelc_uniform_cyl_this.lt.1000) call GetIntPoissonDist(real(real(Nelc_uniform_cyl_this),psn),Nelc_uniform_cyl_this)
		call GetUsedPrtlIndex(off)
		if(off+init_block.gt.prtl_arr_size) call ReshapePrtlArr(new_size=int(prtl_arr_size+1.1*init_block),used_ind=off)
	
		do i=0,Nelc_uniform_cyl_this-1,init_block
			init_block_this=min(init_block,Nelc_uniform_cyl_this-i)
		    off_prev=off
			do j=1,init_block_this
				call random_number(rnd1)
			    call random_number(rnd2)
			    call random_number(rnd3)
				call random_number(rnd_acpt)
				r_global2 = rmin_global2 + rnd1*(rmax_global2 - rmin_global2)
				yglobal= yborders(procyind) + rnd2*(yglobal_max-yborders(procyind))
				zglobal= zborders(proczind) + rnd3*(zborders(proczind+1)-zborders(proczind))
				theta=yglobal*dtheta
				if(rnd_acpt.le.Den(sqrt(r_global2),theta,zglobal)/max_den) then
					 xlocal=sqrt(r_global2)-rborders(procxind)+xmin
		 			 ylocal=yglobal - yborders(procyind) + ymin
		 			 zlocal=zglobal - zborders(proczind) + zmin
					 call InsertParticleAt(off+1,xlocal,ylocal,zlocal,0.0_psn,0.0_psn,0.0_psn,ch,0,flv,0.0_psn)
					 off=off+1
				end if 
		    end do

			if(off+init_block.gt.prtl_arr_size) then 
				estmt_size=prtl_arr_size+int(max(1.1*init_block, 1.1*init_block* ((Nelc_uniform_cyl_this-i)/init_block-1)*(real(off-off_prev)/real(init_block))))
				call ReshapePrtlArr(new_size=estmt_size, used_ind=off)
			end if 

		end do 
	end subroutine InitPrtlPos_cyl
	
	
	
end module cyl_help_setup