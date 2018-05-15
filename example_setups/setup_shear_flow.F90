module setup
     use parameters
     use vars
	 use help_setup
	 implicit none
	 
integer, parameter :: GammaTableLength=10000
real(psn), dimension(GammaTableLength) :: GammaElc,GammaIon,Gamma_PDF_Elc,Gamma_PDF_Ion
integer :: xcond_left,xcond_right,PrtlWallLeft,PrtlWallRight
real(psn) :: Temp_ion,Temp_elc
real(psn) :: Alfven_speed,PlasmaBeta,Btheta,Bphi,BextMag,TeTiRatio,AmpVd,ShearScale

 contains 

subroutine InitUser
	real(psn) :: xlocal,ylocal,zlocal,ubeta,vbeta,wbeta,ugamma,vgamma,wgamma,x2,x1
	real(psn) :: r1,r2,r3,beta_drift
	integer :: i
    !Physical parameters 
	Alfven_speed=0.1_psn
	PlasmaBeta=0.1
	Btheta=(pi/180)*2.0
	Bphi=0.0
	AmpVd=0.5
	ShearScale=3*c_ompe*sqrt(mi/me) 
	TeTiRatio=1
	
	!boundary condition related parameters
    xcond_left=4
	xcond_right=nx-xcond_left
	PrtlWallLeft=32
	PrtlWallRight=nx-PrtlWallLeft 
	
    !set the background magentic field 
    BextMag=Alfven_speed*sqrt(epc*(massi+masse)*c*c)     
    Bz_ext0=BextMag*cos(Btheta)
    Bx_ext0=BextMag*sin(Btheta)*cos(Bphi)
	By_ext0=BextMag*sin(Btheta)*sin(Bphi)
	
    Temp_ion=PlasmaBeta*(Alfven_speed**2)*0.5 !in units of rest mass energy, Temp_ion=kT_i/m_ic^2 
    Temp_elc=Temp_ion*TeTiRatio*(mi/me)! Temp_elc=kT_e/m_ec^2
    
	call InitMaxwellNonRelTable(GammaTableLength,GammaElc,Gamma_PDF_Elc,Temp_elc)
    call InitMaxwellNonRelTable(GammaTableLength,GammaIon,Gamma_PDF_Ion,Temp_ion)

	x1=max(3,PrtlWallLeft-xborders(procxind(proc))+3)
	x2=min(mx-2,PrtlWallRight-xborders(procxind(proc))+3)
#ifndef twoD
          Nelc=epc*(x2-x1)*(my-5)*(mz-5) ! Initial Number of electrons  
#else
          Nelc=epc*(x2-x1)*(my-5)
#endif    
	
	do i=1,Nelc
        !call GetRandomPositionHomogeneous(xlocal,ylocal,zlocal)
       call random_number(r1)
       call random_number(r2)
       call random_number(r3)

       xlocal=x1+(x2-x1)*r1
       ylocal=ymin+(ymax-ymin)*r2
       zlocal=zmin+(zmax-zmin)*r3
	   call GetDriftVelocity(xlocal,ylocal,zlocal,ubeta,vbeta,wbeta)   
       beta_drift=sqrt(ubeta**2+vbeta**2+wbeta**2)	      
       call GetVelGamma_MaxwellJuttner(beta_drift,Temp_ion,ugamma,vgamma,wgamma,direction=(/ubeta,vbeta,wbeta/))

	   call InsertParticleAt(i,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0,0,1,0.0_psn)
   end do
   
   !electrons 
  		  
   do i=1,Nelc
	    call GetDriftVelocity(xp(i),yp(i),zp(i),ubeta,vbeta,wbeta)      
	    beta_drift=sqrt(ubeta**2+vbeta**2+wbeta**2)
	    call GetVelGamma_MaxwellJuttner(beta_drift,Temp_elc,ugamma,vgamma,wgamma,direction=(/ubeta,vbeta,wbeta/))  
        call InsertParticleAt(i+Nelc,xp(i),yp(i),zp(i),ugamma,vgamma,wgamma,-1.0,0,2,0.0_psn) 
   end do   
      
   call FldCondBoundary 
end subroutine InitUser


subroutine GetDriftVelocity(xlocal,ylocal,zlocal,ubeta,vbeta,wbeta)
     real(psn) :: ubeta,vbeta,wbeta
     real(psn) :: xlocal,ylocal,zlocal
     real(dbpsn)::xglobal,yglobal,zglobal
     real(psn) :: upol_drift,vpol_drift,wpol_drift
      
    xglobal=xlocal-3+xborders(procxind(proc))
    yglobal=ylocal-3+yborders(procyind(proc))
#ifdef twoD     
    zglobal=zlocal
#else
    zglobal=zlocal-3+zborders(proczind(proc))
#endif
     
     !EXB drift 
	 vbeta=AmpVd*tanh((xglobal-nx/2)/ShearScale)  
	 ubeta=0.0
	 wbeta=0.0         
end subroutine GetDriftVelocity


subroutine InitEMfield
     integer :: i,j,k
     real(dbpsn):: xglobal,yglobal,zglobal
     real(psn)  :: Ax,Ay,Az
     do k=1,mz
          do j=1,my
               do i=1,mx                                        
                    call GetGlobalPosition_DP(real(i),real(j+0.5),real(k),xglobal,yglobal,zglobal)
                    Ex(i,j,k)=Ex(i,j,k)-AmpVd*Bz_ext0*tanh((xglobal-nx/2)/ShearScale)                                  
               end do
          end do
     end do
end subroutine InitEMfield

subroutine RefPrtl
	integer :: PrtlWall_local
	integer :: n
	
	PrtlWall_local=PrtlWallLeft-xborders(procxind(proc))+3
    if(PrtlWall_local.gt.1) then 
		do n=1,used_prtl_arr_size
			if((xp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) up(n)=-up(n)
		end do
		
		do n=1,used_test_prtl_arr_size
			if((xtp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) utp(n)=-utp(n)
		end do		
    end if 
	
	PrtlWall_local=PrtlWallRight-xborders(procxind(proc))+3
    if(PrtlWall_local.lt.mx) then
		do n=1,used_prtl_arr_size
			if((xp(n).gt.PrtlWall_local).and.(qp(n).ne.0)) up(n)=-up(n)
		end do
		
		do n=1,used_test_prtl_arr_size
			if((xtp(n).gt.PrtlWall_local).and.(qp(n).ne.0)) utp(n)=-utp(n)
		end do		
	end if  


end subroutine RefPrtl



subroutine FldCondBoundary
	integer :: xcond_local
	xcond_local=xcond_left-xborders(procxind(proc))+3
	if(xcond_local.ge.1) then
	     Ey(1:min(mx,xcond_local),:,:)=0.0_psn
	     Ez(1:min(mx,xcond_local),:,:)=0.0_psn
	end if
	
	xcond_local=xcond_right-xborders(procxind(proc))+3
	if(xcond_local.lt.mx) then
	     Ey(max(xcond_local,1):mx,:,:)=0.0_psn
	     Ez(max(xcond_local,1):mx,:,:)=0.0_psn
	end if
	
end subroutine FldCondBoundary



	 
subroutine InitOverride !override the default intial conditons
    call InitEMfield	
end subroutine InitOverride  

!The following subroutine is called after moving particles and depositing the current on grid 
! Any additional step which are setup specific should be implelemnted here 
subroutine PostMovDep

end subroutine PostMovDep

subroutine PreAddCurrent
end subroutine PreAddCurrent

subroutine Finalsubroutines !this subroutine is called at the end of every time step	
	call FldCondBoundary
	call RefPrtl 
end subroutine Finalsubroutines	 


end module setup 