module setup
     use parameters
     use vars
     use help_setup
	 use bc
     implicit none 
	 !---------------------------------------------------------------------------------------
	 !               ||   Physical Parameters  ||
	 !1) All parameters are dimensionless
	 !2) All speeds are in units of c, unless specified otherwise
	 !---------------------------------------------------------------------------------------
	 
	 real(psn) :: Ti, Te !Ion and electron temperature
	 
contains 
	
subroutine InitSetup
	 
	 Ti = 0.01 !in units of rest mass energy, dimensionless ion temp. Ti = kT_i/m_ic^2 
     Te = Ti*(mi/me)! dimensionless elc. temp. Te = kT_e/m_ec^2 
	 
		 
	 call InitPrtl(Flvr1=1, Flvr2=2, Density=Den, Temperature1=TempIon, Temperature2=TempElc, DriftVelocity1=DriftLeft, DriftVelocity2=DriftRight)
     call InitPrtl(Flvr1=1, Flvr2=2, Density=Den, Temperature1=TempIon, Temperature2=TempElc, DriftVelocity1=DriftRight, DriftVelocity2=DriftLeft)
				
end subroutine InitSetup

subroutine BoundaryConditions
end subroutine BoundaryConditions

subroutine DriftLeft(x,y,z,vx,vy,vz)
	real(dbpsn) :: x,y,z
	real(psn)   :: vx,vy,vz
	vx = -g0
	vy = 0
	vz = 0
end subroutine DriftLeft 
subroutine DriftRight(x,y,z,vx,vy,vz)
	real(dbpsn) :: x,y,z
	real(psn)   :: vx,vy,vz
	vx = g0
	vy = 0
	vz = 0
end subroutine DriftRight 
	
real(psn) function Den(x,y,z)
     real(dbpsn) :: x,y,z
	 Den=0.5_psn
end function Den



real(psn) function TempIon(x,y,z)
     real(dbpsn) :: x,y,z
	 TempIon=Ti
end function TempIon

real(psn) function TempElc(x,y,z)
     real(dbpsn) :: x,y,z
	 TempElc=Te
end function TempElc






!================================================================================
!The following subroutines are meant for further customizations and must be part of every setup module
!Any additional actions which are setup-specific should be implelemnted in the following 
!================================================================================
subroutine SaveSetupParameters ! the parameters defined in this setup can be saved in the "param" file 
	call SaveParam('Te', Te)
	call SaveParam('Ti', Ti)
end subroutine SaveSetupParameters

subroutine InitOverride !if needed, override some default initial conditons 
end subroutine InitOverride  

!This subroutine is called right after moving particles and depositing current on the grid
subroutine PostMovDep
end subroutine PostMovDep

!This subroutine is called right before depositing current on the grid
subroutine PreAddCurrent
end subroutine PreAddCurrent

!This subroutine is called at the end of every time step
subroutine FinalSubroutines 
end subroutine FinalSubroutines


end module setup
