module boundary
     use parameters
     use vars
 contains 
	 subroutine ConductingWallLeft(XWallWave,XWallPrtl)
	    integer :: n
	 	integer :: j,k,xcond_local,PrtlWall_local
	 	real(psn), intent(in) :: XWallWave,XWallPrtl
		real(psn) :: r1,r2
	 	real(psn) :: mag,cos_theta,sin_theta,cos_phi,sin_phi,phi,mag_max
	
	 	xcond_local=XWallWave-xborders(procxind(proc))+3
	 	PrtlWall_local=XWallPrtl-xborders(procxind(proc))+3 !particle's wall is generally few cells ahead of the EM reflector
		
	 	if(PrtlWall_local.lt.0) return 	
		
		do n=1,used_prtl_arr_size
			if((xp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) then 
				xp(n)=PrtlWall_local+(PrtlWall_local-xp(n))
				up(n)=-up(n)
			end if
		end do
		do n=1,used_test_prtl_arr_size
			if((xtp(n).lt.PrtlWall_local).and.(qp(n).ne.0)) then 
				xtp(n)=PrtlWall_local+(PrtlWall_local-xtp(n))
				utp(n)=-utp(n)
			end if
		end do
	
	 	if(xcond_local.ge.1) then 
	 	     Ey(1:min(mx,xcond_local),:,:)=0.0_psn
	 	     Ez(1:min(mx,xcond_local),:,:)=0.0_psn
	 	end if
	 end subroutine ConductingWallLeft
	 
	 subroutine DriftFldRight_Attenuate(x1,Ncell,driftEx,driftEy,driftEz)
	 	integer :: x1
	 	integer :: Ncell !Number of cells in the dissipating zone 
	    integer :: i,j,k,imin
		real(psn), intent(IN) :: driftEx,driftEy,driftEz
	 	real(psn) :: AttenFactor
	
	 	imin=max(1,x1)
	 	if(x1.le.mx) then
	 		 do k=1,mz
	 	          do j=1,my
	 	               do i=imin,mx  
	 					    if(i.ge.x1+Ncell) then
	 		                    Ex(i,j,k)=driftEx
	 		                    Ey(i,j,k)=driftEy
	 		                    Ez(i,j,k)=driftEz        
	 					    else 
	 						   	AttenFactor=(real(abs(x1-i),psn)/real(Ncell,psn))**3							
	 	   	                    Ex(i,j,k)=Ex(i,j,k)-(Ex(i,j,k)-driftEx)*AttenFactor
	 	   	                    Ey(i,j,k)=Ey(i,j,k)-(Ey(i,j,k)-driftEy)*AttenFactor
	 	   	                    Ez(i,j,k)=Ez(i,j,k)-(Ez(i,j,k)-driftEz)*AttenFactor
	 					   end if 
	 	               end do
	 	          end do
	 	     end do
	 	end if 	
	 end subroutine DriftFldRight_Attenuate
	 
	 subroutine ReflCurrentPrtlWall(XWallPrtl)
		real(psn) :: XWallPrtl 
	 	integer :: PrtlWall_local
	
	 	PrtlWall_local=XWallPrtl-xborders(procxind(proc))+3
 		if(PrtlWall_local.ge.3.and.PrtlWall_local.le.mx-2) then
 			!To  ensure that the current is deposited on right place for reflected particles
 			Jx(PrtlWall_local,:,:)=Jx(PrtlWall_local,:,:)-Jx(PrtlWall_local-1,:,:)
 		 	Jy(PrtlWall_local,:,:)=Jy(PrtlWall_local,:,:)+Jy(PrtlWall_local-1,:,:)
 		 	Jz(PrtlWall_local,:,:)=Jz(PrtlWall_local,:,:)+Jz(PrtlWall_local-1,:,:)
 			Jx(PrtlWall_local-1,:,:)=0.0_psn
 			Jy(PrtlWall_local-1,:,:)=0.0_psn
 			Jz(PrtlWall_local-1,:,:)=0.0_psn
 		end if
	 end subroutine ReflCurrentPrtlWall

end module boundary