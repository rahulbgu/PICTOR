module bc_adjust_speed
	use parameters
	use vars
contains 





!--------------------------------------------------------------------------------------------
!
!        Adjust Position of the injector automatically in the shock simulation
!
!--------------------------------------------------------------------------------------------	
	subroutine  AdjustInjectorPositionFromShock(Start, Freq, MinDist, MaxDist, MinInjSpeed, BpeakXmin)
		integer :: Start, Freq, MinDist, MaxDist, BpeakXmin
		real(psn) :: MinInjSpeed, x1, x2, dist, x3
		integer :: Xshock 
	
		if(modulo(t,Freq).ne.0) return
		if(t.lt.Start) return
	
		call FindPosX_BmagPeak(Xshock, BpeakXmin)
	
		!print*,'Xshock is',Xshock
	
		dist = bc_face(2)%pos_prtl - Xshock
		x3 = (MaxDist - MinDist)/3.0
		x1 = MinDist + x3
		x2 = MaxDist - x3
	
		if(dist .lt. x1) bc_face(2)%speed = bc_face(2)%speed + (x1 - dist)/x3 
		if(dist .gt. x2) bc_face(2)%speed = bc_face(2)%speed - (dist - x2)/x3
	
		bc_face(2)%speed = min(bc_face(2)%speed, 1.0_psn)
		bc_face(2)%speed = max(bc_face(2)%speed, MinInjSpeed) 
	end subroutine  AdjustInjectorPositionFromShock


!--------------------------------------------------------------------------------------------
!
!        Adjust speed of the left and right boundaries according to the prescribed distance from the laser pulse
!        Position of the laser is defined by maximum of the magentic field magnitude
!--------------------------------------------------------------------------------------------	
	subroutine  AdjustBoundarySpeedWakeField(Axis, Freq, DistLeft, DistRight, MinSpeedLeft, MinSpeedRight, minEmag2)
		integer :: Freq
		real(psn) :: DistLeft, DistRight, MinSpeedLeft, MinSpeedRight
		real(psn) :: minEmag2
		real(psn) :: x1, x2, dist, x3, MinDist, MaxDist
		integer :: xpeak, axis_ind
		character (len=*) :: Axis

		if(modulo(t,Freq).ne.0) return
		
		axis_ind = 0 
		if(Axis.eq.'x') then 
			call FirstOvershootFromRightX(Ex,Ey,Ez,minEmag2,xpeak)
		end if
	    if(Axis.eq.'z') then 
			axis_ind = 2
			call FirstOvershootFromRightZ(Ex,Ey,Ez,minEmag2,xpeak)
		end if

	
		
		!right boundary
		MinDist = DistRight
		MaxDist = DistRight + c_ompe

		dist = bc_face(2*axis_ind+2)%pos_prtl - xpeak
		x3 = (MaxDist - MinDist)/3.0
		x1 = MinDist + x3
		x2 = MaxDist - x3

		if(dist .lt. x1) bc_face(2*axis_ind+2)%speed = min(bc_face(2*axis_ind+2)%speed + (x1 - dist)/x3, 1.05_psn) 
		if(dist .gt. x2) bc_face(2*axis_ind+2)%speed = max(bc_face(2*axis_ind+2)%speed - (dist - x2)/x3, MinSpeedRight)

		if(proc.eq.0) print*,'laser pulse front is:',xpeak,'bc speed 2 ',bc_face(2*axis_ind+2)%speed,'dist',dist,'MinDist',MinDist

		!left bundary
		MinDist = DistLeft
		MaxDist = DistLeft + c_ompe
		dist = xpeak - bc_face(2*axis_ind+1)%pos_prtl 
		x3 = (MaxDist - MinDist)/3.0
		x1 = MinDist + x3
		x2 = MaxDist - x3
		
		if(dist .lt. x1) bc_face(2*axis_ind+1)%speed = max(bc_face(2*axis_ind+1)%speed - (x1 - dist)/x3, MinSpeedLeft) 
		if(dist .gt. x2) bc_face(2*axis_ind+1)%speed = min(bc_face(2*axis_ind+1)%speed + (dist - x2)/x3, 1.05_psn)
		
		if(proc.eq.0) print*,'laser pulse front is:',xpeak,'bc speed 1 ',bc_face(2*axis_ind+1)%speed,'dist',dist,'MinDist',MinDist
	end subroutine  AdjustBoundarySpeedWakeField


!--------------------------------------------------------------------------------------------
!
!        Find x-position of the global magnetic peak (averaged along y,z direction)
!
!--------------------------------------------------------------------------------------------

	!the following subroutines tries to find the current location of the shock by following total magnetic energy vs. x
	subroutine FindPosX_BmagPeak(Xpeak, BpeakXmin)
		implicit none
		integer :: Xpeak, BpeakXmin
		integer :: i
		real(dbpsn) :: mag_peak
		real(dbpsn), dimension(:), allocatable :: mag1D

		Xpeak = 0

#ifdef gpu
		Bx = Bx_gpu; By=By_gpu; Bz=Bz_gpu;
#endif
        allocate( mag1D( xborders(nSubDomainsX)-xborders(0)) )

        call SumFldX1D(Bx,By,Bz,BpeakXmin,mag1D)
		
		if(kproc.eq.0) then

			 !find the peak position
			 Xpeak = 0
			 mag_peak = 0
			 do i=1,xborders(nSubDomainsX)-xborders(0)
				 if(mag1D(i).gt.mag_peak) then 
					 mag_peak = mag1D(i)
					 xpeak = (i-1) + xborders(0)
				 end if
			 end do 

		 end if
		 
		 deallocate(mag1D)
	 
		 !Step 3 : communicate the xpeak along indepLBaxis
		 call MPI_Bcast( Xpeak,1, MPI_INTEGER, 0, comm_indepLBaxis)
	 
	end subroutine FindPosX_BmagPeak

!--------------------------------------------------------------------------------------------
!
!   global 1D array - sum all elements along y,z direction
!   Note :: 1) indepLBaxis must be along the x-axis to use this routine
!           2) the following routine produces on global 1D arrays on krpoc = 0 only
!              further processing should be done only on kproc=0 and then comm. to other proc.      
!--------------------------------------------------------------------------------------------	

	subroutine SumFldX1D(Fldx,Fldy,Fldz,sum_xmin,fld1D)
		real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
	    real(dbpsn), dimension(mx) :: fld1D_this_proc
		real(dbpsn), dimension(:)              :: fld1D
		real(dbpsn), dimension(:), allocatable :: fld1D_recv
	    integer :: i,j,k, k1
		integer :: sum_xmin, sum_xmin_local, mx_recv, off 
		type(MPI_Status) :: mpi_stat
		
		
		!first compute the local peak and its location
		fld1D_this_proc=0
	    do i=3,mx-3
#ifndef twoD
	          do k=3,mz-3
#else
	          do k=1,1
#endif
	                do j=3,my-3
	                    fld1D_this_proc(i)=fld1D_this_proc(i)+(Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2)
	                end do
	          end do
		 end do

		 ! make sure that any x less than sum_xmin are set to zero
		 sum_xmin_local = sum_xmin -xborders(procxind) +3
		 do i=3,mx-3
			 if(i.lt.sum_xmin_local) fld1D_this_proc(i) = 0
 		 end do

		 allocate( fld1D_recv( xborders(nSubDomainsX)-xborders(0) + 5) )
	 
		 fld1D = 0.0_dbpsn

	     ! Step 1:  comm. along indepLBaxis :: send local sum to kproc=0
		 if(kproc.eq.0) then

			 !add mag1D of kproc=0 into the larger 1D array
			 do i=3,mx-3
				 fld1D(i-2) = fld1D(i-2) +  fld1D_this_proc(i)
			 end do

			 do k1 = 1,nSubDomainsX-1
				mx_recv = xborders(k1+1) - xborders(k1) +5
			 	call MPI_RECV(fld1D_recv(1:mx_recv),mx_recv,MPI_DOUBLE_PRECISION,ProcGrid(iproc,jproc)%procs(k1),k1,MPI_COMM_WORLD,mpi_stat)
				!add mag1D of kproc=k1 into the larger 1D array
				off = xborders(k1)-xborders(0) 
			
				do i=3,mx_recv-3
					fld1D(i-2 + off) = fld1D(i-2 + off) + fld1D_recv(i)
				end do
		
			 end do
		 else
			 call MPI_SEND(fld1D_this_proc,mx,MPI_DOUBLE_PRECISION,ProcGrid(iproc,jproc)%procs(0),kproc,MPI_COMM_WORLD)
		 end if
		 
		 !Step 2 : comm among kproc=0 to get global mag1D on each proc with kproc=0
		 if(kproc.eq.0) then
			 !print*,'before mag1D',mag1D
             call MPI_AllREDUCE(MPI_IN_PLACE,fld1D,xborders(nSubDomainsX)-xborders(0),MPI_DOUBLE_PRECISION,mpi_sum, comm_kproc0)
			 
		 end if
		 
		 deallocate(fld1D_recv)

		
	end subroutine SumFldX1D
	
!--------------------------------------------------------------------------------------------
!
!   find location of the first place from right where the magnitude^2 of vector Fldx,Fldy,Fldz exceed "mag2_min" thresold
!     
!--------------------------------------------------------------------------------------------		

	subroutine FirstOvershootFromRightX(Fldx,Fldy,Fldz,mag2_min,pos)
		real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
		real(psn) :: mag2_min
		integer   :: pos
		integer   :: pos_this
	    integer   :: i,j,k
	
		
		pos = -huge(1) 
		pos_this = -huge(1) 
		!local scan
	    do i=mx-3,3,-1
#ifndef twoD
	          do k=3,mz-3
#else
	          do k=1,1
#endif
	                do j=3,my-3
						if((Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2).gt.mag2_min) pos_this = i + xborders(procxind) -3
	                end do
	          end do
			  if(pos_this.gt.-huge(1)) exit
		 end do

         !global overshoot position "pos" : global maximium of the local overshoot positions
         call MPI_AllREDUCE(pos_this,pos,1,MPI_INTEGER,mpi_max, MPI_COMM_WORLD)

	end subroutine FirstOvershootFromRightX
	
	subroutine FirstOvershootFromRightZ(Fldx,Fldy,Fldz,mag2_min,pos)
		real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
		real(psn) :: mag2_min
		integer   :: pos
		integer   :: pos_this
	    integer   :: i,j,k
	
		
		pos = -huge(1) 
		pos_this = -huge(1) 
		!local scan
	    do k=mz-3,3,-1
	          do j=3,my-3
	                do i=3,mx-3
						if((Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2).gt.mag2_min) pos_this = k + zborders(proczind) -3
	                end do
	          end do
			  if(pos_this.gt.-huge(1)) exit
		 end do

         !global overshoot position "pos" : global maximium of the local overshoot positions
         call MPI_AllREDUCE(pos_this,pos,1,MPI_INTEGER,mpi_max, MPI_COMM_WORLD)

	end subroutine FirstOvershootFromRightZ
		
	
end module bc_adjust_speed