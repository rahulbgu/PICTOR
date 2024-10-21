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
		integer :: Start, Freq
		real(psn) :: MinDist, MaxDist, BpeakXmin
		real(psn) :: MinInjSpeed, x1, x2, dist, x3
		integer :: Xshock 
	
		if(modulo(t,Freq).ne.0) return
		if(t.lt.Start) return
	
		call FindPosX_BmagPeak(Xshock, int(BpeakXmin))
	
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

        call SumFld1D(Bx,By,Bz,BpeakXmin,mag1D,0)
		
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

	subroutine SumFld1D(Fldx,Fldy,Fldz,sum_xmin,fld1D,axis)
		integer, intent(in)            :: axis ! sum transverse to the axis
		real(psn), dimension(mx,my,mz) :: Fldx, Fldy, Fldz
	    real(dbpsn), dimension(:), allocatable :: fld1D_this_proc
		real(dbpsn), dimension(:)              :: fld1D
		real(dbpsn), dimension(:), allocatable :: fld1D_recv
	    integer :: i,j,k, k1 , b1,b2, nx_axis, mx_axis
		integer :: sum_xmin, sum_xmin_local, mx_recv, off 
		type(MPI_Status) :: mpi_stat

		!set local variables based on the axis
		select case(axis)
		case(0)
			b1 = xborders(procxind)
			b2 = xborders(procxind+1)
			nx_axis = xborders(nSubDomainsX) -xborders(0)
		case(1)
			b1 = yborders(procyind)
			b2 = yborders(procyind+1)
			nx_axis = yborders(nSubDomainsY) -yborders(0)
		case(2)
			b1 = zborders(proczind)
			b2 = zborders(proczind+1)
			nx_axis = zborders(nSubDomainsZ) -zborders(0)
		end select
		mx_axis = b2 - b1 +5
		
		
		!first compute the local peak and its location
		allocate(fld1D_this_proc(b2-b1+5))
		fld1D_this_proc=0
#ifndef twoD
	    do k=3,mz-3
#else
	    do k=1,1
#endif
	        do j=3,my-3
				do i=3,mx-3
					select case(axis)
					case(0)
	                	fld1D_this_proc(i)=fld1D_this_proc(i)+(Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2)
					case(1)
						fld1D_this_proc(j)=fld1D_this_proc(j)+(Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2)
					case(2)
	                	fld1D_this_proc(k)=fld1D_this_proc(k)+(Fldx(i,j,k)**2+Fldy(i,j,k)**2+Fldz(i,j,k)**2)
					end select
	            end do
	        end do
		end do

		 ! make sure that any x less than sum_xmin are set to zero
		 sum_xmin_local = sum_xmin - b1 +3
		 do i=3,mx_axis-3
			 if(i.lt.sum_xmin_local) fld1D_this_proc(i) = 0
 		 end do

		 allocate( fld1D_recv( nx_axis + 5) )
	 
		 fld1D = 0.0_dbpsn

	     ! Step 1:  comm. along indepLBaxis :: send local sum to kproc=0
		 if(kproc.eq.0) then

			 !add mag1D of kproc=0 into the larger 1D array
			 do i=3,mx_axis-3
				 fld1D(i-2) = fld1D(i-2) +  fld1D_this_proc(i)
			 end do

			 do k1 = 1,nSubDomainsX-1
				select case(axis)
				case(0)
					mx_recv = xborders(k1+1) - xborders(k1) +5
					off = xborders(k1)-xborders(0)
				case(1)
					mx_recv = yborders(procyind+1) - yborders(procyind) +5
					off = yborders(procyind)-yborders(0)
				case(2)
					mx_recv = zborders(proczind+1) - zborders(proczind) +5
					off = zborders(proczind)-zborders(0)
				end select

			 	call MPI_RECV(fld1D_recv(1:mx_recv),mx_recv,MPI_DOUBLE_PRECISION,ProcGrid(iproc,jproc)%procs(k1),k1,MPI_COMM_WORLD,mpi_stat)
				
				!add mag1D of kproc=k1 into the larger 1D array			
				do i=3,mx_recv-3
					fld1D(i-2 + off) = fld1D(i-2 + off) + fld1D_recv(i)
				end do
		
			 end do
		 else
			 call MPI_SEND(fld1D_this_proc,mx_axis,MPI_DOUBLE_PRECISION,ProcGrid(iproc,jproc)%procs(0),kproc,MPI_COMM_WORLD)
		 end if
		 
		 !Step 2 : comm among kproc=0 to get global mag1D on each proc with kproc=0
		 if(kproc.eq.0) then
			 !print*,'before mag1D',mag1D
             call MPI_AllREDUCE(MPI_IN_PLACE,fld1D,nx_axis,MPI_DOUBLE_PRECISION,mpi_sum, comm_kproc0)
		 end if
		 
		 deallocate(fld1D_recv)
		 deallocate(fld1D_this_proc)

	end subroutine SumFld1D
	
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