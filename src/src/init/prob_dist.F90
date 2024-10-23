module prob_dist
	use parameters
	use vars
	implicit none	
	real(psn), dimension(:), allocatable :: MB_Table, MB_PDF_Table
    logical :: init_maxbolt=.false.
		
contains  
	
	!The following subroutine returns three-velocity (\gamma*v/c) corresponding to a drifting Maxwell-Boltzman distribution 
	!Drift= drift speed/c (if Drift<1) OR Lorentz factor (if Drift>1)
	!Temp= (2kT/mc^2), where T is in the physical units (Temp is dimensionless)
	subroutine GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma) 
		 integer                :: TableSize=10000    
	     real(psn), intent(in)  :: Temp
	     real(psn), intent(out) :: ugamma,vgamma,wgamma !\gamma * v/c
	     real(psn) :: r1,gamma,beta
	     integer :: index
		 
	     if(init_maxbolt.eqv..false.) call Init_MaxBolt_Table(TableSize)
		 call random_number(r1)     
	     call BinarySearch(TableSize,MB_PDF_Table,r1,index)  
		 call InvPDF(TableSize,MB_Table,MB_PDF_Table,r1,index,beta)
		
		 beta=sqrt(beta*Temp)
	     if(beta.ge.1.0_psn) beta=0.0 !in some very rare cases beta can be equal to 1
	       									   
		 gamma=1.0_psn/sqrt((1.0_psn+beta)*(1.0_psn-beta))
		 call iso3vel(gamma,beta,ugamma,vgamma,wgamma)
	end subroutine GetVelGamma_MaxBolt

	
	!Numerical values of the dimensionless quantity (v/c)^2 / (2kT/mc^2) are stored in MB_Table 
	subroutine Init_MaxBolt_Table(N)
		 integer :: N
	     real(psn) :: dX, PDF_sum
	     integer :: i
		 
		 if(.not.allocated(MB_Table)) allocate(MB_Table(N),MB_PDF_Table(N))
	     dX=(15.0_psn-0.0_psn)/(N-1) !Min X = 0.0, Max X = 5.0
 
	     do i=1,N
	          MB_Table(i)=dX*(i-1)
	     end do  
	    
	     PDF_sum=0
	     MB_PDF_Table(1)=0
	     do i=2,N
	       PDF_sum=PDF_sum+sqrt(MB_Table(i))*exp(-MB_Table(i))
	       MB_PDF_Table(i)=PDF_sum
	     end do		 
	     MB_PDF_Table=MB_PDF_Table/PDF_sum !Normalisation
		 init_maxbolt=.true.
	end subroutine Init_MaxBolt_Table
	
	subroutine iso3vel(gamma,beta,ugamma,vgamma,wgamma)
		real(psn) :: gamma,beta,ugamma,vgamma,wgamma
		real(psn) :: r2, r3
	    real(psn) :: phi,CosTheta,SinTheta
		
	    call random_number(r2)
	    call random_number(r3)
	    CosTheta=2*r2-1 
	    SinTheta=sqrt(1-CosTheta**2)
	    phi=2*pi*r3
		
	    ugamma=gamma*beta*CosTheta
	    vgamma=gamma*beta*SinTheta*cos(phi)
	    wgamma=gamma*beta*SinTheta*sin(phi)
	end subroutine iso3vel
	
	subroutine InvPDF(N,X_Table,PDF_Table,PDF,index,X)
		integer :: N,index
		real(psn) :: PDF,X
		real(psn), dimension(N) :: X_Table, PDF_Table
		
		if(PDF_Table(index+1).eq.PDF_Table(index)) then 
			 X=X_Table(index)
		else    
	         X=X_Table(index) + ( X_Table(index+1)-X_Table(index) )* ( (PDF-PDF_Table(index)) / (PDF_Table(index+1)-PDF_Table(index)) ) 
		end if	
		
	end subroutine InvPDF
	
	subroutine BinarySearch(N,Arr,val,index)
	     integer :: N
	     real(psn),intent(in), dimension(N):: Arr
	     real(psn),intent(in) :: val
	     integer :: index,imin,imax,imid,imid_prev
	     logical :: search
     
	     imin=1
	     imax=N
		 imid=(imin+imax)/2
	     search=.true. 
		 do while(search)
			  imid_prev=imid 
	          if(Arr(imid).ge.val) imax=imid
	          if(Arr(imid).le.val) imin=imid
			  imid=(imin+imax)/2
	          if(imin+1.ge.imax) search=.false.
			  if(imid.eq.imid_prev) search=.false. !to ensure exit from the loop if Arr(imin)=Arr(imax) and imin+1<imax 		  
	     end do 
	     index=imin	 
	end subroutine BinarySearch 
	
	!--------------------------------------------------------------------------------------------
	!                    Initialise a PDF table from user defined function
	!--------------------------------------------------------------------------------------------
	subroutine InitPDFTable(N,X_Table,PDF_Table,PDF_FUNC,vmax)
		integer :: N 
		real(psn), dimension(:), allocatable :: PDF_Table, X_Table
		procedure(func1D) :: PDF_FUNC
		real(psn) :: vmax
		real(psn) :: dX, PDF_sum
		integer   :: i 
		
		if(allocated(X_Table)) deallocate(X_Table,PDF_Table)
		allocate(X_Table(N),PDF_Table(N))
        
		dX=(vmax-0.0_psn)/(N-1) !Min X = 0.0, Max X = 5.0

        do i=1,N
       		X_Table(i)=dX*(i-1)
        end do  
		
	    PDF_sum=0
	    PDF_Table(1)=0
	    do i=2,N
	    	PDF_sum=PDF_sum + PDF_FUNC(X_Table(i))
	        PDF_Table(i)=PDF_sum
	    end do		 
	    PDF_Table=PDF_Table/PDF_sum !Normalisation
		
	end subroutine InitPDFTable
	
	!--------------------------------------------------------------------------------------------
	!                    Use a speed distribution PDF table to get isotropic 3-velocity 
	!--------------------------------------------------------------------------------------------
	subroutine GetIsoVelGammaTable(N,X_Table,PDF_Table,ugamma,vgamma,wgamma)     
		 integer :: N
		 real(psn), dimension(N) :: X_Table, PDF_Table
	     real(psn), intent(out) :: ugamma,vgamma,wgamma !\gamma * v/c
	     real(psn) :: r1,gamma,beta
	     integer :: index
		 
		 call random_number(r1)     
	     call BinarySearch(N,PDF_Table,r1,index)  
		 call InvPDF(N,X_Table,PDF_Table,r1,index,beta)
		
	     if(beta.ge.1.0_psn) beta=0.0 !in some very rare cases beta can be equal to 1
	       									   
		 gamma=1.0_psn/sqrt((1.0_psn+beta)*(1.0_psn-beta))
		 call iso3vel(gamma,beta,ugamma,vgamma,wgamma)
	end subroutine GetIsoVelGammaTable
	
	!--------------------------------------------------------------------------------------------
	!      Get isotropic 3-velocity using an arbitrary speed dist. function
	!--------------------------------------------------------------------------------------------
	subroutine GetVelGammaNonThermal(VelDist,ugamma,vgamma,wgamma,vmax)
		procedure(scalar_local) :: VelDist
		real(psn) :: ugamma,vgamma,wgamma,vmax
		real(psn) :: vx,vy,vz,f,gamma
		real(psn) :: r1,r2,r3,rnd_acpt
		logical :: search
	
		search=.true.
		do while(search)
			call random_number(r1)
			call random_number(r2)
			call random_number(r3)
			call random_number(rnd_acpt)
			vx=(2*r1-1)*vmax
			vy=(2*r2-1)*vmax
			vz=(2*r3-1)*vmax
			f=VelDist(vx,vy,vz)
			if(rnd_acpt.le.f) then
				gamma=1.0_psn/sqrt(1.0_psn-(vx*vx + vy*vy + vz*vz))
				ugamma=gamma*vx
				vgamma=gamma*vy
				wgamma=gamma*vz
				search=.false. 
			end if
	    end do 
	end subroutine GetVelGammaNonThermal
	
	
	!--------------------------------------------------------------------------------------------
	!
	!                    Poisson distribution 
	!
	!--------------------------------------------------------------------------------------------
	
	subroutine GetIntPoissonDist(mean,num)
		real(psn)   :: mean 
		integer     :: num
		real(dbpsn)   :: r1
		integer     :: n
		real(dbpsn) :: sum, sum_logn
		real(dbpsn) :: log_mean
		
		call random_number(r1)  
		sum = exp(-mean)
		sum_logn = 0 
		log_mean = log(mean)
		n=0
		
		do while(sum.lt.r1)
			n=n+1
			sum_logn = sum_logn + log(real(n,dbpsn))
			sum = sum + exp(-Mean + n*log_mean - sum_logn)	
			if(n.gt.10000) exit ! max. 10000 iterations, the mean should be much less than 10000
		end do
		num = n
		
	end subroutine GetIntPoissonDist
	
	
	!------------------------------------------------------------
	!the following subroutine is used to add drift velocity to the thermal velocity
	!------------------------------------------------------------

	subroutine AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
		real(psn),intent(inout):: ugamma,vgamma,wgamma
		real(psn),intent(in)   :: vx,vy,vz
	    real(psn) :: Drift,DriftBeta,DriftGamma,gamma
	    real(psn), dimension(3) :: dirn
	    real(psn) :: dirn_mag,projn

	    Drift=sqrt( vx*vx + vy*vy + vz*vz )
	    if(abs(Drift).lt.1.0_psn) then 
	         DriftBeta=Drift
	         DriftGamma=1.0_psn/sqrt((1.0_psn-Drift)*(1.0_psn+Drift))          
	    else 
	         DriftGamma=abs(Drift)   
	         DriftBeta=sqrt((Drift-1.0_psn)*(Drift+1.0_psn))/Drift          
	    end if


         dirn(1)=vx; dirn(2)=vy; dirn(3)=vz;
         dirn_mag=sqrt(dirn(1)**2+dirn(2)**2+dirn(3)**2)

         if(dirn_mag.eq.0) then 
              dirn=0
         else 
              dirn=dirn/dirn_mag
         end if
         gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)  
         projn=ugamma*dirn(1)+vgamma*dirn(2)+wgamma*dirn(3)
 
         ugamma=ugamma+(DriftGamma-1)*projn*dirn(1) + DriftGamma*DriftBeta*gamma*dirn(1)
         vgamma=vgamma+(DriftGamma-1)*projn*dirn(2) + DriftGamma*DriftBeta*gamma*dirn(2)
         wgamma=wgamma+(DriftGamma-1)*projn*dirn(3) + DriftGamma*DriftBeta*gamma*dirn(3)
	end subroutine AddDriftVel
	
	
	
end module prob_dist