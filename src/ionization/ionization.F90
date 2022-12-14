module ionization
	use parameters
	use vars
	use interpolation
	implicit none
	
	type IonizationProperty
		character (len=4) :: Type
		integer :: FlvrFrom, FlvrTo
		procedure(func1D) :: radn_ionzn_prob
	end type IonizationProperty
	type(IonizationProperty), dimension(100) :: ionization_list 
   	integer :: used_ionization_list_size = 0
	
contains 
	
	subroutine Ionize
		integer :: n
		
		do n =1, used_ionization_list_size
			if(ionization_list(n)%Type.eq.'radn') call RadiationIonization(ionization_list(n))
		end do 
	
	end subroutine Ionize
	
	subroutine RadiationIonization(prpt)
		type(IonizationProperty) :: prpt
		integer :: n, tag, FlvrTo, FlvrFrom
		real(psn) :: r1, r2 
		real(psn) :: e_x,e_y,e_z,b_x,b_y,b_z, e2b2, wt
		
		FlvrTo = prpt%FlvrTo
		FlvrFrom = prpt%FlvrFrom
		
		do n=1,used_prtl_arr_size
			if(flvp(n).eq.0) cycle
			if(flvp(n).eq.FlvrFrom .and. flvp(n).ne.0) then
				if(qp(n).gt.0 .and. qp(n).lt.FlvrCharge(flvp(n))) then 
					call random_number(r1)
					call Interp_VecEM_GridPoints(x,y,z,e_x,e_y,e_z,b_x,b_y,b_z)
					e2b2 = e_x**2 + e_y**2 + b_y**2 + b_z**2
					
					if(r1.lt.radn_ionzn_prob(e2b2)) then 
						wt = var1p(n) !var1 is used to save wt when ionization physics is included
						qp(n) = qp(n) + wt
						
						tag = GetTag(FlvrTo) 
						
						call InsertParticleAt(used_prtl_arr_size+1,xp(n),yp(n),zp(n),up(n),vp(n),wp(n),-wt,tag, FlvrTo, wt) 
		   			    used_prtl_arr_size=used_prtl_arr_size+1
		   			    np=np+1
					end if
				
				end if 
			end do
		end do 	
		
	end subroutine RadiationIonization
	
	
end module ionization