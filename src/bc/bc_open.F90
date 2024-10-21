module bc_open
	use parameters
	use vars
	use mem_prtl
	implicit none
contains 
	
	subroutine RemovePrtl_Left(pos_bc, q,flv,tag,x,y,z,u,v,w)
		integer, dimension(:) :: flv, tag
		real(psn), dimension(:) :: q,x,y,z,u,v,w
		real(psn) :: pos_bc
		integer :: n , count
		
		count = 0		
		do n=1,used_prtl_arr_size
			if(flv(n).eq.0) cycle
			if(x(n).lt.pos_bc) then 
				q(n)=0
				flv(n)=0
				u(n)=0
				v(n)=0
				w(n)=0
				x(n)=3.8_psn
				y(n)=3.8_psn
#ifdef twoD
                z(n)=1.5_psn 
#else
                z(n)=3.8_psn
#endif	
				tag(n) = 0
				count = count +1 
			end if  
		end do
		np=np -count
	end subroutine RemovePrtl_Left
	
	subroutine RemovePrtl_Right(pos_bc, q,flv,tag,x,y,z,u,v,w)
		integer, dimension(:) :: flv, tag
		real(psn), dimension(:) :: q,x,y,z,u,v,w
		real(psn) :: pos_bc
		integer :: n , count
				
		count = 0 
		do n=1,used_prtl_arr_size
			if(flv(n).eq.0) cycle
			if(x(n).gt.pos_bc) then 
				q(n)=0
				flv(n)=0
				u(n)=0
				v(n)=0
				w(n)=0
				x(n)=3.8_psn
				y(n)=3.8_psn
#ifdef twoD
                z(n)=1.5_psn
#else
                z(n)=3.8_psn
#endif				
				tag(n)=0
				count = count +1 
			end if  
		end do
		np = np -count
	end subroutine RemovePrtl_Right
	
	
end module bc_open