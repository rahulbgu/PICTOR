module thrust_module
	interface thrustscan
	subroutine scan_int( input,N) bind(C,name="scan_int_wrapper")
	use iso_c_binding
	integer(c_int),device:: input(*)
	integer(c_int),value:: N
	end subroutine
	end interface
	
end module thrust_module