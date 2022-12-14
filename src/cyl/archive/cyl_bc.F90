module cyl_bc
    use parameters
    use vars
	use cyl_vars
#ifdef gpu 
    use cyl_bc_gpu
#else
    use cyl_bc_cpu 
#endif  	
	implicit none 
contains 
	
!---------------------------------------------------------------------------------
!  Subroutines to implement BCs from the setup file 
!---------------------------------------------------------------------------------
    subroutine SetBC_Fld(Axis,Position,Type,Side)
		 character (len=*) :: Axis,Type,Side
		 real(psn)  :: Position
		 if((Axis.eq.'r').and.(Type.eq.'cond')) then
			 if(Side.eq.'inner') then 
				 BC_Rmin_Fld=Position
				 BC_Rmin_Fld_Type=Type
			 end if 
			 if(Side.eq.'outer') then 
				 BC_Rmax_Fld=Position
				 BC_Rmax_Fld_Type=Type
			 end if 
		 end if 
	end subroutine SetBC_Fld	 		

    subroutine SetBC_Prtl(Axis,Position,Type,Side)
		 character (len=*) :: Axis,Type,Side
		 real(psn)  :: Position
		 if((Axis.eq.'r').and.(Type.eq.'refl')) then
			 if(Side.eq.'inner') then 
				 BC_Rmin_Prtl=Position
				 BC_Rmin_Prtl_Type=Type
			 end if 
			 if(Side.eq.'outer') then 
				 BC_Rmax_Prtl=Position
				 BC_Rmax_Prtl_Type=Type
			 end if 
		 end if 
	end subroutine SetBC_Prtl	
!---------------------------------------------------------------------------------
!  Essential rotuines to enforece BCs from the main loop 
!---------------------------------------------------------------------------------	
	subroutine EnforceBC_PostMovDep
		if(BC_Rmax_Prtl.gt.0) then
			if(BC_Rmax_Prtl_Type.eq.'refl') call RefOuterBC_Prtl(BC_Rmax_Prtl)
		end if 
		if(.not.inc_axis) then
			if((BC_Rmin_Prtl.ge.0).and.(BC_Rmax_Prtl_Type.eq.'refl')) call RefInnerBC_Prtl(BC_Rmin_Prtl)
		end if
		if(inc_axis) call AxisCurrentBC
	end subroutine EnforceBC_PostMovDep
	
	subroutine EnforceBC_Final
		if(BC_Rmax_Fld.gt.0) then
			if(BC_Rmax_Fld_Type.eq.'cond') call ConductingOuterBC_Fld(BC_Rmax_Fld)
		end if 
		
		if(.not.inc_axis) then
			if((BC_Rmin_Fld.ge.0).and.(BC_Rmax_Fld_Type.eq.'cond')) call ConductingInnerBC_Fld(BC_Rmin_Fld) 
		end if
		!if(inc_axis) call UpdateFldAxis
	end subroutine EnforceBC_Final
	
	
end module cyl_bc