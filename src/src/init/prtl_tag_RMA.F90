module prtl_tag
	use parameters
	use vars
	use communication
	use, intrinsic :: ISO_C_BINDING
	implicit none 
	TYPE(C_PTR) :: cptr_TagBlock
	type(MPI_win) :: TagBlock_Win
	integer :: TagBlock_increment
	integer, dimension(:), allocatable :: TagBlock_temp
contains

    ! Prtl tags are unique integers on every MPI rank
    ! Each MPI rank requests distinct integer blocks from proc 0 using one-sided communication


	
	subroutine InitPrtlTag
		integer (KIND=MPI_ADDRESS_KIND) :: lb, size_of_int
		integer (KIND=MPI_ADDRESS_KIND) :: arr_size, target_disp 
		integer :: disp_unit
		
		call MPI_Type_get_extent(MPI_INTEGER, lb, size_of_int)
		arr_size = Nflvr*size_of_int
		disp_unit = size_of_int

		call MPI_Win_allocate(arr_size, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, cptr_TagBlock, TagBlock_win)
 		call C_F_POINTER(cptr_TagBlock, TagBlock, (/Nflvr/) )

 		call MPI_Win_fence(0,  TagBlock_Win)

        TagBlock = 0 ! initialise 
		
	end subroutine InitPrtlTag
	
	subroutine ReallocatePrtlTagArr
		integer :: old_size
		old_size = size(TagBlock)

		allocate(TagBlock_temp(old_size))
		TagBlock_temp(1:old_size) = TagBlock(1:old_size)
		
		call MPI_Win_Free(TagBlock_win)
		call InitPrtlTag
	    
		TagBlock(1:old_size) = TagBlock_temp(1:old_size)
		deallocate(TagBlock_temp)
	end subroutine ReallocatePrtlTagArr
	
	integer function GetTagBlock(flvID)
		integer (KIND=MPI_ADDRESS_KIND) :: target_disp
		integer :: flvID
		integer :: orig, res
		target_disp = flvID-1
		orig=1
		call MPI_Get_Accumulate(orig,1,MPI_INTEGER,res,1,MPI_INTEGER,0,target_disp,1,MPI_INTEGER,MPI_SUM,TagBlock_Win)
		GetTagBlock = res
	end function GetTagBlock
	
	integer function GetTag(flvID)
		integer :: flvID
		real(dbpsn)    :: r1
		GetTag = 0
		call random_number(r1)
		if(FlvrSaveRatio(flvID).gt.0 .and. r1*FlvrSaveRatio(flvID).le.1) then 
			if(mod(CurrentTagID(FlvID),NtagProcLen).eq.0) CurrentTagID(FlvID) = NtagProcLen*GetTagBlock(flvID) ! old block exhausted, request a new block
			GetTag =  CurrentTagID(FlvID) + 1 
			CurrentTagID(FlvID) = CurrentTagID(FlvID) + 1 
		end if
	end function GetTag
	
end module prtl_tag 