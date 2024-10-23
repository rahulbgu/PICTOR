module prtl_tag
	use parameters
	use vars
	use communication

contains

    ! Prtl tags are a pair of unique integers (tag ID and Proc. ID )
    ! Proc. ID is rank +1. Incremently  by nproc everytime tag ID is exchausted 
		
	subroutine GetTag(flvID, tag_id, proc_id)
		integer :: flvID
		integer :: tag_id, proc_id
		real(dbpsn)    :: r1
		
		tag_id  = 0
		proc_id = 0
		

		if(FlvrSaveRatio(flvID).gt.0 ) then
			call random_number(r1)
			
			if ( CurrentTagID(FlvID) .ge. 2147483645) then 
				CurrentTagID(FlvID) = 0
				CurrentTagProcID(FlvID) = CurrentTagProcID(FlvID) + nproc
			end if  

			proc_id = CurrentTagProcID(FlvID)

			if(r1*FlvrSaveRatio(flvID).le.1) then
				tag_id =  CurrentTagID(FlvID) +1 ! data for positive tag particles are always saved
			else 
				tag_id =  - (CurrentTagID(FlvID) +1 )  ! negative tag particles are tracked, but the data is saved only if needed
			end if

			CurrentTagID(FlvID) = CurrentTagID(FlvID) + 1 

		end if
	end subroutine GetTag
	
end module prtl_tag 