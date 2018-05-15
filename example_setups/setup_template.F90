!The template file
!A new problem should be produced by constructing the essential subroutines listed here
module setup

use parameters
use vars
implicit none 

contains

!This subroutine is called just once at the very begining of the simulation. Initialize all the setup related variables here.       
subroutine InitUser

end subroutine InitUser

!This subroutine is called just one right before entering the main interation loop. The default values of the vriables can be changed at this point.       
subroutine InitOverride

end subroutine InitOverride  

! Some setup sepcific actions can be defiend in the folliwng three subroutine
!This subroutine is called once after the end of every time step. Usage example: add some new particles at a boundary. 
subroutine Finalsubroutines 
     
end subroutine Finalsubroutines


!This subroutine is called right after moving the particles. Usage example: scatter particles
subroutine PostMovDep

end subroutine PostMovDep

!This subroutine is called right before including the current in the Maxwell's equations. Usage example: add some external (antena) current
subroutine PreAddCurrent

end subroutine PreAddCurrent


end module setup