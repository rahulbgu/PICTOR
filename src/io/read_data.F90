module read_data
    use parameters
    use vars
    implicit none
    type GridData
        real(psn), dimension(:), allocatable :: arr1D
        real(psn) :: spacing_x !in units of cell-size
    end type GridData


contains

!reads data from FileName and return "GridData"  
function loadGridDataFromFile(FileName, PhysicalUnit, Spacing) result(grid)
    character(len=*), intent(in) :: FileName
    character(len=*), intent(in) :: PhysicalUnit
    real(psn)                    :: Spacing
    type(GridData) :: grid
    integer :: ioStat, funit, n, i
    real(psn) :: temp, spacing_x, skin_depth

    ! Open the file 
    open(newunit=funit, file=FileName, status='old', action='read', iostat=ioStat)
    if (ioStat .ne. 0) then
        if(proc.eq.0) print *, "Failed to open file:", FileName
        stop
    end if


    ! count the number of numbers in the file
    n = 0
    do
        read(funit, *, iostat=ioStat) temp
        if (ioStat .ne. 0) exit
        n = n + 1
    end do
    rewind(funit)  
    allocate(grid%arr1D(n))

    ! read the numbers into arr1D
    do i = 1, n
        read(funit, *) grid%arr1D(i)
    end do

    close(funit)

    !set the spacing

    skin_depth = 5.314*1.0e3/sqrt(ne0) ! in meter
    if(PhysicalUnit.eq.'meter') then 
        spacing_x = (grid_dx*skin_depth / c_ompe) / Spacing
    end if 
    
    grid%spacing_x = spacing_x
end function loadGridDataFromFile


function readFromGridData(x,grid_data) result(val)
    type(GridData), intent(in) :: grid_data
    real(dbpsn), intent(in)    :: x
    real(dbpsn)                :: x1
    integer                    :: ind
    real(psn)                  :: f, val

    x1 = x * grid_data%spacing_x
    ind = x1
    f = x1 - ind 
    ind = ind + 1 

    !print*,'ind is', x , ind, grid_data%spacing_x
    if (ind .ge. 1 .and. ind+1 .le. size(grid_data%arr1D)) then
        !print*,'reading data from here'
        val = grid_data%arr1D(ind) + f * (grid_data%arr1D(ind+1) - grid_data%arr1D(ind))
    else
        val = 0 ! default value for out of bounds
    end if

end function readFromGridData 



end module read_data