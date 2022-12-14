module mem_fld
	use parameters
	use vars
contains 
	
	subroutine InitAuxFld
	    if(allocated(Jx)) deallocate(Jx,Jy,Jz,F0)
	    allocate(Jx(mx,my,mz),Jy(mx,my,mz),Jz(mx,my,mz),F0(mx,my,mz))
        Jx = 0; Jy = 0; Jz=0;
	    	
        ! Filtered electric field is used to move test particles
		if(nMoverEMfilter.gt.0) then 
			 if(allocated(FilteredEx)) deallocate(FilteredEx,FilteredEy,FilteredEz)
			 allocate(FilteredEx(mx,my,mz),FilteredEy(mx,my,mz),FilteredEz(mx,my,mz))
		     FilteredEx=0;FilteredEy=0;FilteredEz=0;
		end if 
		
        !arrays used in mover 
        if(allocated(VecEx)) deallocate(VecEx,VecEy,VecEz,VecBx,VecBy,VecBz,VecJ)
#ifdef twoD	
        allocate(VecEx(4,mx*my),VecEy(4,mx*my),VecEz(4,mx*my))
        allocate(VecBx(4,mx*my),VecBy(4,mx*my),VecBz(4,mx*my))
        allocate(VecJ(8,mx*my,Nthreads))
#else
        allocate(VecEx(8,mx*my*mz),VecEy(8,mx*my*mz),VecEz(8,mx*my*mz))	
        allocate(VecBx(8,mx*my*mz),VecBy(8,mx*my*mz),VecBz(8,mx*my*mz))	
        allocate(VecJ(12,mx*my*mz,Nthreads))
#endif
	    VecEx=0; VecEy=0; VecEz=0; VecBx=0; VecBy=0; VecBz=0;
			 
	end subroutine InitAuxFld
		
	
	subroutine ResizeFld(mx,my,mz,mx_new,my_new,mz_new,xshift,yshift,zshift)
		 integer :: mx, my, mz
		 integer :: mx_new, my_new, mz_new
		 integer :: xshift, yshift, zshift
		 integer, dimension(6) :: ind_new, ind
		 
		 ind_new(1) = max(1,1+xshift)
		 ind_new(2) = min(mx_new,mx+xshift)
		 ind_new(3) = max(1,1+yshift) 
		 ind_new(4) = min(my_new,my+yshift) 
		 ind_new(5) = max(1,1+zshift)
		 ind_new(6) = min(mz_new,mz+zshift) 
#ifdef twoD
         ind_new(5) =1
         ind_new(6) =1 
#endif			 
		 
		 ind(1) = max(1,1-xshift)
		 ind(2) = min(mx,mx_new-xshift)
		 ind(3) = max(1,1-yshift) 
		 ind(4) = min(my,my_new-yshift) 
		 ind(5) = max(1,1-zshift)
		 ind(6) = min(mz,mz_new-zshift) 
#ifdef twoD
         ind(5) =1
         ind(6) =1 
#endif		 
		 call ResizeFldArr(Ex,mx_new,my_new,mz_new,ind_new,ind)
		 call ResizeFldArr(Ey,mx_new,my_new,mz_new,ind_new,ind)
		 call ResizeFldArr(Ez,mx_new,my_new,mz_new,ind_new,ind)
		 call ResizeFldArr(Bx,mx_new,my_new,mz_new,ind_new,ind)
		 call ResizeFldArr(By,mx_new,my_new,mz_new,ind_new,ind)
		 call ResizeFldArr(Bz,mx_new,my_new,mz_new,ind_new,ind)
		 
	end subroutine ResizeFld
	
	subroutine ResizeFldArr(Fld,mx_new,my_new,mz_new,ind_new,ind)
		 integer     :: mx_new, my_new, mz_new
		 integer, dimension(6)     :: ind_new, ind
		 real(psn), dimension(:,:,:), allocatable :: Fld
		 real(psn), dimension(:,:,:), allocatable :: FldTemp
	     
		 allocate(FldTemp(mx_new,my_new,mz_new))

		 FldTemp=0.0_psn 
		 FldTemp(ind_new(1):ind_new(2),ind_new(3):ind_new(4),ind_new(5):ind_new(6)) = Fld(ind(1):ind(2),ind(3):ind(4),ind(5):ind(6))
		   
		 deallocate(Fld)
	     call move_alloc(FldTemp,Fld)
	end subroutine ResizeFldArr
	
end module mem_fld