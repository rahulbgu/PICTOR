module usc_CFDTD
	use parameters
	use vars
#ifdef gpu
    use var_gpu
    use usc_CFDTD_gpu
#endif
contains 
	
	
	!------------------------------------------------------------------------------
	! Conformal FDTD method
	!------------------------------------------------------------------------------
	
	 subroutine Init_USC_CFDTD
  		if(allocated(e_lx)) call Deallocate_USC_CFDTD_Vars
  		
		allocate(e_lx(mx,my,mz),e_ly(mx,my,mz),e_lz(mx,my,mz),b_arx(mx,my,mz),b_ary(mx,my,mz),b_arz(mx,my,mz))
		e_lx = 0; e_ly =0; e_lz=0; b_arx=0; b_ary=0; b_arz=0;
		
 		allocate(usc_db1(mx,my,mz),usc_db2(mx,my,mz))
 		usc_db1=0.0_psn; usc_db2=0.0_psn; 
 		
		allocate(usc_fdtd_bx(mx,my,mz),usc_fdtd_by(mx,my,mz),usc_fdtd_bz(mx,my,mz))	
 		usc_fdtd_bx=0; usc_fdtd_by=0; usc_fdtd_bz=0;
	 end subroutine Init_USC_CFDTD
	 
	 subroutine Deallocate_USC_CFDTD_Vars
		 deallocate(e_lx,e_ly,e_lz, b_arx, b_ary, b_arz) 
		 deallocate(usc_db1,usc_db2)
		 deallocate(usc_fdtd_bx,usc_fdtd_by,usc_fdtd_bz)
		 deallocate(usc_wtr_bx,usc_wtr_by,usc_wtr_bz)
		 deallocate(usc_wtc_bx,usc_wtc_by,usc_wtc_bz)
#ifdef gpu
         call DeallocateCurvedBCVars_gpu
#endif		 
	 end subroutine Deallocate_USC_CFDTD_Vars
	 
 	 subroutine UpdateBfieldHalf_CFDTD ! called directly from the main B update routine 
 		!call UpdateBfieldHalf_CFDTD_cyl
#ifdef gpu
 		call UpdateBfieldHalf_USC_gpu
 		return
#endif 
 		call UpdateBfieldHalf_USC
 	 end subroutine UpdateBfieldHalf_CFDTD	
	 
	 
	 subroutine SetExteriorEfld_USC ! called from BC.F90 routine
		 integer :: i,j,k
#ifdef gpu
         call SetExteriorEfld_USC_gpu
		 return 
#endif		 
		 
		 do k=1,mz
		     do j=1,my
		          do i=1,mx
					  if(e_lx(i,j,k).eq.0) Ex(i,j,k)=0.0_psn
					  if(e_ly(i,j,k).eq.0) Ey(i,j,k)=0.0_psn
					  if(e_lz(i,j,k).eq.0) Ez(i,j,k)=0.0_psn
				  end do 
			 end do 
		 end do 
		 
	 end subroutine SetExteriorEfld_USC
	 
	 subroutine SetExteriorBfld_USC ! called from BC.F90 routine
		 integer :: i,j,k
		 
		 do k=1,mz
		     do j=1,my
		          do i=1,mx
					  if(b_arx(i,j,k).eq.0) Bx(i,j,k)=0.0_psn
					  if(b_ary(i,j,k).eq.0) By(i,j,k)=0.0_psn
					  if(b_arz(i,j,k).eq.0) Bz(i,j,k)=0.0_psn
				  end do 
			 end do 
		 end do 
		 
	 end subroutine SetExteriorBfld_USC	 	 
  	
	!------------------------------------------------------------------------------
  	! Modified FDTD for the USC method; the B field update is modified, E is set to 0 inside the conductor
  	!------------------------------------------------------------------------------
 	 subroutine UpdateBfieldHalf_USC
 		 usc_db1=0
 		 usc_db2=0
		 
 		 call Calc_dbx(usc_db1)
 		 call gather_dbx_bc_cells(usc_db1,usc_db2,usc_wtr_bx)
 		 call gather_dbx_bc_cells(usc_db2,usc_db1,usc_wtc_bx)
 		 Bx=Bx+usc_db1

 		 usc_db1=0
 		 usc_db2=0		 
		 
 		 call Calc_dby(usc_db1)
 		 call gather_dby_bc_cells(usc_db1,usc_db2,usc_wtr_by)
 		 call gather_dby_bc_cells(usc_db2,usc_db1,usc_wtc_by)
 		 By=By+usc_db1
		 
 		 usc_db1=0
 		 usc_db2=0
		 
 		 call Calc_dbz(usc_db1)
 		 call gather_dbz_bc_cells(usc_db1,usc_db2,usc_wtr_bz)
 		 call gather_dbz_bc_cells(usc_db2,usc_db1,usc_wtc_bz)
 		 Bz=Bz+usc_db1	 
		 
		 !call SetExteriorBfld_USC
! 		 Bx=0.0_psn
! 		 By=0.0_psn
! 		 Bz=0.0_psn
 	 end subroutine UpdateBfieldHalf_USC
	 
	 
 	 subroutine Calc_dbx(arr)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:,:) :: arr
#ifdef twoD
          do k=1,1
#else			 		 
          do k=1,mz-1
#endif
               do j=1,my-1
                    do i=1,mx	
#ifdef twoD
                       arr(i,j,k)=-fld_halfc*(e_lz(i,j+1,k)*Ez(i,j+1,k)-e_lz(i,j,k)*Ez(i,j,k))
#else
 				       arr(i,j,k)=-fld_halfc*(e_lz(i,j+1,k)*Ez(i,j+1,k)-e_lz(i,j,k)*Ez(i,j,k))+&
 						           fld_halfc*(e_ly(i,j,k+1)*Ey(i,j,k+1)-e_ly(i,j,k)*Ey(i,j,k))
#endif					   					   
 				   end do 
 			  end do
 		 end do  
 	 end subroutine Calc_dbx
	 
 	 subroutine gather_dbx_bc_cells(arr_in,arr_out,wt)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:) :: wt
 		 real(psn), dimension(:,:,:) :: arr_in,arr_out
		 
#ifdef twoD
          do k=1,1
#else			 		 
          do k=2,mz-1
#endif
               do j=2,my-1
                    do i=1,mx
 					    if(usc_fdtd_bx(i,j,k).gt.0) arr_out(i,j,k) = usc_gather_dbx(i,j,k,usc_fdtd_bx(i,j,k),wt,arr_in)
 				   end do 
 			  end do 
 		  end do 	 
		 
 	 end subroutine gather_dbx_bc_cells
	 
 	 subroutine Calc_dby(arr)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:,:) :: arr
#ifdef twoD
          do k=1,1
#else			 		 
          do k=1,mz-1
#endif
               do j=1,my
                    do i=1,mx-1	
#ifdef twoD
                        arr(i,j,k)=fld_halfc*(e_lz(i+1,j,k)*Ez(i+1,j,k)-e_lz(i,j,k)*Ez(i,j,k))
#else
 				        arr(i,j,k)=-fld_halfc*(e_lx(i,j,k+1)*Ex(i,j,k+1)-e_lx(i,j,k)*Ex(i,j,k))+&
 				                   fld_halfc*(e_lz(i+1,j,k)*Ez(i+1,j,k)-e_lz(i,j,k)*Ez(i,j,k))
#endif					   					   
 				   end do 
 			  end do
 		 end do  
		 
 	 end subroutine Calc_dby
	 
 	 subroutine gather_dby_bc_cells(arr_in,arr_out,wt)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:) :: wt
 		 real(psn), dimension(:,:,:) :: arr_in,arr_out
#ifdef twoD
          do k=1,1
#else			 		 
          do k=2,mz-1
#endif
              do j=1,my
                   do i=2,mx-1
 					  if(usc_fdtd_by(i,j,k).gt.0) arr_out(i,j,k) = usc_gather_dby(i,j,k,usc_fdtd_by(i,j,k),wt,arr_in)
 				  end do 
 			 end do
 		  end do   
 	 end subroutine gather_dby_bc_cells
	 
 	 subroutine Calc_dbz(arr)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:,:) :: arr
          do k=1,mz
               do j=1,my-1
                    do i=1,mx-1	
                        arr(i,j,k)=-fld_halfc*(e_ly(i+1,j,k)*Ey(i+1,j,k)-e_ly(i,j,k)*Ey(i,j,k))+&
 					               fld_halfc*(e_lx(i,j+1,k)*Ex(i,j+1,k)-e_lx(i,j,k)*Ex(i,j,k))
 				   end do 
 			   end do 
 		 end do 	 
 	 end subroutine Calc_dbz
	 
 	 subroutine gather_dbz_bc_cells(arr_in,arr_out,wt)
 		 integer :: i,j,k
 		 real(psn), dimension(:,:) :: wt
 		 real(psn), dimension(:,:,:) :: arr_in,arr_out
          do k=1,mz
               do j=2,my-1
                    do i=2,mx-1	
 					   if(usc_fdtd_bz(i,j,k).gt.0) arr_out(i,j,k) = usc_gather_dbz(i,j,k,usc_fdtd_bz(i,j,k),wt,arr_in)
 				   end do 
 			   end do 
 		 end do 	
 	 end subroutine gather_dbz_bc_cells
	 
	 
 	!---------------------------------------------------------------
 	! Initialise weights in the FDTD notation for the USC method (Igor Zagorodnov et al, JCP 2007)
 	!---------------------------------------------------------------
	
 	subroutine InitWeight_USC
 		integer :: i,j,k
 		integer :: count_bx, count_by, count_bz
 		real(psn) :: a=0.5_psn
		

	
 		allocate(usc_norm_bx(mx,my,mz),usc_norm_by(mx,my,mz),usc_norm_bz(mx,my,mz))
		
		
 		count_bx=0; count_by=0; count_bz=0
		
		
 		!first count number of cell facets that require USC fdtd
#ifdef twoD		
 		do k=1,1
#else
        do k=2,mz-1
#endif			
 			do j=2,my-1
 				do i=2,mx-1
 					if(b_arx(i,j,k).gt.0 .and. unstableCell(i,i,j-1,j+1,k-1,k+1,b_arx,a).eq.1) count_bx=count_bx+1
 					if(b_ary(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j,j,k-1,k+1,b_ary,a).eq.1) count_by=count_by+1
 					if(b_arz(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j-1,j+1,k,k,b_arz,a).eq.1) count_bz=count_bz+1
 				end do
 			end do
 		end do 
		
 		allocate(usc_wtr_bx(count_bx,9), usc_wtr_by(count_by,9), usc_wtr_bz(count_bz,9)) ! r= row 
 		allocate(usc_wtc_bx(count_bx,9), usc_wtc_by(count_by,9), usc_wtc_bz(count_bz,9)) ! c= column

 		count_bx=0 ; count_by=0; count_bz=0;
 		usc_wtr_bx=0.0_psn; usc_wtr_by=0.0_psn; usc_wtr_bz=0.0_psn; 
 		usc_wtc_bx=0.0_psn; usc_wtc_by=0.0_psn; usc_wtc_bz=0.0_psn; 
		
 		usc_norm_bx=1.0_psn; usc_norm_by=1.0_psn; usc_norm_bz=1.0_psn; 
		
#ifdef twoD		
 		do k=1,1
#else
         do k=2,mz-1
#endif		
 			do j=2,my-1
 				do i=2,mx-1

 					if(b_arx(i,j,k).gt.0 .and. unstableCell(i,i,j-1,j+1,k-1,k+1,b_arx,a).eq.1) then
 						count_bx=count_bx+1
 						call set_usc_wtc_bx(i,j,k,count_bx,a)
 						usc_norm_bx(i,j,k)=usc_sum_wt(usc_wtc_bx,count_bx,9)
 						call usc_wt_normalise(usc_wtc_bx,count_bx,9,usc_norm_bx(i,j,k))
 						usc_fdtd_bx(i,j,k)=count_bx
 					end if

 					if(b_ary(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j,j,k-1,k+1,b_ary,a).eq.1) then
 						count_by=count_by+1
 						call set_usc_wtc_by(i,j,k,count_by,a)
 						usc_norm_by(i,j,k)=usc_sum_wt(usc_wtc_by,count_by,9)
 						call usc_wt_normalise(usc_wtc_by,count_by,9,usc_norm_by(i,j,k))
 						usc_fdtd_by(i,j,k)=count_by
 					end if

 					if(b_arz(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j-1,j+1,k,k,b_arz,a).eq.1) then
 						count_bz=count_bz+1
 						call set_usc_wtc_bz(i,j,k,count_bz,a)
 						usc_norm_bz(i,j,k)=usc_sum_wt(usc_wtc_bz,count_bz,9)
 						call usc_wt_normalise(usc_wtc_bz,count_bz,9,usc_norm_bz(i,j,k))
 						usc_fdtd_bz(i,j,k)=count_bz
 					end if

 				end do
 			end do
 		end do 
		
 		count_bx=0 ; count_by=0; count_bz=0;
#ifdef twoD		
 		do k=1,1
#else
         do k=2,mz-1
#endif			
 			do j=2,my-1
 				do i=2,mx-1
 					if(b_arx(i,j,k).gt.0 .and. unstableCell(i,i,j-1,j+1,k-1,k+1,b_arx,a).eq.1) then
 						count_bx=count_bx+1
 						call set_usc_wtr_bx(i,j,k,count_bx,a,usc_norm_bx)
 						call usc_wt_normalise(usc_wtr_bx,count_bx,9,usc_eff_area_bx(i,j,k,count_bx))
 					end if

 					if(b_ary(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j,j,k-1,k+1,b_ary,a).eq.1) then
 						count_by=count_by+1
 						call set_usc_wtr_by(i,j,k,count_by,a,usc_norm_by)
 						call usc_wt_normalise(usc_wtr_by,count_by,9,usc_eff_area_by(i,j,k,count_by))
 					end if

 					if(b_arz(i,j,k).gt.0 .and. unstableCell(i-1,i+1,j-1,j+1,k,k,b_arz,a).eq.1) then
 						count_bz=count_bz+1
 						call set_usc_wtr_bz(i,j,k,count_bz,a,usc_norm_bz)
 						call usc_wt_normalise(usc_wtr_bz,count_bz,9,usc_eff_area_bz(i,j,k,count_bz))
 					end if
 				end do
 			end do
 		end do 
		
 		deallocate(usc_norm_bx, usc_norm_by, usc_norm_bz)

 	end subroutine InitWeight_USC
		
 	integer function unstableCell(i1,i2,j1,j2,k1,k2, arr, a)
 		integer :: i1,i2,j1,j2,k1,k2
 		real(psn), dimension(:,:,:) :: arr
 		real(psn) :: a
 		integer :: i,j,k
		
 		unstableCell=0
#ifdef twoD
        do k=1,1 
#else
        do k=k1,k2
#endif		
	 		do j=j1,j2
	 			do i=i1,i2
	 				if(arr(i,j,k).gt.0 .and. arr(i,j,k).lt.a) unstableCell=1
	 			end do 
	 		end do 
		end do
 	end function unstableCell
	
	
 	subroutine set_usc_wtr_bx(i,j,k,n,a,norm)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:,:) :: norm
 		real(psn) :: a

#ifdef twoD 		
 		usc_wtr_bx(n,1)=1.0_psn /norm(i,j,k)
 		if(e_lz(i,j,k)  .ne.0) usc_wtr_bx(n,2)=max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j,k)) /norm(i,j-1,k)
 		if(e_lz(i,j+1,k).ne.0) usc_wtr_bx(n,3)=max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j+1,k)) /norm(i,j+1,k)
#else		
 		usc_wtr_bx(n,1)=1.0_psn /norm(i,j,k)
 		if(e_ly(i,j,k)  .ne.0) usc_wtr_bx(n,2)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k)) /norm(i,j,k-1)
 		if(e_ly(i,j,k+1).ne.0) usc_wtr_bx(n,3)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k+1))/norm(i,j,k+1)
 		if(e_lz(i,j,k)  .ne.0) usc_wtr_bx(n,4)=max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j,k)) /norm(i,j-1,k)
 		if(e_lz(i,j+1,k).ne.0) usc_wtr_bx(n,5)=max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j+1,k))/norm(i,j+1,k)
 		if(e_ly(i,j,k)  .ne.0 .and. e_lz(i,j,k)  .ne.0) usc_wtr_bx(n,6)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k)) * max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j,k))/norm(i,j-1,k-1)
 		if(e_ly(i,j,k+1).ne.0 .and. e_lz(i,j,k)  .ne.0) usc_wtr_bx(n,7)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k+1)) * max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j,k))/norm(i,j-1,k+1)
 		if(e_ly(i,j,k)  .ne.0 .and. e_lz(i,j+1,k).ne.0) usc_wtr_bx(n,8)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k)) * max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j+1,k))/norm(i,j+1,k-1) 
 		if(e_ly(i,j,k+1).ne.0 .and. e_lz(i,j+1,k).ne.0) usc_wtr_bx(n,9)=max(0.0_psn,a-b_arx(i,j,k)/e_ly(i,j,k+1)) * max(0.0_psn,a-b_arx(i,j,k)/e_lz(i,j+1,k))/norm(i,j+1,k+1)
#endif		
 	end subroutine set_usc_wtr_bx
	
 	subroutine set_usc_wtc_bx(i,j,k,n,a)
 		integer :: i,j,k,n
 		real(psn) :: a
#ifdef twoD		
 		usc_wtc_bx(n,1)=1.0_psn
 		if(e_lz(i,j,k)  .ne.0) usc_wtc_bx(n,2)=max(0.0_psn,a-b_arx(i,j-1,k)/e_lz(i,j,k)) 
 		if(e_lz(i,j+1,k).ne.0) usc_wtc_bx(n,3)=max(0.0_psn,a-b_arx(i,j+1,k)/e_lz(i,j+1,k))
#else 		
 		usc_wtc_bx(n,1)=1.0_psn
 		if(e_ly(i,j,k)  .ne.0) usc_wtc_bx(n,2)=max(0.0_psn,a-b_arx(i,j,k-1)/e_ly(i,j,k)) 
 		if(e_ly(i,j,k+1).ne.0) usc_wtc_bx(n,3)=max(0.0_psn,a-b_arx(i,j,k+1)/e_ly(i,j,k+1))
 		if(e_lz(i,j,k)  .ne.0) usc_wtc_bx(n,4)=max(0.0_psn,a-b_arx(i,j-1,k)/e_lz(i,j,k)) 
 		if(e_lz(i,j+1,k).ne.0) usc_wtc_bx(n,5)=max(0.0_psn,a-b_arx(i,j+1,k)/e_lz(i,j+1,k))
 		if(e_lz(i,j,k-1)  .ne.0 .and. e_ly(i,j-1,k)  .ne.0) usc_wtc_bx(n,6)=max(0.0_psn,a-b_arx(i,j-1,k-1)/e_lz(i,j,k-1)) * max(0.0_psn,a-b_arx(i,j-1,k-1)/e_ly(i,j-1,k)) 
 		if(e_lz(i,j,k+1)  .ne.0 .and. e_ly(i,j-1,k+1).ne.0) usc_wtc_bx(n,7)=max(0.0_psn,a-b_arx(i,j-1,k+1)/e_lz(i,j,k+1)) * max(0.0_psn,a-b_arx(i,j-1,k+1)/e_ly(i,j-1,k+1)) 
 		if(e_lz(i,j+1,k-1).ne.0 .and. e_ly(i,j+1,k)  .ne.0) usc_wtc_bx(n,8)=max(0.0_psn,a-b_arx(i,j+1,k-1)/e_lz(i,j+1,k-1)) * max(0.0_psn,a-b_arx(i,j+1,k-1)/e_ly(i,j+1,k))
 		if(e_lz(i,j+1,k+1).ne.0 .and. e_ly(i,j+1,k+1).ne.0) usc_wtc_bx(n,9)=max(0.0_psn,a-b_arx(i,j+1,k+1)/e_lz(i,j+1,k+1)) * max(0.0_psn,a-b_arx(i,j+1,k+1)/e_ly(i,j+1,k+1))
#endif		
 	end subroutine set_usc_wtc_bx
	
 	function usc_eff_area_bx(i,j,k,n)
 		integer :: i,j,k,n
 		real(psn) :: usc_eff_area_bx

#ifdef twoD		
 		usc_eff_area_bx = usc_wtr_bx(n,1)*b_arx(i,j,k)+&
 		                  usc_wtr_bx(n,2)*b_arx(i,j-1,k)+&
 						  usc_wtr_bx(n,3)*b_arx(i,j+1,k)
#else						  
   		usc_eff_area_bx = usc_wtr_bx(n,1)*b_arx(i,j,k)+&
   		                  usc_wtr_bx(n,2)*b_arx(i,j,k-1)+&
   						  usc_wtr_bx(n,3)*b_arx(i,j,k+1)+&
   						  usc_wtr_bx(n,4)*b_arx(i,j-1,k)+&
   						  usc_wtr_bx(n,5)*b_arx(i,j+1,k)+&
   						  usc_wtr_bx(n,6)*b_arx(i,j-1,k-1)+&
   						  usc_wtr_bx(n,7)*b_arx(i,j-1,k+1)+&
   						  usc_wtr_bx(n,8)*b_arx(i,j+1,k-1)+&
   						  usc_wtr_bx(n,9)*b_arx(i,j+1,k+1)
#endif						  
 	end function usc_eff_area_bx
	
 	function usc_gather_dbx(i,j,k,n,wt,arr)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:) :: wt
 		real(psn), dimension(:,:,:) :: arr
 		real(psn) :: usc_gather_dbx
#ifdef twoD		
 		usc_gather_dbx = wt(n,1)*arr(i,j,k) + wt(n,2)*arr(i,j-1,k) + wt(n,3)*arr(i,j+1,k)
#else		
 		usc_gather_dbx=wt(n,1)*arr(i,j,k) + wt(n,2)*arr(i,j,k-1) + wt(n,3)*arr(i,j,k+1) + wt(n,4)*arr(i,j-1,k) + wt(n,5)*arr(i,j+1,k)+&
 		               wt(n,6)*arr(i,j-1,k-1) + wt(n,7)*arr(i,j-1,k+1) + wt(n,8)*arr(i,j+1,k-1) + wt(n,9)*arr(i,j+1,k+1)
#endif	
 	end function usc_gather_dbx
	
	
 	subroutine set_usc_wtr_by(i,j,k,n,a,norm)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:,:) :: norm
 		real(psn) :: a

#ifdef twoD		
 		usc_wtr_by(n,1)=1.0_psn /norm(i,j,k)
 		if(e_lz(i,j,k)  .ne.0) usc_wtr_by(n,2)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i,j,k)) /norm(i-1,j,k)
 		if(e_lz(i+1,j,k).ne.0) usc_wtr_by(n,3)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i+1,j,k)) /norm(i+1,j,k)
#else		
 		usc_wtr_by(n,1)=1.0_psn /norm(i,j,k)
 		if(e_lz(i,j,k)  .ne.0) usc_wtr_by(n,2)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i,j,k)) /norm(i-1,j,k)
 		if(e_lz(i+1,j,k).ne.0) usc_wtr_by(n,3)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i+1,j,k))/norm(i+1,j,k)
 		if(e_lx(i,j,k)  .ne.0) usc_wtr_by(n,4)=max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k)) /norm(i,j,k-1)
 		if(e_lx(i,j,k+1).ne.0) usc_wtr_by(n,5)=max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k+1))/norm(i,j,k+1)
 		if(e_lz(i,j,k)  .ne.0 .and. e_lx(i,j,k)  .ne.0) usc_wtr_by(n,6)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i,j,k)) * max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k))/norm(i-1,j,k-1)
 		if(e_lz(i+1,j,k).ne.0 .and. e_lx(i,j,k)  .ne.0) usc_wtr_by(n,7)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i+1,j,k)) * max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k))/norm(i+1,j,k-1)
 		if(e_lz(i,j,k)  .ne.0 .and. e_lx(i,j,k+1).ne.0) usc_wtr_by(n,8)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i,j,k)) * max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k+1))/norm(i-1,j,k+1) 
 		if(e_lz(i+1,j,k).ne.0 .and. e_lx(i,j,k+1).ne.0) usc_wtr_by(n,9)=max(0.0_psn,a-b_ary(i,j,k)/e_lz(i+1,j,k)) * max(0.0_psn,a-b_ary(i,j,k)/e_lx(i,j,k+1))/norm(i+1,j,k+1)
#endif		
 	end subroutine set_usc_wtr_by
	
 	subroutine set_usc_wtc_by(i,j,k,n,a)
 		integer :: i,j,k,n
 		real(psn) :: a

#ifdef twoD 		
 		usc_wtc_by(n,1)=1.0_psn
 		if(e_lz(i,j,k)  .ne.0) usc_wtc_by(n,2)=max(0.0_psn,a-b_ary(i-1,j,k)/e_lz(i,j,k)) 
 		if(e_lz(i+1,j,k).ne.0) usc_wtc_by(n,3)=max(0.0_psn,a-b_ary(i+1,j,k)/e_lz(i+1,j,k))
#else				
 		usc_wtc_by(n,1)=1.0_psn
 		if(e_lz(i,j,k)  .ne.0) usc_wtc_by(n,2)=max(0.0_psn,a-b_ary(i-1,j,k)/e_lz(i,j,k)) 
 		if(e_lz(i+1,j,k).ne.0) usc_wtc_by(n,3)=max(0.0_psn,a-b_ary(i+1,j,k)/e_lz(i+1,j,k))
 		if(e_lx(i,j,k)  .ne.0) usc_wtc_by(n,4)=max(0.0_psn,a-b_ary(i,j,k-1)/e_lx(i,j,k)) 
 		if(e_lx(i,j,k+1).ne.0) usc_wtc_by(n,5)=max(0.0_psn,a-b_ary(i,j,k+1)/e_lx(i,j,k+1))
 		if(e_lx(i-1,j,k)  .ne.0 .and. e_lz(i,j,k-1)  .ne.0) usc_wtc_by(n,6)=max(0.0_psn,a-b_ary(i-1,j,k-1)/e_lx(i-1,j,k)) * max(0.0_psn,a-b_ary(i-1,j,k-1)/e_lz(i,j,k-1)) 
 		if(e_lx(i+1,j,k)  .ne.0 .and. e_lz(i+1,j,k-1).ne.0) usc_wtc_by(n,7)=max(0.0_psn,a-b_ary(i+1,j,k-1)/e_lx(i+1,j,k)) * max(0.0_psn,a-b_ary(i+1,j,k-1)/e_lz(i+1,j,k-1)) 
 		if(e_lx(i-1,j,k+1).ne.0 .and. e_lz(i,j,k+1)  .ne.0) usc_wtc_by(n,8)=max(0.0_psn,a-b_ary(i-1,j,k+1)/e_lx(i-1,j,k+1)) * max(0.0_psn,a-b_ary(i-1,j,k+1)/e_lz(i,j,k+1))
 		if(e_lx(i+1,j,k+1).ne.0 .and. e_lz(i+1,j,k+1).ne.0) usc_wtc_by(n,9)=max(0.0_psn,a-b_ary(i+1,j,k+1)/e_lx(i+1,j,k+1)) * max(0.0_psn,a-b_ary(i+1,j,k+1)/e_lz(i+1,j,k+1))
#endif	
 	end subroutine set_usc_wtc_by
	
 	function usc_eff_area_by(i,j,k,n)
 		integer :: i,j,k,n
 		real(psn) :: usc_eff_area_by
#ifdef twoD		
 		usc_eff_area_by = usc_wtr_by(n,1)*b_ary(i,j,k)+&
 		                  usc_wtr_by(n,2)*b_ary(i-1,j,k)+&
 						  usc_wtr_by(n,3)*b_ary(i+1,j,k)
#else 
 		usc_eff_area_by = usc_wtr_by(n,1)*b_ary(i,j,k)+&
 		                  usc_wtr_by(n,2)*b_ary(i-1,j,k)+&
 						  usc_wtr_by(n,3)*b_ary(i+1,j,k)+&
 						  usc_wtr_by(n,4)*b_ary(i,j,k-1)+&
 						  usc_wtr_by(n,5)*b_ary(i,j,k+1)+&
 						  usc_wtr_by(n,6)*b_ary(i-1,j,k-1)+&
 						  usc_wtr_by(n,7)*b_ary(i+1,j,k-1)+&
 						  usc_wtr_by(n,8)*b_ary(i-1,j,k+1)+&
 						  usc_wtr_by(n,9)*b_ary(i+1,j,k+1)
#endif						  
 	end function usc_eff_area_by
	
 	function usc_gather_dby(i,j,k,n,wt,arr)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:) :: wt
 		real(psn), dimension(:,:,:) :: arr
 		real(psn) :: usc_gather_dby
#ifdef twoD		
 		usc_gather_dby = wt(n,1)*arr(i,j,k) + wt(n,2)*arr(i-1,j,k) + wt(n,3)*arr(i+1,j,k)
#else		
 		usc_gather_dby=wt(n,1)*arr(i,j,k) + wt(n,2)*arr(i-1,j,k) + wt(n,3)*arr(i+1,j,k) + wt(n,4)*arr(i,j,k-1) + wt(n,5)*arr(i,j,k+1)+&
 		               wt(n,6)*arr(i-1,j,k-1) + wt(n,7)*arr(i+1,j,k-1) + wt(n,8)*arr(i-1,j,k+1) + wt(n,9)*arr(i+1,j,k+1)
#endif
 	end function usc_gather_dby
	
 	subroutine set_usc_wtr_bz(i,j,k,n,a,norm)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:,:) :: norm
 		real(psn) :: a
		
 		usc_wtr_bz(n,1)=1.0_psn /norm(i,j,k)
 		if(e_ly(i,j,k)  .ne.0) usc_wtr_bz(n,2)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i,j,k)) /norm(i-1,j,k)
 		if(e_ly(i+1,j,k).ne.0) usc_wtr_bz(n,3)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i+1,j,k))/norm(i+1,j,k)
 		if(e_lx(i,j,k)  .ne.0) usc_wtr_bz(n,4)=max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j,k)) /norm(i,j-1,k)
 		if(e_lx(i,j+1,k).ne.0) usc_wtr_bz(n,5)=max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j+1,k))/norm(i,j+1,k)
 		if(e_ly(i,j,k)  .ne.0 .and. e_lx(i,j,k)  .ne.0) usc_wtr_bz(n,6)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i,j,k)) * max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j,k))/norm(i-1,j-1,k)
 		if(e_ly(i+1,j,k).ne.0 .and. e_lx(i,j,k)  .ne.0) usc_wtr_bz(n,7)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i+1,j,k)) * max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j,k))/norm(i+1,j-1,k)
 		if(e_ly(i,j,k)  .ne.0 .and. e_lx(i,j+1,k).ne.0) usc_wtr_bz(n,8)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i,j,k)) * max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j+1,k))/norm(i-1,j+1,k) 
 		if(e_ly(i+1,j,k).ne.0 .and. e_lx(i,j+1,k).ne.0) usc_wtr_bz(n,9)=max(0.0_psn,a-b_arz(i,j,k)/e_ly(i+1,j,k)) * max(0.0_psn,a-b_arz(i,j,k)/e_lx(i,j+1,k))/norm(i+1,j+1,k)
		
 	end subroutine set_usc_wtr_bz
	
 	subroutine set_usc_wtc_bz(i,j,k,n,a)
 		integer :: i,j,k,n
 		real(psn) :: a
		
 		usc_wtc_bz(n,1)=1.0_psn
 		if(e_ly(i,j,k)  .ne.0) usc_wtc_bz(n,2)=max(0.0_psn,a-b_arz(i-1,j,k)/e_ly(i,j,k)) 
 		if(e_ly(i+1,j,k).ne.0) usc_wtc_bz(n,3)=max(0.0_psn,a-b_arz(i+1,j,k)/e_ly(i+1,j,k))
 		if(e_lx(i,j,k)  .ne.0) usc_wtc_bz(n,4)=max(0.0_psn,a-b_arz(i,j-1,k)/e_lx(i,j,k)) 
 		if(e_lx(i,j+1,k).ne.0) usc_wtc_bz(n,5)=max(0.0_psn,a-b_arz(i,j+1,k)/e_lx(i,j+1,k))
 		if(e_lx(i-1,j,k)  .ne.0 .and. e_ly(i,j-1,k)  .ne.0) usc_wtc_bz(n,6)=max(0.0_psn,a-b_arz(i-1,j-1,k)/e_lx(i-1,j,k)) * max(0.0_psn,a-b_arz(i-1,j-1,k)/e_ly(i,j-1,k)) 
 		if(e_lx(i+1,j,k)  .ne.0 .and. e_ly(i+1,j-1,k).ne.0) usc_wtc_bz(n,7)=max(0.0_psn,a-b_arz(i+1,j-1,k)/e_lx(i+1,j,k)) * max(0.0_psn,a-b_arz(i+1,j-1,k)/e_ly(i+1,j-1,k)) 
 		if(e_lx(i-1,j+1,k).ne.0 .and. e_ly(i,j+1,k)  .ne.0) usc_wtc_bz(n,8)=max(0.0_psn,a-b_arz(i-1,j+1,k)/e_lx(i-1,j+1,k)) * max(0.0_psn,a-b_arz(i-1,j+1,k)/e_ly(i,j+1,k))
 		if(e_lx(i+1,j+1,k).ne.0 .and. e_ly(i+1,j+1,k).ne.0) usc_wtc_bz(n,9)=max(0.0_psn,a-b_arz(i+1,j+1,k)/e_lx(i+1,j+1,k)) * max(0.0_psn,a-b_arz(i+1,j+1,k)/e_ly(i+1,j+1,k))
		
 	end subroutine set_usc_wtc_bz
	
 	function usc_eff_area_bz(i,j,k,n)
 		integer :: i,j,k,n
 		real(psn) :: usc_eff_area_bz 
 		usc_eff_area_bz = usc_wtr_bz(n,1)*b_arz(i,j,k)+&
 		                  usc_wtr_bz(n,2)*b_arz(i-1,j,k)+&
 						  usc_wtr_bz(n,3)*b_arz(i+1,j,k)+&
 						  usc_wtr_bz(n,4)*b_arz(i,j-1,k)+&
 						  usc_wtr_bz(n,5)*b_arz(i,j+1,k)+&
 						  usc_wtr_bz(n,6)*b_arz(i-1,j-1,k)+&
 						  usc_wtr_bz(n,7)*b_arz(i+1,j-1,k)+&
 						  usc_wtr_bz(n,8)*b_arz(i-1,j+1,k)+&
 						  usc_wtr_bz(n,9)*b_arz(i+1,j+1,k)
 	end function usc_eff_area_bz
	
	
 	function usc_gather_dbz(i,j,k,n,wt,arr)
 		integer :: i,j,k,n
 		real(psn), dimension(:,:) :: wt
 		real(psn), dimension(:,:,:) :: arr
 		real(psn) :: usc_gather_dbz
		
 		usc_gather_dbz=wt(n,1)*arr(i,j,k) + wt(n,2)*arr(i-1,j,k) + wt(n,3)*arr(i+1,j,k) + wt(n,4)*arr(i,j-1,k) + wt(n,5)*arr(i,j+1,k)+&
 		               wt(n,6)*arr(i-1,j-1,k) + wt(n,7)*arr(i+1,j-1,k) + wt(n,8)*arr(i-1,j+1,k) + wt(n,9)*arr(i+1,j+1,k)        
 	end function usc_gather_dbz
	
		
 	subroutine usc_wt_normalise(arr,n,size,norm)
 		real(psn), dimension(:,:) :: arr
 		integer :: n,size,i
 		real(psn) :: norm
		
 		do i=1,size
 			arr(n,i)=arr(n,i)/norm
 		end do 
 	end subroutine usc_wt_normalise
	
 	function usc_sum_wt(arr,n,size)
 		real(psn), dimension(:,:) :: arr
 		integer :: n, size, i
 		real(psn) :: usc_sum_wt
		
 		usc_sum_wt=0.0_psn 
 		do i=1,size
 			usc_sum_wt=usc_sum_wt+arr(n,i)
 		end do 
     end function usc_sum_wt 
	 
	 ! 	!--------------------------------------------------------------
	 ! 	! update the edge length as described in Igor Zagorodnov et al, JCP 2007 (SC method)
	 ! 	!--------------------------------------------------------------
	 ! 	subroutine update_lxy_cyl_SC
	 ! 		integer :: i,j
	 !
	 ! 		do j=2,my
	 ! 			do i=2,mx
	 ! 				e_lx(i,j,1)=min(2.0_psn*min(b_arz(i,j,1),b_arz(i,j-1,1),b_ary(i,j,1)),e_lx(i,j,1))
	 ! 				e_ly(i,j,1)=min(2.0_psn*min(b_arz(i,j,1),b_arz(i-1,j,1),b_arx(i,j,1)),e_ly(i,j,1))
	 ! 				e_lz(i,j,1)=min(2.0_psn*min(b_arx(i,j,1),b_arx(i,j-1,1),b_ary(i,j,1),b_ary(i-1,j,1)),e_lz(i,j,1))
	 ! 			end do
	 ! 		end do
	 !
	 ! 	end subroutine update_lxy_cyl_SC
	
	 
end module usc_CFDTD