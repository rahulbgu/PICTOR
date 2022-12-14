module movdep_routines
	use parameters
	use vars
contains 
	 
	 subroutine GatherVecEfld(VecFx,VecFy,VecFz,Fx,Fy,Fz)
#ifdef twoD		 
		 real(psn), dimension(4,mx*my) :: VecFx,VecFy,VecFz
#else		 
		 real(psn), dimension(8,mx*my*mz) :: VecFx,VecFy,VecFz
#endif
		 real(psn), dimension(mx,my,mz) :: Fx,Fy,Fz
		 integer :: i1,j1,k1,short_arr_ind
		 integer :: imin, imax, jmin, jmax, kmin, kmax

		 call PrtlDomainIndMinMax(imin,imax,jmin,jmax,kmin,kmax)
		 
#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
			do k1=kmin,kmax
#else 
            do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 					
				do j1=jmin,jmax
!$OMP SIMD
					do i1=imin,imax
					  short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
			  		  VecFx(1,short_arr_ind)=(Fx(i1-1,j1,k1)+Fx(i1,j1,k1))*0.5_psn
			  		  VecFx(2,short_arr_ind)=(Fx(i1,j1,k1)+Fx(i1+1,j1,k1))*0.5_psn
			  		  VecFx(3,short_arr_ind)=(Fx(i1-1,j1+1,k1)+Fx(i1,j1+1,k1))*0.5_psn
			  		  VecFx(4,short_arr_ind)=(Fx(i1,j1+1,k1)+Fx(i1+1,j1+1,k1))*0.5_psn
#ifndef twoD 		  
                      VecFx(5,short_arr_ind)=(Fx(i1-1,j1,k1+1)+Fx(i1,j1,k1+1))*0.5_psn
                      VecFx(6,short_arr_ind)=(Fx(i1,j1,k1+1)+Fx(i1+1,j1,k1+1))*0.5_psn
                      VecFx(7,short_arr_ind)=(Fx(i1-1,j1+1,k1+1)+Fx(i1,j1+1,k1+1))*0.5_psn
                      VecFx(8,short_arr_ind)=(Fx(i1,j1+1,k1+1)+Fx(i1+1,j1+1,k1+1))*0.5_psn
#endif 
					end do 
				end do    
			end do	
			
			
								
#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
			do k1=kmin,kmax
#else 
            do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 					
				do j1=jmin,jmax
!$OMP SIMD					
					do i1=imin,imax
				     short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
	                 VecFy(1,short_arr_ind)=(Fy(i1,j1-1,k1)+Fy(i1,j1,k1))*0.5_psn
	                 VecFy(2,short_arr_ind)=(Fy(i1+1,j1-1,k1)+Fy(i1+1,j1,k1))*0.5_psn
	                 VecFy(3,short_arr_ind)=(Fy(i1,j1,k1)+Fy(i1,j1+1,k1))*0.5_psn
	                 VecFy(4,short_arr_ind)=(Fy(i1+1,j1,k1)+Fy(i1+1,j1+1,k1))*0.5_psn
#ifndef twoD 		  
	                 VecFy(5,short_arr_ind)=(Fy(i1,j1-1,k1+1)+Fy(i1,j1,k1+1))*0.5_psn
	                 VecFy(6,short_arr_ind)=(Fy(i1+1,j1-1,k1+1)+Fy(i1+1,j1,k1+1))*0.5_psn
	                 VecFy(7,short_arr_ind)=(Fy(i1,j1,k1+1)+Fy(i1,j1+1,k1+1))*0.5_psn
	                 VecFy(8,short_arr_ind)=(Fy(i1+1,j1,k1+1)+Fy(i1+1,j1+1,k1+1))*0.5_psn
#endif		
					end do 
				end do    
			end do	
			
	
#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
			do k1=kmin,kmax
#else 
            do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 					
				do j1=jmin,jmax
!$OMP SIMD					
					do i1=imin,imax			
					  short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
#ifndef twoD
					  VecFz(1,short_arr_ind)=(Fz(i1,j1,k1-1)+Fz(i1,j1,k1))*0.5_psn
					  VecFz(2,short_arr_ind)=(Fz(i1+1,j1,k1-1)+Fz(i1+1,j1,k1))*0.5_psn
					  VecFz(3,short_arr_ind)=(Fz(i1,j1+1,k1-1)+Fz(i1,j1+1,k1))*0.5_psn
					  VecFz(4,short_arr_ind)=(Fz(i1+1,j1+1,k1-1)+Fz(i1+1,j1+1,k1))*0.5_psn	  
					  VecFz(5,short_arr_ind)=(Fz(i1,j1,k1)+Fz(i1,j1,k1+1))*0.5_psn
					  VecFz(6,short_arr_ind)=(Fz(i1+1,j1,k1)+Fz(i1+1,j1,k1+1))*0.5_psn
					  VecFz(7,short_arr_ind)=(Fz(i1,j1+1,k1)+Fz(i1,j1+1,k1+1))*0.5_psn
					  VecFz(8,short_arr_ind)=(Fz(i1+1,j1+1,k1)+Fz(i1+1,j1+1,k1+1))*0.5_psn
#else
					  VecFz(1,short_arr_ind)=Fz(i1,j1,k1)
					  VecFz(2,short_arr_ind)=Fz(i1+1,j1,k1)
					  VecFz(3,short_arr_ind)=Fz(i1,j1+1,k1)
					  VecFz(4,short_arr_ind)=Fz(i1+1,j1+1,k1)		 
#endif 			
	
					end do 
				end do    
			end do	
	 end subroutine GatherVecEfld
	 
	 subroutine GatherVecBfld
		 integer :: i1,j1,k1,short_arr_ind
         integer :: imin, imax, jmin, jmax, kmin, kmax
		 
		 call PrtlDomainIndMinMax(imin,imax,jmin,jmax,kmin,kmax)
		 
#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
		 do k1=kmin,kmax
#else 
         do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 			
             do j1=jmin,jmax
!$OMP SIMD
		       do i1=imin,imax	
				      short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1

#ifndef twoD 		  
					  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1-1)+Bx(i1,j1,k1-1)+Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.25_psn	+Bx_ext0
					  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1-1)+Bx(i1+1,j1,k1-1)+Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.25_psn+Bx_ext0  	  
					  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1-1)+Bx(i1,j1+1,k1-1)+Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.25_psn	+Bx_ext0
					  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1-1)+Bx(i1+1,j1+1,k1-1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.25_psn	+Bx_ext0
					  VecBx(5,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1)+Bx(i1,j1-1,k1+1)+Bx(i1,j1,k1+1))*0.25_psn	+Bx_ext0
					  VecBx(6,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1)+Bx(i1+1,j1-1,k1+1)+Bx(i1+1,j1,k1+1))*0.25_psn +Bx_ext0 	  
					  VecBx(7,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1)+Bx(i1,j1,k1+1)+Bx(i1,j1+1,k1+1))*0.25_psn	+Bx_ext0
					  VecBx(8,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1)+Bx(i1+1,j1,k1+1)+Bx(i1+1,j1+1,k1+1))*0.25_psn	+Bx_ext0	
#else
					  VecBx(1,short_arr_ind)=(Bx(i1,j1-1,k1)+Bx(i1,j1,k1))*0.5_psn	+Bx_ext0
					  VecBx(2,short_arr_ind)=(Bx(i1+1,j1-1,k1)+Bx(i1+1,j1,k1))*0.5_psn +Bx_ext0 	  
					  VecBx(3,short_arr_ind)=(Bx(i1,j1,k1)+Bx(i1,j1+1,k1))*0.5_psn	+Bx_ext0
					  VecBx(4,short_arr_ind)=(Bx(i1+1,j1,k1)+Bx(i1+1,j1+1,k1))*0.5_psn+Bx_ext0		
#endif 
                  end do 
				 end do 
			end do 


#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
			do k1=kmin,kmax
#else 
            do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 			
             do j1=jmin,jmax	
!$OMP SIMD				 	
		       do i1=imin,imax	
				      short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
#ifndef twoD 		  
					  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1-1)+By(i1-1,j1,k1)+By(i1,j1,k1-1)+By(i1,j1,k1))*0.25_psn +By_ext0
					  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1-1)+By(i1,j1,k1)+By(i1+1,j1,k1-1)+By(i1+1,j1,k1))*0.25_psn +By_ext0
					  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1-1)+By(i1-1,j1+1,k1)+By(i1,j1+1,k1-1)+By(i1,j1+1,k1))*0.25_psn +By_ext0
					  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1-1)+By(i1,j1+1,k1)+By(i1+1,j1+1,k1-1)+By(i1+1,j1+1,k1))*0.25_psn+By_ext0
					  	  VecBy(5,short_arr_ind)=(By(i1-1,j1,k1)+By(i1-1,j1,k1+1)+By(i1,j1,k1)+By(i1,j1,k1+1))*0.25_psn+By_ext0
					  	  VecBy(6,short_arr_ind)=(By(i1,j1,k1)+By(i1,j1,k1+1)+By(i1+1,j1,k1)+By(i1+1,j1,k1+1))*0.25_psn+By_ext0
					  	  VecBy(7,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1-1,j1+1,k1+1)+By(i1,j1+1,k1)+By(i1,j1+1,k1+1))*0.25_psn+By_ext0
					  	  VecBy(8,short_arr_ind)=(By(i1,j1+1,k1)+By(i1,j1+1,k1+1)+By(i1+1,j1+1,k1)+By(i1+1,j1+1,k1+1))*0.25_psn+By_ext0
#else
					  	  VecBy(1,short_arr_ind)=(By(i1-1,j1,k1)+By(i1,j1,k1))*0.5_psn+By_ext0
					  	  VecBy(2,short_arr_ind)=(By(i1,j1,k1)+By(i1+1,j1,k1))*0.5_psn+By_ext0
					  	  VecBy(3,short_arr_ind)=(By(i1-1,j1+1,k1)+By(i1,j1+1,k1))*0.5_psn+By_ext0
					  	  VecBy(4,short_arr_ind)=(By(i1,j1+1,k1)+By(i1+1,j1+1,k1))*0.5_psn+By_ext0
#endif 

                  end do 
				 end do 
			end do 
			

#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
			do k1=kmin,kmax
#else 
            do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 			
             do j1=jmin,jmax
!$OMP SIMD				 
		       do i1=imin,imax	
				      short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1

					  VecBz(1,short_arr_ind)=(Bz(i1-1,j1-1,k1)+Bz(i1-1,j1,k1)+Bz(i1,j1-1,k1)+Bz(i1,j1,k1))*0.25_psn+Bz_ext0
					  VecBz(2,short_arr_ind)=(Bz(i1,j1-1,k1)+Bz(i1,j1,k1)+Bz(i1+1,j1-1,k1)+Bz(i1+1,j1,k1))*0.25_psn+Bz_ext0
					  VecBz(3,short_arr_ind)=(Bz(i1-1,j1,k1)+Bz(i1-1,j1+1,k1)+Bz(i1,j1,k1)+Bz(i1,j1+1,k1))*0.25_psn+Bz_ext0
					  VecBz(4,short_arr_ind)=(Bz(i1,j1,k1)+Bz(i1,j1+1,k1)+Bz(i1+1,j1,k1)+Bz(i1+1,j1+1,k1))*0.25_psn+Bz_ext0
#ifndef twoD 		  
					  VecBz(5,short_arr_ind)=(Bz(i1-1,j1-1,k1+1)+Bz(i1-1,j1,k1+1)+Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1))*0.25_psn+Bz_ext0
					  VecBz(6,short_arr_ind)=(Bz(i1,j1-1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1+1,j1-1,k1+1)+Bz(i1+1,j1,k1+1))*0.25_psn+Bz_ext0
					  VecBz(7,short_arr_ind)=(Bz(i1-1,j1,k1+1)+Bz(i1-1,j1+1,k1+1)+Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1))*0.25_psn+Bz_ext0
					  VecBz(8,short_arr_ind)=(Bz(i1,j1,k1+1)+Bz(i1,j1+1,k1+1)+Bz(i1+1,j1,k1+1)+Bz(i1+1,j1+1,k1+1))*0.25_psn+Bz_ext0
#endif 
                  end do 
				 end do 
			end do 	
		 
	 end subroutine GatherVecBfld
	 
	 subroutine PrtlDomainIndMinMax(imin,imax,jmin,jmax,kmin,kmax)
		 integer :: imin,imax,jmin,jmax,kmin,kmax
		 imin = 2; imax = mx-1; jmin = 2; jmax = my-1; kmin =2; kmax = mz-1;
#ifdef twoD
		 kmin =1; kmax =1; 
#endif		 
		 if(bc_face(1)%pos_prtl.gt.xborders(0))            imin = max(2,   floor(bc_face(1)%pos_prtl - xborders(procxind) + 3)-2)
		 if(bc_face(2)%pos_prtl.lt.xborders(nSubDomainsX)) imax = min(mx-1,ceiling(bc_face(2)%pos_prtl - xborders(procxind) + 3)+2)
		 if(bc_face(3)%pos_prtl.gt.yborders(0))            jmin = max(2,   floor(bc_face(3)%pos_prtl - yborders(procyind) + 3)-2)
		 if(bc_face(4)%pos_prtl.lt.yborders(nSubDomainsY)) jmax = min(my-1,ceiling(bc_face(4)%pos_prtl - yborders(procyind) + 3)+2)
#ifndef twoD
         if(bc_face(5)%pos_prtl.gt.zborders(0))            kmin = max(2,   floor(bc_face(5)%pos_prtl - zborders(proczind) + 3)-2)
         if(bc_face(6)%pos_prtl.lt.zborders(nSubDomainsZ)) kmax = min(mz-1,ceiling(bc_face(6)%pos_prtl - zborders(proczind) + 3)+2)
#endif 		 
	 end subroutine PrtlDomainIndMinMax
	 
	 
	 
	 subroutine ReduceCurrentMatrix
		 integer :: i1,j1,k1, short_arr_ind
		 

!-------------Collect Current From All threads ----------!
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(j1,k1)
#endif
do j1=1,size(VecJ,2)
do k1=2,Nthreads
VecJ(:,j1,1)=VecJ(:,j1,1)+VecJ(:,j1,k1)
end do 
end do 
!-------------------------------------------------------!




! Now deposit the current into the main array
#ifndef twoD
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,k1,short_arr_ind)
#endif				
	 do k1=1,mz-1
#else 
      do k1=1,1
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i1,j1,short_arr_ind)
#endif
#endif 
                do j1=1,my-1
!$OMP SIMD					   
		          do i1=1,mx-1
#ifndef twoD
				      short_arr_ind= (k1-1)*my*mx + (j1-1)*mx + i1
			          Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind,1)
					  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind,1)
					  Jx(i1,j1,k1+1)  =Jx(i1,j1,k1+1)  +VecJ(3,short_arr_ind,1)
					  Jx(i1,j1+1,k1+1)=Jx(i1,j1+1,k1+1)+VecJ(4,short_arr_ind,1)

			          Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(5,short_arr_ind,1)
					  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(6,short_arr_ind,1)
					  Jy(i1,j1,k1+1)  =Jy(i1,j1,k1+1)  +VecJ(7,short_arr_ind,1)
					  Jy(i1+1,j1,k1+1)=Jy(i1+1,j1,k1+1)+VecJ(8,short_arr_ind,1)

	                  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(9,short_arr_ind,1)
	                  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(10,short_arr_ind,1)
	                  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(11,short_arr_ind,1)
	                  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(12,short_arr_ind,1)
#else
					  short_arr_ind=  (j1-1)*mx + i1
					  Jx(i1,j1,k1)    =Jx(i1,j1,k1)    +VecJ(1,short_arr_ind,1)
					  Jx(i1,j1+1,k1)  =Jx(i1,j1+1,k1)  +VecJ(2,short_arr_ind,1)

					  Jy(i1,j1,k1)    =Jy(i1,j1,k1)    +VecJ(3,short_arr_ind,1)
					  Jy(i1+1,j1,k1)  =Jy(i1+1,j1,k1)  +VecJ(4,short_arr_ind,1)

					  Jz(i1,j1,k1)    =Jz(i1,j1,k1)    +VecJ(5,short_arr_ind,1)
					  Jz(i1+1,j1,k1)  =Jz(i1+1,j1,k1)  +VecJ(6,short_arr_ind,1)
					  Jz(i1,j1+1,k1)  =Jz(i1,j1+1,k1)  +VecJ(7,short_arr_ind,1)
					  Jz(i1+1,j1+1,k1)=Jz(i1+1,j1+1,k1)+VecJ(8,short_arr_ind,1)
#endif

				  end do
			   end do
		    end do

		 
	 end subroutine ReduceCurrentMatrix
end module movdep_routines