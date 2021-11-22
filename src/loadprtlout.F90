module loadprtlout
use parameters 
use vars
use memory
implicit none 

contains 
!-------------------------
!Note : Reshaping of transfer array size is not in the scan loop 
! For partial safety safety check for the usage elsewhere
!-------------------------
	 
#ifdef twoD	 
     subroutine LoadPrtlOutliers
		 implicit none 
         integer:: i,ind,off
		  
		  do off=0,used_prtl_arr_size-1,outp_arr_block_size
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif			  
          do i=1+off,min(off+outp_arr_block_size,used_prtl_arr_size)	  
               if(qp(i).eq.0) cycle
               if(xp(i).lt.xmin) then                    
                    if((yp(i).le.(ymax+xmin-xp(i))).and.(yp(i).ge.(ymin-xmin+xp(i)))) then  
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=lcross
lcross=lcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							                       
                         call LoadPrtl(loutp,loutp_size,ind+1,i)
                         call DeletePrtl(i)                               
                    end if               
               else if(xp(i).gt.xmax) then
                    if((yp(i).le.(ymax+xp(i)-xmax)).and.(yp(i).ge.(ymin+xmax-xp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=rcross
rcross=rcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							
						 call LoadPrtl(routp,routp_size,ind+1,i)
                         call DeletePrtl(i)                                                            
                     end if     
               end if
               if(yp(i).lt.ymin) then          
                    if((xp(i).gt.(xmin-ymin+yp(i))).and.(xp(i).lt.(xmax+ymin-yp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=bcross
bcross=bcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif						
						 call LoadPrtl(boutp,boutp_size,ind+1,i)
                         call DeletePrtl(i)
                    end if
               else if(yp(i).gt.ymax) then
                    if((xp(i).gt.(xmin-yp(i)+ymax)).and.(xp(i).lt.(xmax+yp(i)-ymax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=tcross
tcross=tcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							
						 call LoadPrtl(toutp,toutp_size,ind+1,i)
                         call DeletePrtl(i)
                    end if
               end if			   
         end do
		 call UpdateTransferOutSize
	     end do 
		 
     end subroutine LoadPrtlOutliers
	 
     subroutine LoadTestPrtlOutliers
		 implicit none 
          integer:: i,ind,off
		  
		  do off=0,used_test_prtl_arr_size-1,outp_arr_block_size
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif		  
          do i=1+off,min(off+outp_arr_block_size,used_test_prtl_arr_size)	
               if(qtp(i).eq.0) cycle
               if(xtp(i).lt.xmin) then                    
                    if((ytp(i).le.(ymax+xmin-xtp(i))).and.(ytp(i).ge.(ymin-xmin+xtp(i)))) then 
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=lcross
lcross=lcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif						                        
                         call LoadTestPrtl(loutp,loutp_size,ind+1,i)
                         call DeleteTestPrtl(i)                               
                    end if               
               else if(xtp(i).gt.xmax) then
                    if((ytp(i).le.(ymax+xtp(i)-xmax)).and.(ytp(i).ge.(ymin+xmax-xtp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=rcross
rcross=rcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif						
						 call LoadTestPrtl(routp,routp_size,ind+1,i)
                         call DeleteTestPrtl(i)                                                            
                end if     
               end if
               if(ytp(i).lt.ymin) then          
                    if((xtp(i).gt.(xmin-ymin+ytp(i))).and.(xtp(i).lt.(xmax+ymin-ytp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=bcross
bcross=bcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif						
						 call LoadTestPrtl(boutp,boutp_size,ind+1,i)
                         call DeleteTestPrtl(i)
                    end if
               else if(ytp(i).gt.ymax) then
                    if((xtp(i).gt.(xmin-ytp(i)+ymax)).and.(xtp(i).lt.(xmax+ytp(i)-ymax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=tcross
tcross=tcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							
						 call LoadTestPrtl(toutp,toutp_size,ind+1,i)
                         call DeleteTestPrtl(i)
                    end if
               end if
         end do
		 call UpdateTransferOutSize
	     end do
     end subroutine LoadTestPrtlOutliers
	 
	 
	 
	 
#else	 
     subroutine LoadPrtlOutliers
		 implicit none
         integer:: i,ind,off
		  
		 do off=0,used_prtl_arr_size-1,outp_arr_block_size 
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
         do i=1+off,min(off+outp_arr_block_size,used_prtl_arr_size)
               if(qp(i).eq.0) cycle
               if(xp(i).lt.xmin) then                    
                    if((yp(i).le.(ymax+xmin-xp(i))).and.(yp(i).ge.(ymin-xmin+xp(i)))) then     
                         if(zp(i).le.(zmax+xmin-xp(i)).and.(zp(i).ge.(zmin-xmin+xp(i)))) then 
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=lcross
lcross=lcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif						                     
                             call LoadPrtl(loutp,loutp_size,ind+1,i)
                             call DeletePrtl(i) 
                        end if                              
                    end if               
               else if(xp(i).gt.xmax) then
                    if((yp(i).le.(ymax+xp(i)-xmax)).and.(yp(i).ge.(ymin+xmax-xp(i)))) then
                         if((zp(i).le.(zmax+xp(i)-xmax)).and.(zp(i).ge.(zmin+xmax-xp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=rcross
rcross=rcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							 
							 call LoadPrtl(routp,routp_size,ind+1,i)
                             call DeletePrtl(i)     
                        end if                                                       
                    end if     
               end if
               if(yp(i).lt.ymin) then          
                    if((xp(i).gt.(xmin-ymin+yp(i))).and.(xp(i).lt.(xmax+ymin-yp(i)))) then
                         if((zp(i).ge.(zmin-ymin+yp(i))).and.(zp(i).le.(zmax+ymin-yp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=bcross
bcross=bcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif								 
							 call LoadPrtl(boutp,boutp_size,ind+1,i)
                             call DeletePrtl(i)
                        end if
                    end if
               else if(yp(i).gt.ymax) then
                    if((xp(i).gt.(xmin-yp(i)+ymax)).and.(xp(i).lt.(xmax+yp(i)-ymax))) then
                         if((zp(i).ge.(zmin-yp(i)+ymax)).and.(zp(i).le.(zmax+yp(i)-ymax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=tcross
tcross=tcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							 
							 call LoadPrtl(toutp,toutp_size,ind+1,i)
                             call DeletePrtl(i)
                        end if
                    end if
               end if
               
               
               if(zp(i).lt.zmin) then          
                    if((xp(i).gt.(xmin-zmin+zp(i))).and.(xp(i).lt.(xmax+zmin-zp(i)))) then
                         if((yp(i).gt.(ymin-zmin+zp(i))).and.(yp(i).lt.(ymax+zmin-zp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=dcross
dcross=dcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif								 
							 call LoadPrtl(doutp,doutp_size,ind+1,i)
                             call DeletePrtl(i)
                        end if
                    end if
               else if(zp(i).gt.zmax) then
                    if((xp(i).gt.(xmin-zp(i)+zmax)).and.(xp(i).lt.(xmax+zp(i)-zmax))) then
                         if((yp(i).gt.(ymin-zp(i)+zmax)).and.(yp(i).lt.(ymax+zp(i)-zmax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=ucross
ucross=ucross+1
#ifdef OPEN_MP
!$omp end atomic
#endif	
							 call LoadPrtl(uoutp,uoutp_size,ind+1,i)
                             call DeletePrtl(i)
                        end if
                    end if
               end if
         end do
		 call UpdateTransferOutSize
	     end do
     end subroutine LoadPrtlOutliers
	 
     subroutine LoadTestPrtlOutliers
		 implicit none
          integer:: i,ind,off
		  
		  do off=0,used_test_prtl_arr_size-1,outp_arr_block_size
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif	
          do i=1+off,min(off+outp_arr_block_size,used_test_prtl_arr_size)
               if(qtp(i).eq.0) cycle
               if(xtp(i).lt.xmin) then                    
                    if((ytp(i).le.(ymax+xmin-xtp(i))).and.(ytp(i).ge.(ymin-xmin+xtp(i)))) then     
                         if(ztp(i).le.(zmax+xmin-xtp(i)).and.(ztp(i).ge.(zmin-xmin+xtp(i)))) then 
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=lcross
lcross=lcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							                     
                             call LoadTestPrtl(loutp,loutp_size,ind+1,i)
                             call DeleteTestPrtl(i) 
                        end if                              
                    end if               
               else if(xtp(i).gt.xmax) then
                    if((ytp(i).le.(ymax+xtp(i)-xmax)).and.(ytp(i).ge.(ymin+xmax-xtp(i)))) then
                         if((ztp(i).le.(zmax+xtp(i)-xmax)).and.(ztp(i).ge.(zmin+xmax-xtp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=rcross
rcross=rcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							 
							 call LoadTestPrtl(routp,routp_size,ind+1,i)
                             call DeleteTestPrtl(i)     
                        end if                                                       
                    end if     
               end if
               if(ytp(i).lt.ymin) then          
                    if((xtp(i).gt.(xmin-ymin+ytp(i))).and.(xtp(i).lt.(xmax+ymin-ytp(i)))) then
                         if((ztp(i).ge.(zmin-ymin+ytp(i))).and.(ztp(i).le.(zmax+ymin-ytp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=bcross
bcross=bcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif	
							 call LoadTestPrtl(boutp,boutp_size,ind+1,i)
                             call DeleteTestPrtl(i)
                        end if
                    end if
               else if(ytp(i).gt.ymax) then
                    if((xtp(i).gt.(xmin-ytp(i)+ymax)).and.(xtp(i).lt.(xmax+ytp(i)-ymax))) then
                         if((ztp(i).ge.(zmin-ytp(i)+ymax)).and.(ztp(i).le.(zmax+ytp(i)-ymax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=tcross
tcross=tcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif	
							 call LoadTestPrtl(toutp,toutp_size,ind+1,i)
                             call DeleteTestPrtl(i)
                        end if
                    end if
               end if
               
               
               if(ztp(i).lt.zmin) then          
                    if((xtp(i).gt.(xmin-zmin+ztp(i))).and.(xtp(i).lt.(xmax+zmin-ztp(i)))) then
                         if((ytp(i).gt.(ymin-zmin+ztp(i))).and.(ytp(i).lt.(ymax+zmin-ztp(i)))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=dcross
dcross=dcross+1
#ifdef OPEN_MP
!$omp end atomic
#endif							 
							 call LoadTestPrtl(doutp,doutp_size,ind+1,i)
                             call DeleteTestPrtl(i)
                        end if
                    end if
               else if(ztp(i).gt.zmax) then
                    if((xtp(i).gt.(xmin-ztp(i)+zmax)).and.(xtp(i).lt.(xmax+ztp(i)-zmax))) then
                         if((ytp(i).gt.(ymin-ztp(i)+zmax)).and.(ytp(i).lt.(ymax+ztp(i)-zmax))) then
#ifdef OPEN_MP
!$omp atomic capture
#endif
ind=ucross
ucross=ucross+1
#ifdef OPEN_MP
!$omp end atomic
#endif								 
							 call LoadTestPrtl(uoutp,uoutp_size,ind+1,i)
                             call DeleteTestPrtl(i)
                        end if
                    end if
               end if
         end do
		 call UpdateTransferOutSize
	     end do 
     end subroutine LoadTestPrtlOutliers
#endif
end module loadprtlout