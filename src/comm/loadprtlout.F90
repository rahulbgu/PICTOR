module loadprtlout
	use parameters 
	use vars
	use mem_prtl
implicit none 

contains 
!-------------------------
!Note : Reshaping of transfer array size is not in the scan loop 
! For partial safety safety check for the usage elsewhere
!-------------------------
	 
     subroutine LoadPrtlOutliers(q,x,y,z,u,v,w,var1,flv,tag,min_ind,max_ind)
		 implicit none
		 real(psn), dimension(:) :: q,x,y,z,u,v,w,var1
		 integer, dimension(:) :: flv,tag
		 integer :: min_ind, max_ind
         integer:: i,off

		 do off=min_ind-1,max_ind-1,outp_arr_block_size 
		 
		 call check_comm_prtl_size(ngbr_send)
			 
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i)
#endif
         do i=1+off,min(off+outp_arr_block_size,max_ind)
               if(flv(i).eq.0) cycle
               if(x(i).lt.xmin) then                    
                    if((y(i).le.(ymax+xmin-x(i))).and.(y(i).ge.(ymin-xmin+x(i)))) then  
#ifndef twoD						   
                         if(z(i).le.(zmax+xmin-x(i)).and.(z(i).ge.(zmin-xmin+x(i)))) then 
#endif						 					                     
                             call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(1),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i) 
#ifndef twoD							 
                        end if
#endif						                              
                    end if               
               else if(x(i).gt.xmax) then
                    if((y(i).le.(ymax-xmax+x(i))).and.(y(i).ge.(ymin+xmax-x(i)))) then
#ifndef twoD						
                         if((z(i).le.(zmax-xmax+x(i))).and.(z(i).ge.(zmin+xmax-x(i)))) then
#endif
							 call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(2),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i)  
#ifndef twoD							    
                        end if 
#endif						                                                      
                    end if     
               end if
               if(y(i).lt.ymin) then          
                    if((y(i).lt.(ymin-xmin+x(i))).and.(y(i).lt.(ymin+xmax-x(i)))) then 
#ifndef twoD						
                         if((z(i).ge.(zmin-ymin+y(i))).and.(z(i).le.(zmax+ymin-y(i)))) then
#endif
							 call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(3),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i)
#ifndef twoD							 
                        end if
#endif						
                    end if
               else if(y(i).gt.ymax) then
                    if((y(i).gt.(ymax+xmin-x(i))).and.(y(i).gt.(ymax-xmax+x(i)))) then 
#ifndef twoD						
                         if((z(i).ge.(zmin+ymax-y(i))).and.(z(i).le.(zmax-ymax+y(i)))) then
#endif							 
							 call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(4),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i)
#ifndef twoD							 
                        end if
#endif						
                    end if
               end if
               
#ifndef twoD               
               if(z(i).lt.zmin) then  
                    if((z(i).lt.(zmin-xmin+x(i))).and.(z(i).lt.(zmin+xmax-x(i)))) then  
						 if((z(i).lt.(zmin-ymin+y(i))).and.(z(i).lt.(zmin+ymax-y(i)))) then	 				 
							
							 call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(5),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i)
                        end if
                    end if
               else if(z(i).gt.zmax) then
                    if(z(i).gt.(zmax+xmin-x(i)).and.(z(i).gt.(zmax-xmax+x(i)))) then 
                         if((z(i).gt.(zmax+ymin-y(i))).and.(z(i).gt.(zmax-ymax+y(i)))) then

							 call CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,ngbr_send(6),i)
                             call RemovePrtl(q,flv,tag,x,y,z,u,v,w,i)
                        end if
                    end if
               end if
#endif			   
         end do
		 
	     
		 
		 end do
		 
	
    

    end subroutine LoadPrtlOutliers

	recursive subroutine CopyTransferPrtl(q,x,y,z,u,v,w,var1,flv,tag,list,m)
	  implicit none 
	  real(psn), dimension(:) :: q,x,y,z,u,v,w,var1
	  integer, dimension(:) :: flv,tag
	  real(psn) :: pos
	  type(ngbr_list) :: list
	  integer :: m,n,k
	  
	  if(indepLBaxis.eq.0) pos = x(m)
	  if(indepLBaxis.eq.1) pos = y(m)
	  if(indepLBaxis.eq.2) pos = z(m)
	  
	  call SendPrtlNgbrInd(list%num_ngbrs,list%edges,pos,k)
	  
	  n = list%ngbr(k)%pcount +1 
	  list%ngbr(k)%pcount = n
	  
	   
	  list%ngbr(k)%p(n)%q=q(m)
	  list%ngbr(k)%p(n)%x=x(m) - list%ngbr(k)%xshift
	  list%ngbr(k)%p(n)%y=y(m) - list%ngbr(k)%yshift
	  list%ngbr(k)%p(n)%z=z(m) - list%ngbr(k)%zshift
	  list%ngbr(k)%p(n)%u=u(m)
	  list%ngbr(k)%p(n)%v=v(m)
	  list%ngbr(k)%p(n)%w=w(m)
	  list%ngbr(k)%p(n)%flv=flv(m) 
	  list%ngbr(k)%p(n)%tag=tag(m) 
	  list%ngbr(k)%p(n)%var1=var1(m)
	end subroutine CopyTransferPrtl
	
	subroutine SendPrtlNgbrInd(size,edges,pos,k)
		integer, intent(IN) :: size
		integer, dimension(size+1), intent(IN) :: edges
		real(psn), intent(IN) :: pos
		integer, intent(INOUT) :: k
		integer :: n
		
		k=1 !default
		do n=1,size
			
			if(pos.ge.edges(n) .and. pos.lt.edges(n+1)) then 
				k=n
				exit
			end if 
		end do
	end subroutine SendPrtlNgbrInd
	
	subroutine RemovePrtl(q,flv,tag,x,y,z,u,v,w,ind)
		real(psn), dimension(:) :: q,x,y,z,u,v,w
		integer, dimension(:) :: flv,tag
		integer :: ind
		q(ind)=0
		flv(ind)=0
		tag(ind)=0
		u(ind)=0
		v(ind)=0
		w(ind)=0
		x(ind)=xmin
		y(ind)=ymin
		z(ind)=zmin
	end subroutine RemovePrtl
	
	 subroutine check_comm_prtl_size(list)
		 type(ngbr_list), dimension(:) :: list
		 integer :: n, m 
 
		 do n=1, size(list)
			 do m=1,list(n)%num_ngbrs
				 if( list(n)%ngbr(m)%pcount+outp_arr_block_size .gt. size(list(n)%ngbr(m)%p) ) then 
					 call ResizeCommPrtl( list(n)%ngbr(m)%p , list(n)%ngbr(m)%pcount, int(1.1*size(list(n)%ngbr(m)%p)) + 2*outp_arr_block_size )
				 end if 
			 end do 
		 end do 
 
	 end subroutine check_comm_prtl_size

	 subroutine ResizeCommPrtl(p,used_size,new_size)
		 type(particle), dimension(:), allocatable :: p
		 type(particle), dimension(:), allocatable :: p_temp
		 integer :: used_size, new_size 
		 integer :: n
		 allocate(p_temp(new_size))
		 do n=1,used_size
			 p_temp(n)%q = p(n)%q
			 p_temp(n)%x = p(n)%x
			 p_temp(n)%y = p(n)%y
			 p_temp(n)%z = p(n)%z
			 p_temp(n)%u = p(n)%u
			 p_temp(n)%v = p(n)%v
			 p_temp(n)%w = p(n)%w
			 p_temp(n)%flv = p(n)%flv
			 p_temp(n)%tag = p(n)%tag
			 p_temp(n)%var1 = p(n)%var1
		 end do 	 
		 call move_alloc(p_temp,p)
	 end subroutine ResizeCommPrtl
	

end module loadprtlout