module comm_prtl_gpu
use parameters 
use vars
use var_gpu
implicit none 

contains 
!-------------------------
!Note : Reshaping of transfer array size is not in the scan loop 
! For partial safety safety check for the usage elsewhere
!-------------------------

recursive subroutine LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size1,outp,size2,m,n)
  implicit none 
  integer :: size1,size2
  real(psn), dimension(size1) :: q,x,y,z,u,v,w,var1
  integer, dimension(size1) :: flv,tag
  type(particle), dimension(size2) :: outp
  integer :: m,n
	outp(n)%q=q(m)
	outp(n)%x=x(m)
	outp(n)%y=y(m)
	outp(n)%z=z(m)
	outp(n)%u=u(m)
	outp(n)%v=v(m)
	outp(n)%w=w(m)
	outp(n)%flv=flv(m) 
	outp(n)%tag=tag(m)
	outp(n)%var1=var1(m)
end subroutine LoadTransferPrtlGPU	 

recursive subroutine UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size1,outp,size2,m,n)
  implicit none 
  integer :: size1,size2
  real(psn), dimension(size1) :: q,x,y,z,u,v,w,var1
  integer, dimension(size1) :: flv,tag
  type(particle), dimension(size2) :: outp
  integer, intent(IN) :: m,n
	q(n)=outp(m)%q
	x(n)=outp(m)%x
	y(n)=outp(m)%y
	z(n)=outp(m)%z
	u(n)=outp(m)%u
	v(n)=outp(m)%v
	w(n)=outp(m)%w
	flv(n)=outp(m)%flv
	tag(n)=outp(m)%tag
	var1(n)=outp(m)%var1
end subroutine UnloadTransferPrtlGPU 


#ifdef twoD	 
     subroutine LoadPrtlOutliersGPU(q,x,y,z,u,v,w,flv,var1,tag,size,ind_full,lindex,rindex,bindex,tindex,dindex,uindex,lc,rc,bc,tc,dc,uc,buff_size)
		 implicit none 
          integer:: ind_full,size,buff_size
		  real(psn), dimension(size) :: q,x,y,z,u,v,w,var1
		  integer, dimension(size) :: flv,tag
		  integer, dimension(buff_size,Nthreads) :: lindex,rindex,bindex,tindex,dindex,uindex
		  integer, dimension(Nthreads):: lc,rc,bc,tc,dc,uc
		  integer :: tid,i,j
		  tid=1
		  lc=0 
		  rc=0 
		  bc=0
		  tc=0
		  
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,tid)
#endif			  
          do i=1,ind_full
#ifdef OPEN_MP			  
			   tid=OMP_GET_THREAD_NUM()+1
#endif 			   			  
               if(x(i).lt.xmin) then                    
                    if((y(i).le.(ymax+xmin-x(i))).and.(y(i).ge.(ymin-xmin+x(i)))) then  
                        lc(tid)=lc(tid)+1
						lindex(lc(tid),tid)=i
                    end if               
               else if(x(i).gt.xmax) then
                    if((y(i).le.(ymax+x(i)-xmax)).and.(y(i).ge.(ymin+xmax-x(i)))) then
                        rc(tid)=rc(tid)+1
						rindex(rc(tid),tid)=i
                    end if     
               end if
               if(y(i).lt.ymin) then          
                    if((x(i).gt.(xmin-ymin+y(i))).and.(x(i).lt.(xmax+ymin-y(i)))) then						
                        bc(tid)=bc(tid)+1
						bindex(bc(tid),tid)=i
                    end if
               else if(y(i).gt.ymax) then
                    if((x(i).gt.(xmin-y(i)+ymax)).and.(x(i).lt.(xmax+y(i)-ymax))) then
                        tc(tid)=tc(tid)+1
						tindex(tc(tid),tid)=i
                    end if
               end if			   
         end do
		 
		 
		 
		 do j=1,Nthreads
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i)
#endif			 
			 do i=1,lc(j)
				 call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,loutp,loutp_size,lindex(i,j),i+lcross)
			 end do
			 lcross=lcross+lc(j) 
		 end do 
		 
		 do j=1,Nthreads
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i)
#endif			 
			 do i=1,rc(j)
				 call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,routp,routp_size,rindex(i,j),i+rcross)
			 end do
			 rcross=rcross+rc(j) 
		 end do 
		 
		 do j=1,Nthreads
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i)
#endif			 
			 do i=1,tc(j)
				 call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,toutp,toutp_size,tindex(i,j),i+tcross)
			 end do
			 tcross=tcross+tc(j) 
		 end do 
		 
		 do j=1,Nthreads
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i)
#endif			 
			 do i=1,bc(j)
				 call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,boutp,boutp_size,bindex(i,j),i+bcross)
			 end do
			 bcross=bcross+bc(j) 
		 end do 
		 
     end subroutine LoadPrtlOutliersGPU
	 	 
	 
	 
#else	 
     subroutine LoadPrtlOutliersGPU(q,x,y,z,u,v,w,flv,var1,tag,size,ind_full,lindex,rindex,bindex,tindex,dindex,uindex,lc,rc,bc,tc,dc,uc,buff_size)
		 implicit none
          integer:: ind_full,size,buff_size
		  real(psn), dimension(size) :: q,x,y,z,u,v,w,var1
		  integer, dimension(size) :: flv,tag
		  integer, dimension(buff_size,Nthreads) :: lindex,rindex,bindex,tindex,dindex,uindex
		  integer, dimension(Nthreads) :: lc,rc,bc,tc,dc,uc
		  integer :: tid,i,j
		  integer :: ind
		  
		  integer :: l0,r0,b0,t0,d0,u0 		  
		  
		  tid=1
		  lc=0 
		  rc=0 
		  bc=0
		  tc=0
		  dc=0
		  uc=0
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,tid)
#endif			  
          do i=1,ind_full
			  
#ifdef OPEN_MP			  
		      tid=OMP_GET_THREAD_NUM()+1
#endif 		
                 
		    
               if(x(i).lt.xmin) then                    
                    if(   ((x(i)+y(i)).le.(xmin+ymax)) .and. ((x(i)-y(i)).le.(xmin-ymin))   ) then     
                         if( ((z(i)+x(i)).le.(zmax+xmin)) .and. ((z(i)-x(i)).ge.(zmin-xmin))  ) then 
							lc(tid)=lc(tid)+1
	 						lindex(lc(tid),tid)=i
                        end if                              
                    end if               
               else if(x(i).gt.xmax) then
                    if(  ((x(i)-y(i)).ge.(xmax-ymax)) .and. ((x(i)+y(i)).ge.(xmax+ymin))    ) then
                         if(  ((z(i)-x(i)).le.(zmax-xmax)) .and. ((z(i)+x(i)).ge.(zmin+xmax)) ) then
	                        rc(tid)=rc(tid)+1
	 						rindex(rc(tid),tid)=i
                         end if                                                       
                    end if     
               end if
               if(y(i).lt.ymin) then          
                    if(   ((x(i)-y(i)).gt.(xmin-ymin)) .and. ((x(i)+y(i)).lt.(xmax+ymin))     ) then
                        if( ((y(i)-z(i)).le.(ymin-zmin)) .and. ((y(i)+z(i)).le.(ymin+zmax)) ) then
 	                        bc(tid)=bc(tid)+1
 	 						bindex(bc(tid),tid)=i
                        end if
                    end if
               else if(y(i).gt.ymax) then
                    if(   ((x(i)+y(i)).gt.(xmin+ymax)) .and. ((x(i)-y(i)).lt.(xmax-ymax))     ) then
                        if( ((y(i)+z(i)).ge.(zmin+ymax)) .and. ((y(i)-z(i)).ge.(ymax-zmax)) ) then
 	                        tc(tid)=tc(tid)+1
 	 						tindex(tc(tid),tid)=i
                        end if
                    end if
               end if
               
               
               if(z(i).lt.zmin) then          
                    if( ((z(i)-x(i)).lt.(zmin-xmin)) .and. ((z(i)+x(i)).lt.(zmin+xmax))     ) then
                         if( ((y(i)-z(i)).gt.(ymin-zmin)) .and. ((y(i)+z(i)).lt.(ymax+zmin))  ) then
  	                        dc(tid)=dc(tid)+1
  	 						dindex(dc(tid),tid)=i
                        end if
                    end if
               else if(z(i).gt.zmax) then
                    if(  ((z(i)+x(i)).gt.(zmax+xmin)) .and. ((z(i)-x(i)).gt.(zmax-xmax))    ) then
                         if( ((y(i)+z(i)).gt.(ymin+zmax)) .and. ((y(i)-z(i)).lt.(ymax-zmax))  ) then
   	                        uc(tid)=uc(tid)+1
   	 						uindex(uc(tid),tid)=i
                        end if
                    end if
               end if
			   
			   
			   
! 			   if((l0+r0+b0+t0+d0+u0).ne.1) then
! 				   print*,'A particle was send to two proc. and the cordinates of the particles is : ',x(i),y(i),z(i)
! 				   print*,'crossing bolleans are:',l0,r0,b0,t0,d0,u0
! 			   end if
			   
			   
		    
!                if(x(i).lt.xmin) then
!                     if((y(i).le.(xmin+ymax-x(i))).and.(y(i).ge.(ymin-xmin+x(i)))) then
!                          if(z(i).le.(zmax+xmin-x(i)).and.(z(i).ge.(zmin-xmin+x(i)))) then
! 							l0=1
! 							lc(tid)=lc(tid)+1
! 	 						lindex(lc(tid),tid)=i
!                         end if
!                     end if
!                else if(x(i).gt.xmax) then
!                     if((y(i).le.(ymax-xmax+x(i))).and.(y(i).ge.(ymin+xmax-x(i)))) then
!                          if((z(i).le.(zmax-xmax+x(i))).and.(z(i).ge.(zmin+xmax-x(i)))) then
! 							 r0=1
! 	                        rc(tid)=rc(tid)+1
! 	 						rindex(rc(tid),tid)=i
!                          end if
!                     end if
!                end if
!                if(y(i).lt.ymin) then
!                     if((x(i).gt.(xmin-ymin+y(i))).and.(x(i).lt.(xmax+ymin-y(i)))) then
!                          if((z(i).ge.(zmin-ymin+y(i))).and.(z(i).le.(zmax+ymin-y(i)))) then
! 							b0=1
!  	                        bc(tid)=bc(tid)+1
!  	 						bindex(bc(tid),tid)=i
!                         end if
!                     end if
!                else if(y(i).gt.ymax) then
!                     if((x(i).gt.(xmin+ymax-y(i))).and.(x(i).lt.(xmax-ymax+y(i)))) then
!                         if((z(i).ge.(zmin+ymax-y(i))).and.(z(i).le.(zmax-ymax+y(i)))) then
! 							t0=1
!  	                        tc(tid)=tc(tid)+1
!  	 						tindex(tc(tid),tid)=i
!                         end if
!                     end if
!                end if
!
!
!                if(z(i).lt.zmin) then
!                     if((x(i).gt.(xmin-zmin+z(i))).and.(x(i).lt.(xmax+zmin-z(i)))) then
!                          if((y(i).gt.(ymin-zmin+z(i))).and.(y(i).lt.(ymax+zmin-z(i)))) then
! 							d0=1
!   	                        dc(tid)=dc(tid)+1
!   	 						dindex(dc(tid),tid)=i
!                         end if
!                     end if
!                else if(z(i).gt.zmax) then
!                     if((x(i).gt.(xmin+zmax-z(i))).and.(x(i).lt.(xmax-zmax+z(i)))) then
!                          if((y(i).gt.(ymin+zmax-z(i))).and.(y(i).lt.(ymax-zmax+z(i)))) then
! 							 u0=1
!    	                        uc(tid)=uc(tid)+1
!    	 						uindex(uc(tid),tid)=i
!                         end if
!                     end if
!                end if
			   
			   
			   
         end do
		 
		 !print*,'Prtl cross count at CPU',lc,rc,bc,tc,dc,uc
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif			 
			 do i=1,lc(j)
				!call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,loutp,loutp_size,lindex(i,j),i+lcross)
				ind=lindex(i,j)
			 	loutp(i+lcross)%q=q(ind)
			 	loutp(i+lcross)%x=x(ind)
			 	loutp(i+lcross)%y=y(ind)
			 	loutp(i+lcross)%z=z(ind)
			 	loutp(i+lcross)%u=u(ind)
			 	loutp(i+lcross)%v=v(ind)
			 	loutp(i+lcross)%w=w(ind)
			 	loutp(i+lcross)%flv=flv(ind) 
			 	loutp(i+lcross)%tag=tag(ind)
			 	loutp(i+lcross)%var1=var1(ind)
			 end do
			 lcross=lcross+lc(j) 
		 end do 
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif		 
			 do i=1,rc(j)
				 !call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,routp,routp_size,rindex(i,j),i+rcross)
				    ind=rindex(i,j) 
	 			 	routp(i+rcross)%q=q(ind)
	 			 	routp(i+rcross)%x=x(ind)
	 			 	routp(i+rcross)%y=y(ind)
	 			 	routp(i+rcross)%z=z(ind)
	 			 	routp(i+rcross)%u=u(ind)
	 			 	routp(i+rcross)%v=v(ind)
	 			 	routp(i+rcross)%w=w(ind)
	 			 	routp(i+rcross)%flv=flv(ind) 
	 			 	routp(i+rcross)%tag=tag(ind)
	 			 	routp(i+rcross)%var1=var1(ind)
			 end do
			 rcross=rcross+rc(j) 
		 end do 
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif		 
			 do i=1,tc(j)
				!call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,toutp,toutp_size,tindex(i,j),i+tcross)
				ind=tindex(i,j)
 			 	toutp(i+tcross)%q=q(ind)
 			 	toutp(i+tcross)%x=x(ind)
 			 	toutp(i+tcross)%y=y(ind)
 			 	toutp(i+tcross)%z=z(ind)
 			 	toutp(i+tcross)%u=u(ind)
 			 	toutp(i+tcross)%v=v(ind)
 			 	toutp(i+tcross)%w=w(ind)
 			 	toutp(i+tcross)%flv=flv(ind) 
 			 	toutp(i+tcross)%tag=tag(ind)
 			 	toutp(i+tcross)%var1=var1(ind)
			 end do
			 tcross=tcross+tc(j) 
		 end do 
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif			 
			 do i=1,bc(j)
				!call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,boutp,boutp_size,bindex(i,j),i+bcross)
				ind=bindex(i,j)
 			 	boutp(i+bcross)%q=q(ind)
 			 	boutp(i+bcross)%x=x(ind)
 			 	boutp(i+bcross)%y=y(ind)
 			 	boutp(i+bcross)%z=z(ind)
 			 	boutp(i+bcross)%u=u(ind)
 			 	boutp(i+bcross)%v=v(ind)
 			 	boutp(i+bcross)%w=w(ind)
 			 	boutp(i+bcross)%flv=flv(ind) 
 			 	boutp(i+bcross)%tag=tag(ind)
 			 	boutp(i+bcross)%var1=var1(ind)
			 end do
			 bcross=bcross+bc(j) 
		 end do 
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif	 
			 do i=1,dc(j)
				!call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,doutp,doutp_size,dindex(i,j),i+dcross)
 			 	ind=dindex(i,j)
				doutp(i+dcross)%q=q(ind)
 			 	doutp(i+dcross)%x=x(ind)
 			 	doutp(i+dcross)%y=y(ind)
 			 	doutp(i+dcross)%z=z(ind)
 			 	doutp(i+dcross)%u=u(ind)
 			 	doutp(i+dcross)%v=v(ind)
 			 	doutp(i+dcross)%w=w(ind)
 			 	doutp(i+dcross)%flv=flv(ind) 
 			 	doutp(i+dcross)%tag=tag(ind)
 			 	doutp(i+dcross)%var1=var1(ind)
			 end do
			 dcross=dcross+dc(j) 
		 end do 
		 
		 do j=1,Nthreads
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i,ind)
! #endif			 
			 do i=1,uc(j)
				!call LoadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,uoutp,uoutp_size,uindex(i,j),i+ucross)
 			 	ind=uindex(i,j)
				uoutp(i+ucross)%q=q(ind)
 			 	uoutp(i+ucross)%x=x(ind)
 			 	uoutp(i+ucross)%y=y(ind)
 			 	uoutp(i+ucross)%z=z(ind)
 			 	uoutp(i+ucross)%u=u(ind)
 			 	uoutp(i+ucross)%v=v(ind)
 			 	uoutp(i+ucross)%w=w(ind)
 			 	uoutp(i+ucross)%flv=flv(ind) 
 			 	uoutp(i+ucross)%tag=tag(ind)
 			 	uoutp(i+ucross)%var1=var1(ind)
			 end do
			 ucross=ucross+uc(j) 
		 end do 
		 
     end subroutine LoadPrtlOutliersGPU
#endif

subroutine UnloadPrtlOutliersGPU(q,x,y,z,u,v,w,flv,var1,tag,size,li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf)
    integer:: size
    real(psn), dimension(size) :: q,x,y,z,u,v,w,var1
    integer, dimension(size) :: flv,tag
	integer :: li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf
	integer :: i
	integer :: ind
	integer :: off

!print*,'I:',li,ri,bi,ti,di,ui
!print*,'F:',lf,rf,bf,tf,df,uf

off=0
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=li,lf
		ind=i-li+1 + off
	   	q(ind)=linp(i)%q
	   	x(ind)=linp(i)%x-aint(linp(i)%x)+xmin
	   	y(ind)=linp(i)%y
	   	z(ind)=linp(i)%z
	   	u(ind)=linp(i)%u
	   	v(ind)=linp(i)%v
	   	w(ind)=linp(i)%w
	   	flv(ind)=linp(i)%flv
	   	tag(ind)=linp(i)%tag
	   	var1(ind)=linp(i)%var1
    end do
	
off=off+(lf-li+1)	
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=ri,rf
		ind=i-ri+1+off
	   	q(ind)=rinp(i)%q
	   	x(ind)=rinp(i)%x+xlen
	   	y(ind)=rinp(i)%y
	   	z(ind)=rinp(i)%z
	   	u(ind)=rinp(i)%u
	   	v(ind)=rinp(i)%v
	   	w(ind)=rinp(i)%w
	   	flv(ind)=rinp(i)%flv
	   	tag(ind)=rinp(i)%tag
	   	var1(ind)=rinp(i)%var1
    end do
	
off=off+(rf-ri+1)		
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=bi,bf
		ind=i-bi+1+off
	   	q(ind)=binp(i)%q
	   	x(ind)=binp(i)%x
	   	y(ind)=binp(i)%y-aint(binp(i)%y)+ymin
	   	z(ind)=binp(i)%z
	   	u(ind)=binp(i)%u
	   	v(ind)=binp(i)%v
	   	w(ind)=binp(i)%w
	   	flv(ind)=binp(i)%flv
	   	tag(ind)=binp(i)%tag
	   	var1(ind)=binp(i)%var1
    end do
	
off=off+(bf-bi+1)			
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=ti,tf
		ind=i-ti+1+off
	   	q(ind)=tinp(i)%q
	   	x(ind)=tinp(i)%x
	   	y(ind)=tinp(i)%y+ylen
	   	z(ind)=tinp(i)%z
	   	u(ind)=tinp(i)%u
	   	v(ind)=tinp(i)%v
	   	w(ind)=tinp(i)%w
	   	flv(ind)=tinp(i)%flv
	   	tag(ind)=tinp(i)%tag
	   	var1(ind)=tinp(i)%var1
    end do
	
	
#ifndef twoD

off=off+(tf-ti+1)			
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=di,df
		ind=i-di+1+off
	   	q(ind)=dinp(i)%q
	   	x(ind)=dinp(i)%x
	   	y(ind)=dinp(i)%y
	   	z(ind)=dinp(i)%z-aint(dinp(i)%z)+zmin
	   	u(ind)=dinp(i)%u
	   	v(ind)=dinp(i)%v
	   	w(ind)=dinp(i)%w
	   	flv(ind)=dinp(i)%flv
	   	tag(ind)=dinp(i)%tag
	   	var1(ind)=dinp(i)%var1
    end do

off=off+(df-di+1)				
#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(i,ind)
#endif
	do i=ui,uf
		ind=i-ui+1+off
	   	q(ind)=uinp(i)%q
	   	x(ind)=uinp(i)%x
	   	y(ind)=uinp(i)%y
	   	z(ind)=uinp(i)%z+zlen
	   	u(ind)=uinp(i)%u
	   	v(ind)=uinp(i)%v
	   	w(ind)=uinp(i)%w
	   	flv(ind)=uinp(i)%flv
	   	tag(ind)=uinp(i)%tag
	   	var1(ind)=uinp(i)%var1
    end do
#endif
end subroutine UnloadPrtlOutliersGPU


! subroutine UnloadPrtlOutliersGPU(q,x,y,z,u,v,w,flv,var1,tag,size,li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf)
!     integer:: size
!     real(psn), dimension(size) :: q,x,y,z,u,v,w,var1
!     integer, dimension(size) :: flv,tag
! 	integer :: li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf
! 	integer :: i
!
! print*,'I:',li,ri,bi,ti,di,ui
! print*,'F:',lf,rf,bf,tf,df,uf
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=li,lf
! 	   	q(i-li+1)=1.0
! 	   	x(i-li+1)=4
! 	   	y(i-li+1)=4
! 	   	z(i-li+1)=4
! 	   	u(i-li+1)=0.1
! 	   	v(i-li+1)=0.1
! 	   	w(i-li+1)=0.1
! 	   	flv(i-li+1)=1
! 	   	tag(i-li+1)=0
! 	   	var1(i-li+1)=0
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ri,rf
! 	   	q(i-ri+1)=1.0
! 	   	x(i-ri+1)=4
! 	   	y(i-ri+1)=4
! 	   	z(i-ri+1)=4
! 	   	u(i-ri+1)=0.1
! 	   	v(i-ri+1)=0.1
! 	   	w(i-ri+1)=0.1
! 	   	flv(i-ri+1)=1
! 	   	tag(i-ri+1)=0
! 	   	var1(i-ri+1)=0
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=bi,bf
! 	   	q(i-bi+1)=1.0
! 	   	x(i-bi+1)=4
! 	   	y(i-bi+1)=4
! 	   	z(i-bi+1)=4
! 	   	u(i-bi+1)=0.1
! 	   	v(i-bi+1)=0.1
! 	   	w(i-bi+1)=0.1
! 	   	flv(i-bi+1)=1
! 	   	tag(i-bi+1)=0
! 	   	var1(i-bi+1)=0
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ti,tf
! 	   	q(i-ti+1)=1.0
! 	   	x(i-ti+1)=4
! 	   	y(i-ti+1)=4
! 	   	z(i-ti+1)=4
! 	   	u(i-ti+1)=0.1
! 	   	v(i-ti+1)=0.1
! 	   	w(i-ti+1)=0.1
! 	   	flv(i-ti+1)=1
! 	   	tag(i-ti+1)=0
! 	   	var1(i-ti+1)=0
!     end do
! #ifndef twoD
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=di,df
! 	   	q(i-di+1)=1.0
! 	   	x(i-di+1)=4
! 	   	y(i-di+1)=4
! 	   	z(i-di+1)=4
! 	   	u(i-di+1)=0.1
! 	   	v(i-di+1)=0.1
! 	   	w(i-di+1)=0.1
! 	   	flv(i-di+1)=1
! 	   	tag(i-di+1)=0
! 	   	var1(i-di+1)=0
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ui,uf
! 	   	q(i-ui+1)=1.0
! 	   	x(i-ui+1)=4
! 	   	y(i-ui+1)=4
! 	   	z(i-ui+1)=4
! 	   	u(i-ui+1)=0.1
! 	   	v(i-ui+1)=0.1
! 	   	w(i-ui+1)=0.1
! 	   	flv(i-ui+1)=1
! 	   	tag(i-ui+1)=0
! 	   	var1(i-ui+1)=0
!     end do
! #endif
! end subroutine UnloadPrtlOutliersGPU



! subroutine UnloadPrtlOutliersGPU(q,x,y,z,u,v,w,flv,var1,tag,size,li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf)
!     integer:: size
!     real(psn), dimension(size) :: q,x,y,z,u,v,w,var1
!     integer, dimension(size) :: flv,tag
! 	integer :: li,lf,ri,rf,bi,bf,ti,tf,di,df,ui,uf
! 	integer :: i
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=li,lf
! 	   linp(i)%x=linp(i)%x-aint(linp(i)%x)+xmin
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,linp,linp_size,i,i-li+1)
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ri,rf
! 	   rinp(i)%x=rinp(i)%x+xlen
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,rinp,rinp_size,i,i-ri+1)
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=bi,bf
! 	   binp(i)%y=binp(i)%y-aint(binp(i)%y)+ymin
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,binp,binp_size,i,i-bi+1)
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ti,tf
! 	   tinp(i)%y=tinp(i)%y+ylen
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,tinp,tinp_size,i,i-ti+1)
!     end do
! #ifndef twoD
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=di,df
! 	   dinp(i)%z=dinp(i)%z-aint(dinp(i)%z)+zmin
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,dinp,dinp_size,i,i-di+1)
!     end do
! #ifdef OPEN_MP
! !$OMP PARALLEL DO PRIVATE(i)
! #endif
! 	do i=ui,uf
! 	   uinp(i)%z=uinp(i)%z+zlen
! 	   call UnloadTransferPrtlGPU(q,x,y,z,u,v,w,flv,var1,tag,size,uinp,uinp_size,i,i-ui+1)
!     end do
! #endif
! end subroutine UnloadPrtlOutliersGPU
	

end module comm_prtl_gpu