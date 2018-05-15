module memory
     use parameters
     use vars
     implicit none 
     
contains 
     recursive subroutine DeletePrtl(n) !free the slot of particle at i if particle is leaving this proc
          integer :: n
                    qp(n)=0.0_psn
                    xp(n)=xmin+0.5_psn
                    yp(n)=ymin+0.5_psn
                    zp(n)=zmin+0.5_psn
                    up(n)=0.0_psn
                    vp(n)=0.0_psn
                    wp(n)=0.0_psn
                    tagp(n)=0
                    flvp(n)=0
					var1p(n)=0.0_psn					
     end subroutine DeletePrtl 
     recursive subroutine DeleteTestPrtl(n) !free the slot of particle at i if particle is leaving this proc
          integer :: n
                    qtp(n)=0.0_psn
                    xtp(n)=xmin+0.5_psn
                    ytp(n)=ymin+0.5_psn
                    ztp(n)=zmin+0.5_psn
                    utp(n)=0.0_psn
                    vtp(n)=0.0_psn
                    wtp(n)=0.0_psn
                    tagtp(n)=0
                    flvtp(n)=0
					var1tp(n)=0.0_psn
     end subroutine DeleteTestPrtl 
	 
	recursive subroutine LoadPrtl(p,size,n,m)
	 	integer :: size,n,m 
	 	type(particle), dimension(size) :: p
	 		p(n)%q=qp(m)
	 		p(n)%x=xp(m)
	 		p(n)%y=yp(m)
	 		p(n)%z=zp(m)
	 		p(n)%u=up(m)
	 		p(n)%v=vp(m)
	 		p(n)%w=wp(m)
	 		p(n)%flv=flvp(m) 
	 		p(n)%tag=tagp(m)
			p(n)%var1=var1p(m)
	 end subroutine LoadPrtl

	 recursive subroutine LoadTestPrtl(p,size,n,m)
	 	integer :: size,n,m
	 	type(particle), dimension(size) :: p
	 		p(n)%q=qtp(m)
	 		p(n)%x=xtp(m)
	 		p(n)%y=ytp(m)
	 		p(n)%z=ztp(m)
	 		p(n)%u=utp(m)
	 		p(n)%v=vtp(m)
	 		p(n)%w=wtp(m)
	 		p(n)%flv=flvtp(m) 
	 		p(n)%tag=tagtp(m)
			p(n)%var1=var1tp(m)
	 end subroutine LoadTestPrtl
	 
     !the following subroutine is called when the main prtl arrays needs to be reallocated
     subroutine ReshapePrtlArr(new_size)
          integer, intent(IN) :: new_size
		  integer :: n
          allocate(qp_temp(new_size),xp_temp(new_size),yp_temp(new_size),zp_temp(new_size),up_temp(new_size),vp_temp(new_size),wp_temp(new_size),tagp_temp(new_size),flvp_temp(new_size),var1p_temp(new_size))
		  do n=1,used_prtl_arr_size
				    qp_temp(n)=qp(n)
					xp_temp(n)=xp(n)
					yp_temp(n)=yp(n)
					zp_temp(n)=zp(n)
					up_temp(n)=up(n)
					vp_temp(n)=vp(n)
					wp_temp(n)=wp(n)
					tagp_temp(n)=tagp(n)
					flvp_temp(n)=flvp(n)
					var1p_temp(n)=var1p(n)
		  end do
		  deallocate(qp,xp,yp,zp,up,vp,wp,tagp,flvp,var1p)
          call move_alloc(qp_temp,qp)
		  call move_alloc(xp_temp,xp)
		  call move_alloc(yp_temp,yp)
		  call move_alloc(zp_temp,zp)
		  call move_alloc(up_temp,up)
		  call move_alloc(vp_temp,vp)
		  call move_alloc(wp_temp,wp)
		  call move_alloc(tagp_temp,tagp)
		  call move_alloc(flvp_temp,flvp)
		  call move_alloc(var1p_temp,var1p)		            		            
          prtl_arr_size=new_size
		  
		  do n=used_prtl_arr_size+1,prtl_arr_size
			  qp(n)=0
			  flvp(n)=0
		  end do
     end subroutine ReshapePrtlArr
	 
     subroutine ReshapeTestPrtlArr(new_size)
          integer, intent(IN) :: new_size
		  integer :: n

          allocate(qp_temp(new_size),xp_temp(new_size),yp_temp(new_size),zp_temp(new_size),up_temp(new_size),vp_temp(new_size),wp_temp(new_size),tagp_temp(new_size),flvp_temp(new_size),var1p_temp(new_size))
          do n=1,used_test_prtl_arr_size
				    qp_temp(n)=qtp(n)
					xp_temp(n)=xtp(n)
					yp_temp(n)=ytp(n)
					zp_temp(n)=ztp(n)
					up_temp(n)=utp(n)
					vp_temp(n)=vtp(n)
					wp_temp(n)=wtp(n)
					tagp_temp(n)=tagtp(n)
					flvp_temp(n)=flvtp(n)
					var1p_temp(n)=var1tp(n)
		  end do
		  deallocate(qtp,xtp,ytp,ztp,utp,vtp,wtp,tagtp,flvtp,var1tp)		  
          call move_alloc(qp_temp,qtp)
		  call move_alloc(xp_temp,xtp)
		  call move_alloc(yp_temp,ytp)
		  call move_alloc(zp_temp,ztp)
		  call move_alloc(up_temp,utp)
		  call move_alloc(vp_temp,vtp)
		  call move_alloc(wp_temp,wtp)
		  call move_alloc(tagp_temp,tagtp)
		  call move_alloc(flvp_temp,flvtp)	
		  call move_alloc(var1p_temp,var1tp)		            	            
          test_prtl_arr_size=new_size
		  
		  do n=used_test_prtl_arr_size+1,test_prtl_arr_size
			  qtp(n)=0
			  flvtp(n)=0
		  end do		  
     end subroutine ReshapeTestPrtlArr
	 
	 
	 !-------------------------------------------------------------------------------------------------
	 ! Subroutine to help fill particle and test particle arrays
	 !-------------------------------------------------------------------------------------------------
	 ! use the following two subroutines to directly insert particles when all particles are inserted at the first time step  
	 subroutine InsertParticleAt(ind1,x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)
	      integer    :: ind1
	      real(psn)  :: x1,y1,z1
	      real(psn)  :: u1,v1,w1
	      real(psn)  :: q1
		  real(psn)  :: var1_this
	      integer    :: a1
	      integer    :: flv1
     
	      xp(ind1)=x1
	      yp(ind1)=y1
	      zp(ind1)=z1
	      up(ind1)=u1
	      vp(ind1)=v1
	      wp(ind1)=w1
	      qp(ind1)=q1
	      tagp(ind1)=a1
	      flvp(ind1)=flv1
		  var1p(ind1)=var1_this
	 end subroutine InsertParticleAt 

	 subroutine InsertTestParticleAt(ind1,x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)
	      integer    :: ind1
	      real(psn)  :: x1,y1,z1
	      real(psn)  :: u1,v1,w1
	      real(psn)  :: q1
		  real(psn)  :: var1_this
	      integer    :: a1
	      integer    :: flv1
     
	      xtp(ind1)=x1
	      ytp(ind1)=y1
	      ztp(ind1)=z1
	      utp(ind1)=u1
	      vtp(ind1)=v1
	      wtp(ind1)=w1
	      qtp(ind1)=q1
		  tagtp(ind1)=a1
	      flvtp(ind1)=flv1
		  var1tp(ind1)=var1_this
	 end subroutine InsertTestParticleAt 

	 !use the following two subroutines if particles are to be inserted at some random time steps
	 !New particles are inserted in a zone where ordering is not assumed 
	 subroutine InsertNewPrtl(x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)
	      real(psn)  :: x1,y1,z1
	      real(psn)  :: u1,v1,w1
	      real(psn)  :: q1
		  real(psn)  :: var1_this
	      integer    :: a1
	      integer     :: flv1
	 	  if(prtl_random_insert_index.gt.prtl_arr_size) call ReshapePrtlArr(int(1.1*prtl_arr_size+100)) ! make sure that the prtl array is large enough
	 	  call InsertParticleAt(prtl_random_insert_index,x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)	
		  used_prtl_arr_size=prtl_random_insert_index
	      prtl_random_insert_index=prtl_random_insert_index+1
	 	  np=np+1
	 end subroutine InsertNewPrtl

	 subroutine InsertNewTestPrtl(x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)
	      real(psn)  :: x1,y1,z1
	      real(psn)  :: u1,v1,w1
	      real(psn)  :: q1
		  real(psn)  :: var1_this
	      integer    :: a1
	      integer     :: flv1
	 	  if(test_prtl_random_insert_index.gt.test_prtl_arr_size) call ReshapeTestPrtlArr(int(1.1*test_prtl_arr_size+100)) !make sure that the prtl array is large enough
	 	  call InsertTestParticleAt(test_prtl_random_insert_index,x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)	
		  used_test_prtl_arr_size=test_prtl_random_insert_index
	      test_prtl_random_insert_index=test_prtl_random_insert_index+1
	 	  ntp=ntp+1
	 end subroutine InsertNewTestPrtl
	 
     subroutine ReshapeTransferOutArr(arr,curr_size,new_size)
          integer :: n,curr_size,new_size
          type(particle), dimension(:),allocatable :: arr 
          allocate(ptemp(new_size))
          do n=1,curr_size
              ptemp(n)=arr(n)
          end do 
		  deallocate(arr)
          call move_alloc(ptemp,arr)
          curr_size=new_size     
     end subroutine ReshapeTransferOutArr
     subroutine ReshapeTransferInArr(arr,curr_size,new_size)
          integer :: curr_size,new_size
          type(particle),dimension(:),allocatable,intent(inout) :: arr 
          deallocate(arr)
          allocate(arr(new_size))
          curr_size=new_size
     end subroutine ReshapeTransferInArr
     subroutine ReorderPrtl
          if(modulo(t,prtl_reorder_period).eq.0) then
             call ReorderPrtlArr
			 call ReorderTestPrtlArr
         end if
     end subroutine ReorderPrtl
     subroutine ReorderPrtlArr
		  implicit none 
          integer :: n,i,j,k,offset,pcount_this
          integer, dimension (mx,my,mz) :: pcount 
		  integer :: count, np_cpu
		  integer :: ind
          pcount=0
		  count=0
		  
		  if((size(SortedPrtlCountYZ,1).ne.my).or.(size(SortedPrtlCountYZ,2).ne.mz)) then 
			  deallocate(SortedPrtlCountYZ)
			  allocate(SortedPrtlCountYZ(my,mz))
		  end if 
		  SortedPrtlCountYZ=0
          
		  do n=1,used_prtl_arr_size
               if(qp(n).ne.0) then
                 i=floor(xp(n))
                 j=floor(yp(n))
#ifndef twoD                 
                 k=floor(zp(n))
#else
                 k=1
#endif
                 pcount(i,j,k)=pcount(i,j,k)+1
				 count=count+1
              end if
          end do
		  np_cpu=count
		  
		  
       ! keep in mind that some particles end up right at the boundary 
        offset=1
  	    count=0		  
#ifdef twoD
          do k=1,1
#else               
          do k=1,mz
#endif
               do j=1,my
                    do i=1,mx
                         pcount_this=pcount(i,j,k)
                         pcount(i,j,k)=offset
                         offset=offset+pcount_this
                         count=count+pcount_this
                    end do
					SortedPrtlCountYZ(j,k)=count					
               end do 
          end do
		            
          allocate(qp_temp(prtl_arr_size))
          allocate(xp_temp(prtl_arr_size))
          allocate(yp_temp(prtl_arr_size))
          allocate(zp_temp(prtl_arr_size))
          allocate(up_temp(prtl_arr_size))
          allocate(vp_temp(prtl_arr_size))
          allocate(wp_temp(prtl_arr_size))
          allocate(tagp_temp(prtl_arr_size))
          allocate(flvp_temp(prtl_arr_size))
          allocate(var1p_temp(prtl_arr_size))		  

#ifdef OPEN_MP
!$OMP PARALLEL DO PRIVATE(n,i,j,k,ind)
#endif		  
          do n=1,used_prtl_arr_size
               if(qp(n).ne.0) then
                 i=floor(xp(n))
                 j=floor(yp(n))
#ifndef twoD
                 k=floor(zp(n))
#else
                 k=1
#endif    
#ifdef OPEN_MP
!$omp atomic capture
#endif
                 ind=pcount(i,j,k)
                 pcount(i,j,k)=pcount(i,j,k)+1
#ifdef OPEN_MP
!$omp end atomic
#endif	       
                 qp_temp(ind)=qp(n)
                 xp_temp(ind)=xp(n)
                 yp_temp(ind)=yp(n)
                 zp_temp(ind)=zp(n)
                 up_temp(ind)=up(n)
                 vp_temp(ind)=vp(n)
                 wp_temp(ind)=wp(n)
                 tagp_temp(ind)=tagp(n)
                 flvp_temp(ind)=flvp(n)	
                 var1p_temp(ind)=var1p(n)				 				 			 				 
              end if
          end do
		  
		  do n=np_cpu+1,prtl_arr_size
			  qp_temp(n)=0
			  flvp_temp(n)=0
		  end do 
		  
		  deallocate(qp,xp,yp,zp,up,vp,wp,tagp,flvp,var1p)		 
		  call move_alloc(qp_temp,qp)
		  call move_alloc(xp_temp,xp)
          call move_alloc(yp_temp,yp)
          call move_alloc(zp_temp,zp)
          call move_alloc(up_temp,up)
          call move_alloc(vp_temp,vp)
          call move_alloc(wp_temp,wp)
          call move_alloc(tagp_temp,tagp)
          call move_alloc(flvp_temp,flvp)
          call move_alloc(var1p_temp,var1p)
		  
		  used_prtl_arr_size=np_cpu!update the actual size of occupied part of prtl array  
		  prtl_random_insert_index=np_cpu+1
end subroutine ReorderPrtlArr


subroutine ReorderTestPrtlArr
		  implicit none 
          integer :: n,i,j,k,offset,pcount_this
          integer, dimension (mx,my,mz) :: pcount 
		  integer :: count,ntp_cpu
          pcount=0
		  count=0		  
          do n=1,used_test_prtl_arr_size
               if(qtp(n).ne.0) then
                 i=floor(xtp(n))
                 j=floor(ytp(n))
#ifndef twoD                 
                 k=floor(ztp(n))
#else
                 k=1
#endif
				 pcount(i,j,k)=pcount(i,j,k)+1
				 count=count+1
              end if
          end do
		  ntp_cpu=count
       ! keep in mind that some particles end up right at the boundary 
        offset=1
#ifdef twoD
          do k=1,1
#else               
          do k=1,mz
#endif
               do j=1,my
                    do i=1,mx
                         pcount_this=pcount(i,j,k)
                         pcount(i,j,k)=offset
                         offset=offset+pcount_this
					end do
               end do 
          end do
		   
          allocate(qp_temp(test_prtl_arr_size))
          allocate(xp_temp(test_prtl_arr_size))
          allocate(yp_temp(test_prtl_arr_size))
          allocate(zp_temp(test_prtl_arr_size))
          allocate(up_temp(test_prtl_arr_size))
          allocate(vp_temp(test_prtl_arr_size))
          allocate(wp_temp(test_prtl_arr_size))
          allocate(tagp_temp(test_prtl_arr_size))
          allocate(flvp_temp(test_prtl_arr_size))
          allocate(var1p_temp(test_prtl_arr_size))
		  
		  
          do n=1,used_test_prtl_arr_size
               if(qtp(n).ne.0) then
                 i=floor(xtp(n))
                 j=floor(ytp(n))
#ifndef twoD
                 k=floor(ztp(n))
#else
                 k=1
#endif           
                 qp_temp(pcount(i,j,k))=qtp(n)
                 xp_temp(pcount(i,j,k))=xtp(n)
                 yp_temp(pcount(i,j,k))=ytp(n)
                 zp_temp(pcount(i,j,k))=ztp(n)
                 up_temp(pcount(i,j,k))=utp(n)
                 vp_temp(pcount(i,j,k))=vtp(n)
                 wp_temp(pcount(i,j,k))=wtp(n)
                 tagp_temp(pcount(i,j,k))=tagtp(n)
                 flvp_temp(pcount(i,j,k))=flvtp(n)	
                 var1p_temp(pcount(i,j,k))=var1tp(n)				 				 			 				 
                 pcount(i,j,k)=pcount(i,j,k)+1
              end if
          end do
		  
		  do n=ntp_cpu+1,test_prtl_arr_size
			  qp_temp(n)=0
			  flvp_temp(n)=0
		  end do 

		  deallocate(qtp,xtp,ytp,ztp,utp,vtp,wtp,tagtp,flvtp,var1tp)		 
		  call move_alloc(qp_temp,qtp)
		  call move_alloc(xp_temp,xtp)
          call move_alloc(yp_temp,ytp)
          call move_alloc(zp_temp,ztp)
          call move_alloc(up_temp,utp)
          call move_alloc(vp_temp,vtp)
          call move_alloc(wp_temp,wtp)
          call move_alloc(tagp_temp,tagtp)
          call move_alloc(flvp_temp,flvtp)
          call move_alloc(var1p_temp,var1tp)
		  used_test_prtl_arr_size=ntp_cpu !update the actual size of occupied part of prtl array  
		  test_prtl_random_insert_index=ntp_cpu+1
		  
end subroutine ReorderTestPrtlArr

end module memory 