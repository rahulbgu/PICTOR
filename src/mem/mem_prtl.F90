module mem_prtl
     use parameters
     use vars
     implicit none 
     
contains 
     recursive subroutine DeletePrtl(n) !free the slot of particle at i if particle is leaving this proc
          integer :: n
                    qp(n)=0.0_psn
                    xp(n)=xmin+0.8_psn
                    yp(n)=ymin+0.8_psn
                    zp(n)=zmin+0.8_psn
                    up(n)=0.0_psn
                    vp(n)=0.0_psn
                    wp(n)=0.0_psn
                    tagp(n)=0
                    flvp(n)=0
					!var1p(n)=0.0_psn					
     end subroutine DeletePrtl 	 
	 
     !the following subroutine is called when the main prtl arrays needs to be resized
     subroutine ReshapePrtlArr(new_size,used_ind)
          integer, intent(IN) :: new_size
		  integer, optional :: used_ind
		  integer :: n,used_ind_this
		  
		  if(present(used_ind)) then
			  used_ind_this=used_ind 
		  else 
			  used_ind_this=used_prtl_arr_size
		  end if 

          allocate(qp_temp(new_size),xp_temp(new_size),yp_temp(new_size),zp_temp(new_size),up_temp(new_size),vp_temp(new_size),wp_temp(new_size),tagp_temp(new_size),flvp_temp(new_size),var1p_temp(new_size))
		  do n=1,used_ind_this
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
		  
		  do n=used_ind_this+1,prtl_arr_size
			  qp(n)=0
			  flvp(n)=0
		  end do
     end subroutine ReshapePrtlArr
		 
	 !-------------------------------------------------------------------------------------------------
	 ! Subroutine to help fill particle arrays
	 !-------------------------------------------------------------------------------------------------
	 ! use the following subroutine to directly insert particles into the particle arrays 
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

	 !use the following subroutine to insert an new particle into the particle arrays
	 subroutine InsertNewPrtl(x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)
	      real(psn)  :: x1,y1,z1
	      real(psn)  :: u1,v1,w1
	      real(psn)  :: q1
		  real(psn)  :: var1_this
	      integer    :: a1
	      integer     :: flv1
	 	  if(used_prtl_arr_size.ge.prtl_arr_size) call ReshapePrtlArr(int(1.1*prtl_arr_size+100)) ! make sure that the prtl array is large enough
	 	  call InsertParticleAt(used_prtl_arr_size+1,x1,y1,z1,u1,v1,w1,q1,a1,flv1,var1_this)	
		  used_prtl_arr_size=used_prtl_arr_size+1
	 	  np=np+1
	 end subroutine InsertNewPrtl
	 
	 subroutine InsertPrtls(p,xshift,yshift,zshift,nprtl)
		 type(particle), dimension(:) :: p
		 real(psn) :: xshift, yshift, zshift
		 integer :: nprtl, n 
		 
		 do n = 1,nprtl
			 qp(used_prtl_arr_size+n) = p(n)%q
			 xp(used_prtl_arr_size+n) = p(n)%x +xshift
			 yp(used_prtl_arr_size+n) = p(n)%y +yshift
			 zp(used_prtl_arr_size+n) = p(n)%z +zshift
			 up(used_prtl_arr_size+n) = p(n)%u
			 vp(used_prtl_arr_size+n) = p(n)%v
			 wp(used_prtl_arr_size+n) = p(n)%w
			 var1p(used_prtl_arr_size+n) = p(n)%var1
			 flvp(used_prtl_arr_size+n) = p(n)%flv
			 tagp(used_prtl_arr_size+n) = p(n)%tag
		 end do	
		 used_prtl_arr_size = used_prtl_arr_size+nprtl
		 np = np +nprtl 
	 end subroutine InsertPrtls
	 
     subroutine ReorderPrtl
         if(modulo(t,prtl_reorder_period).eq.0) call ReorderPrtlArr
     end subroutine ReorderPrtl
	 
     subroutine ReorderPrtlArr
		  implicit none 
          integer :: n,i,j,k,offset,pcount_this
          integer, dimension (mx,my,mz) :: pcount 
		  integer :: count, np_cpu
		  integer :: ind
          pcount=0
		  count=0
		  
      
		  do n=1,used_prtl_arr_size
               if(flvp(n).ne.0) then
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
               if(flvp(n).ne.0) then
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
end subroutine ReorderPrtlArr


subroutine InitPrtlArr(size)
	integer :: size
	
    prtl_arr_size = size
	
	if(allocated(qp)) deallocate(qp,xp,yp,zp,up,vp,wp,tagp,flvp,var1p)
	
	allocate(qp(prtl_arr_size))
    allocate(xp(prtl_arr_size)) 
    allocate(yp(prtl_arr_size)) 
    allocate(zp(prtl_arr_size)) 
    allocate(up(prtl_arr_size)) 
    allocate(vp(prtl_arr_size)) 
    allocate(wp(prtl_arr_size)) 
    allocate(tagp(prtl_arr_size)) 
    allocate(flvp(prtl_arr_size))
    allocate(var1p(prtl_arr_size)) 
  
    qp=0
    used_prtl_arr_size=0
    np=0 
end subroutine InitPrtlArr 


end module mem_prtl