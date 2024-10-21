module cyl_savedata_routines
     use parameters
     use vars
	 use cyl_common
	 implicit none
 contains 
	 
	 subroutine CalcFldEnergy_cyl(Bfld_energy,Efld_energy)
		 real(dbpsn) :: Bfld_energy,Efld_energy
		 integer :: i1,i,j,k
		 real(psn) :: r,rmin,rp

		 i1=3
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) i1=4     
         rmin=xborders(procxind) + rshift
		 
#ifndef twoD 
         do k=3,mz-3
#else
         do k=1,1
#endif
	         do j=3,my-3
	            do i=i1,mx-3
					 r=(i-3.0_psn)+rmin
					 rp=r+0.5_psn
	                 Bfld_energy=Bfld_energy+( r*Bx(i,j,k)**2  +rp*By(i,j,k)**2 +rp*Bz(i,j,k)**2 )*dtheta
	                 Efld_energy=Efld_energy+( rp*Ex(i,j,k)**2 + r*Ey(i,j,k)**2  +r*Ez(i,j,k)**2 )*dtheta
	            end do
	         end do
         end do
	 end subroutine CalcFldEnergy_cyl
	 
	 subroutine CollectFld_cyl(Fin,fvid,kstart,kend)
         integer, intent(in) :: fvid, kstart, kend
         real(psn),dimension(mx,my,mz), intent(in):: Fin !field to be saved
		 integer :: j,k,j1,k1
		 !fdata(1,:,:) is the value at r=0 which is not 3 but 3.5 in local cord.
         k1=1
#ifndef twoD 
         do k=kstart,kend,fsave_ratio
#else
         do k=1,1
#endif
               j1=1
               do j=3,my-2,fsave_ratio
					 select case (fvid)
				     case(1)
					    fdata(1,j1,k1)=Fin(3,j,k) !Ex
					 case(3)
					    fdata(1,j1,k1)=0.5*(Fin(3,j,k) + Fin(4,j,k)) !Ey
					 case(5)
					    fdata(1,j1,k1)=0.5*(Fin(3,j,k) + Fin(4,j,k)) !Ez
					    !call MeanAxisFldOutput(Fin,kstart,kend) 
					 case(7)
					    fdata(1,j1,k1)=0.5*(Fin(3,j,k) + Fin(4,j,k)) !Bx
					 case(9)
					    fdata(1,j1,k1)=Fin(3,j,k) !By
					 end select
					 j1=j1+1
			   end do
			   k1=k1+1
		  end do  
	 end subroutine CollectFld_cyl
	 
	 subroutine MeanAxisFldOutput(Fld,kstart,kend)
		 real(psn),dimension(mx,my,mz), intent(in):: Fld
		 integer :: k1,kstart,kend
		 real :: sum
		 integer :: i,j,k
		 k1=1 
		 do k=kstart,kend,fsave_ratio
			 sum=0 
			 do j=3,my-3
				 sum=sum+Fld(4,j,k) 
			 end do 
			 sum=sum/ny
			 fdata(1,:,k1)=sum
			 k1=k1+1
		 end do 
	 end subroutine MeanAxisFldOutput
	 
!subroutines to nomalised the qunaitties to per unit cell      
     subroutine NormaliseFldDensity3(Fx,Fy,Fz)
		 real(psn), dimension(mx,my,mz) :: Fx,Fy,Fz 
		 call NormaliseFldDensity_cyl(Fx)
		 call NormaliseFldDensity_cyl(Fy)
		 call NormaliseFldDensity_cyl(Fz)
     end subroutine NormaliseFldDensity3

     subroutine NormaliseFldDensity1(Fx)
		 real(psn), dimension(mx,my,mz) :: Fx
		 call NormaliseFldDensity_cyl(Fx)
     end subroutine NormaliseFldDensity1
	 
	 subroutine FoldInDensityAxis(Fld)
		 real(psn), dimension(mx,my,mz) :: Fld
		 if(inc_axis) then 
		     if(procxind.eq.0) then 
				  Fld(4,:,:)=Fld(4,:,:)-Fld(3,:,:)
				  Fld(3,:,:)=0.0_psn  
			 end if
		 end if
	 end subroutine FoldInDensityAxis
	 
	 subroutine NormaliseFldDensity_cyl(Fld)
		 real(psn), dimension(mx,my,mz) :: Fld
		 integer :: i1,i,j,k
		 real(psn) :: r,rmin, vol
		 
		 i1=1
		 if((inc_axis.eqv..true.).and.(procxind.eq.0)) i1=4
		 rmin=xborders(procxind) + rshift
         do k=1,mz
            do j=1,my
               do i=i1,mx
				   r=(i-3.0_psn)+rmin
#ifdef twoD				   
				   vol=r*dtheta*(fsave_ratio**2)
#else 
                   vol=r*dtheta*(fsave_ratio**3)
#endif				   
				   Fld(i,j,k)=Fld(i,j,k)/vol
			   end do 
		    end do
		  end do  
		  if(inc_axis) call NormaliseFldDensityAxis(Fld)
	 end subroutine NormaliseFldDensity_cyl
	 
	 subroutine NormaliseFldDensityAxis(Fld)
		 real(psn), dimension(mx,my,mz) :: Fld
		 integer :: j,k
		 real(psn) :: sum
		 if(procxind.ne.0) return
		 

#ifndef twoD 
          do k=3,mz-2
#else
          do k=1,1
#endif
             sum=0
			 do j=3,my-2
				 sum=sum+Fld(3,j,k)
			 end do 
	  		   if(fsave_ratio.gt.1) then
#ifdef twoD				    
	  		       sum=sum/(pi*((fsave_ratio/2.0_psn)-0.5_psn)**2) 
#else
                   sum=sum/(pi*fsave_ratio*((fsave_ratio/2.0_psn)-0.5_psn)**2)
#endif				   
	  		   else 
	  			   sum=0
	  		   end if 
			   Fld(3,3,k)=sum
		 end do
		 
		 do k=1,mz
			 Fld(3,1,k)=0
			 Fld(3,2,k)=0
			 do j=3,my-3
				 Fld(3,j,k)=Fld(3,3,k) 
			 end do
			 do j=my-2,my
				 Fld(3,j,k)=0
			 end do  
		 end do 
		 
		   
	 end subroutine NormaliseFldDensityAxis
	 	 
	 subroutine AxisPrtlWt(x,Wx,Wxp)
		 real(psn) :: x
		 real(psn) :: Wx,Wxp
		   if(procxind.eq.0) then
			   if(x.lt.4.0) then
				   if(fsave_ratio.eq.1) then
					   Wxp=Wxp-Wx
					   Wx=0.0
				   end if
				   if(fsave_ratio.gt.1) then
					   Wx=Wx-2*(4.0-x)
				   end if
			   end if 
		   end if 
	 end subroutine AxisPrtlWt

#ifdef gpu	 
	 attributes(device) subroutine AxisPrtlWtGPU(x,Wx,Wxp,procx,fsave_ratio)
		 integer :: procx, fsave_ratio
		 real :: x
		 real :: Wx,Wxp
		   if(procx.eq.0) then
			   if(x.lt.4.0) then
				   if(fsave_ratio.eq.1) then
					   Wxp=Wxp-Wx
					   Wx=0.0
				   end if
				   if(fsave_ratio.gt.1) then
					   Wx=Wx-2*(4.0-x)
				   end if
			   end if 
		   end if 
	 end subroutine AxisPrtlWtGPU
#endif	 


end module cyl_savedata_routines 