module savedata_routines
     use parameters
     use vars
#ifdef cyl
     use cyl_savedata_routines  
#endif 	 
#ifdef gpu
     use savedata_gpu
#endif
     use interpolation
     implicit none 
     integer :: istart,jstart,kstart,iend,jend,kend  
	 real(psn), dimension(6) :: psave_lims_local  
     
contains
     
     subroutine GetSizeofSavePrtl
          integer :: i
          
		  nprtl_save_this=0
          
		  do i=1,used_prtl_arr_size
               if(flvp(i).ne.0 .and. tagp(i).ne.0 .and. within_lims(xp(i),yp(i),zp(i))) nprtl_save_this=nprtl_save_this+1
          end do
		  
          !print*,'nprtl_save_this is',nprtl_save_this
     end subroutine GetSizeofSavePrtl
	 
	 
	 subroutine SetSaveLimsLocal(lims_local, save_lims)
		 real(dbpsn), dimension(6) :: save_lims
		 real(psn), dimension(6) :: lims_local
		 
		 lims_local(1) = save_lims(1) - xborders(procxind) + 3 
		 lims_local(2) = save_lims(2) - xborders(procxind) + 3 
		 lims_local(3) = save_lims(3) - yborders(procyind) + 3 
		 lims_local(4) = save_lims(4) - yborders(procyind) + 3 
		 lims_local(5) = save_lims(5) - zborders(proczind) + 3 
		 lims_local(6) = save_lims(6) - zborders(proczind) + 3 
#ifdef twoD
         lims_local(5) = save_lims(5)  
         lims_local(6) = save_lims(6)
#endif		 
	 end subroutine SetSaveLimsLocal
	 
	 !-------------------------------------------------------------------------------------
	 ! Extract particle attributes to be saved
	 !-------------------------------------------------------------------------------------
	 subroutine GatherSavePrtl(size,q,x,y,z,u,v,w,var1,tag,flv)
		 integer :: size
		 real, dimension(size) :: x,y,z,u,v,w,q,var1
		 integer, dimension(size) :: tag,flv
		 integer :: i, j 
		 
		 j=1
		 do i = 1, used_prtl_arr_size
			  if(flvp(i).ne.0 .and. tagp(i).ne.0 .and. within_lims(xp(i),yp(i),zp(i))) then 
				  q(j) = real(qp(i))
				  x(j) = real(xp(i)) 
				  y(j) = real(yp(i))
				  z(j) = real(zp(i)) 
				  u(j) = real(up(i))
				  v(j) = real(vp(i))
			      w(j) = real(wp(i))
				  tag(j) = tagp(i)
				  flv(j) = flvp(i)
				  j=j+1
			  end if
		 end do 
	 end subroutine GatherSavePrtl
	 
	 !-------------------------------------------------------------------------------------
	 ! Compute local field (E,B, or J) for particles 
	 !-------------------------------------------------------------------------------------
	 subroutine GatherPrtlLocalField(size,x,y,z,pFx,pFy,pFz,VecFx,VecFy,VecFz)
		 integer :: size
		 real, dimension(size) :: x,y,z
		 real, dimension(size) :: pFx,pFy,pFz
		 real(psn), dimension(:,:) :: VecFx,VecFy,VecFz
		 integer :: n 
		 
		 do n=1,size
			 
			 call Interp_VecFld_GridPoints(x(n),y(n),z(n),pFx(n),pFy(n),pFz(n),vecFx,vecFy,vecFz)
			 
		 end do
		   
	 end subroutine GatherPrtlLocalField
	 
     
     subroutine GetSizeofCollectFld
          
          if( (fdataxi.lt.xborders(procxind+1)) .and. (fdataxf.ge.xborders(procxind)) ) then 
               if(fdataxi.ge.xborders(procxind)) then 
                    istart=fdataxi-xborders(procxind)+3
               else 
                    istart=ceil2(xborders(procxind)-fdataxi,fsave_ratio)-(xborders(procxind)-fdataxi)+3
               end if
               iend=fsave_ratio*int((fdataxf-fdataxi)/fsave_ratio)+fdataxi-xborders(procxind)+3
               if(procxind.eq.nSubDomainsX-1) then !this is to include the right edge on right most proc
                     iend=min(mx-2,iend) 
               else 
                     iend=min(mx-3,iend)
               end if
               fdatax=max(0,floor2(iend-istart,fsave_ratio)/fsave_ratio+1)
          else 
               fdatax=0
          end if
          
          if( (fdatayi.lt.yborders(procyind+1)) .and. (fdatayf.ge.yborders(procyind)) ) then 
               if(fdatayi.ge.yborders(procyind)) then 
                    jstart=fdatayi-yborders(procyind)+3
               else 
                    jstart=ceil2(yborders(procyind)-fdatayi,fsave_ratio)-(yborders(procyind)-fdatayi)+3
               end if
               jend=fsave_ratio*int((fdatayf-fdatayi)/fsave_ratio)+fdatayi-yborders(procyind)+3
               if(procyind.eq.nSubDomainsY-1) then !this is to include the right edge on top most proc
                     jend=min(my-2,jend) 
               else 
                     jend=min(my-3,jend)
               end if		   
               fdatay=max(0,floor2(jend-jstart,fsave_ratio)/fsave_ratio+1)               
          else 
               fdatay=0
          end if          
   
#ifndef twoD   
          if( (fdatazi.lt.zborders(proczind+1)) .and. (fdatazf.ge.zborders(proczind)) ) then 
               if(fdatazi.ge.zborders(proczind)) then 
                    kstart=fdatazi-zborders(proczind)+3
               else 
                    kstart=ceil2(zborders(proczind)-fdatazi,fsave_ratio)-(zborders(proczind)-fdatazi)+3
               end if
               kend=fsave_ratio*int((fdatazf-fdatazi)/fsave_ratio)+fdatazi-zborders(proczind)+3
               if(proczind.eq.nSubDomainsZ-1) then !this is to include the right edge on top most proc
                     kend=min(mz-2,kend) 
               else 
                     kend=min(mz-3,kend)
               end if
               fdataz=max(0,floor2(kend-kstart,fsave_ratio)/fsave_ratio+1)               
          else 
               fdataz=0
          end if
#else 
          kstart=1
          kend=1
          fdataz=1
#endif                    
          allocate(fdata(fdatax,fdatay,fdataz))
     end subroutine GetSizeofCollectFld
	 
     
	 
	 integer function ceil2(n1,n2)
          integer :: n1,n2
          if(modulo(n1,n2).eq.0) then 
               ceil2=n1
          else 
               if(n1.gt.0) then 
                  ceil2=n2*int(n1/n2)+n2
               else 
                  ceil2=n2*int(n1/n2)
               end if
          end if          
     end function ceil2
     integer function floor2(n1,n2)
          integer :: n1,n2
          if(n1.ge.0) then 
               floor2=n2*int(n1/n2)
          else 
               floor2=n2*int(n1/n2)-n2
          end if
     end function floor2
	 
	 logical function within_lims(x,y,z)
	 	  real(psn), intent(in) :: x,y,z
		  within_lims = .false.
		  within_lims = (x.gt.psave_lims_local(1).and.x.lt.psave_lims_local(2)) & 
		                .and. (y.gt.psave_lims_local(3).and.y.lt.psave_lims_local(4)) & 
						.and. (z.gt.psave_lims_local(5).and.z.lt.psave_lims_local(6))
 	 end function within_lims 
	 
	 
	 
     subroutine CollectFld(Fin,fvid)
          integer, intent(in) :: fvid
          real(psn),dimension(mx,my,mz), intent(in):: Fin !field to be saved
          integer :: i,j,k,i1,j1,k1
		  
          if((fdatax.eq.0).or.(fdatay.eq.0).or.(fdataz.eq.0)) return 
               
          k1=1
#ifndef twoD 
          do k=kstart,kend,fsave_ratio
#else
          do k=1,1
#endif
                j1=1
               do j=jstart,jend,fsave_ratio
                    i1=1               
                    do i=istart,iend,fsave_ratio
						                       
                         select case (fvid)
                         case(1)
                           fdata(i1,j1,k1)=real((Fin(i-1,j,k)+Fin(i,j,k))*0.5_psn) !Ex
                         case(3)
                           fdata(i1,j1,k1)=real((Fin(i,j-1,k)+Fin(i,j,k))*0.5_psn) !Ey
#ifndef twoD 
                         case(5)
                           fdata(i1,j1,k1)=real((Fin(i,j,k-1)+Fin(i,j,k))*0.5_psn) !Ez
                         case(7)
                           fdata(i1,j1,k1)=real((Fin(i,j,k)+Fin(i,j-1,k)+Fin(i,j,k-1)+Fin(i,j-1,k-1))*0.25_psn) !Bx
                         case(9)
                           fdata(i1,j1,k1)=real((Fin(i,j,k)+Fin(i-1,j,k)+Fin(i,j,k-1)+Fin(i-1,j,k-1))*0.25_psn) !By
#else
	                     case(5)
	                       fdata(i1,j1,k1)=real(Fin(i,j,1)) !Ez
	                     case(7)
	                       fdata(i1,j1,k1)=real((Fin(i,j,k)+Fin(i,j-1,k))*0.5_psn) !Bx
	                     case(9)
	                       fdata(i1,j1,k1)=real((Fin(i,j,k)+Fin(i-1,j,k))*0.5_psn) !By
#endif
                         case(11)
                           fdata(i1,j1,k1)=real((Fin(i,j,k)+Fin(i-1,j,k)+Fin(i,j-1,k)+Fin(i-1,j-1,k))*0.25_psn) !Bz
                         case(13)
                           fdata(i1,j1,k1)=real(Fin(i,j,k)) !charge density
                         end select 
                         i1=i1+1
                    end do
                    j1=j1+1
               end do          
               k1=k1+1
          end do
		  
#ifdef cyl
		  if(inc_axis .and. procxind.eq.0) call CollectFld_cyl(Fin,fvid,kstart,kend)
#endif		   
     
	 end subroutine CollectFld
	


     subroutine CalcEnergy(energy_this) !Calcualtes the total energy in EM field and particles
          real(dbpsn), dimension(Nflvr+2) :: energy_this
          real(dbpsn) :: Bfld_energy,Efld_energy 
          
		  energy_this = 0; Bfld_energy=0; Efld_energy=0
#ifdef cyl
          call CalcFldEnergy_cyl(Bfld_energy,Efld_energy)
#else 
          call CalcFldEnergy(Bfld_energy,Efld_energy)
#endif  

          energy_this(1)=Bfld_energy*0.5_psn
          energy_this(2)=Efld_energy*0.5_psn

		  call CalcPrtlEnergy(energy_this(3:Nflvr+2))
		 
     end subroutine CalcEnergy
	 
	 subroutine CalcPrtlEnergy(energy)
		 integer :: n
		 real(dbpsn), dimension(1:Nflvr) :: energy 
		 real(psn)   :: gamma
		 real(psn)   :: mic2, mec2
		 mic2 = massi*c*c
		 mec2 = masse*c*c
! #ifdef gpu
!          call CalcPrtlEnergyGPU(elc_energy,ion_energy)
! 		 return
! #endif	 
         do n=1,used_prtl_arr_size
              gamma=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
			  !Note : the following two estimates are not applicable in the case of active "ionization", 
			  !since q can't be used to esitamte mass from q/m, and particle's weighting factor should also be accounted for
		      !which is not represented by "q" alone
              if(qp(n).gt.0)  energy(flvp(n))=energy(flvp(n))+(gamma-1.0_psn)*abs(qp(n))*(qmi/flvrqm(flvp(n)))*mic2
              if(qp(n).lt.0)  energy(flvp(n))=energy(flvp(n))+(gamma-1.0_psn)*abs(qp(n))*(qme/flvrqm(flvp(n)))*mec2
         end do
	 end subroutine CalcPrtlEnergy
	 
	 subroutine CalcFldEnergy(Bfld_energy,Efld_energy)
		 real(dbpsn) :: Bfld_energy,Efld_energy
		 integer :: i,j,k
		 
#ifndef twoD 
         do k=3,mz-3
#else
         do k=1,1
#endif
         do j=3,my-3
            do i=3,mx-3
                 Bfld_energy=Bfld_energy+(Bx(i,j,k)**2+By(i,j,k)**2+Bz(i,j,k)**2)
                 Efld_energy=Efld_energy+(Ex(i,j,k)**2+Ey(i,j,k)**2+Ez(i,j,k)**2)
            end do
         end do
         end do
	 end subroutine CalcFldEnergy
	 
	 
     subroutine CalcDivE(div,Fx,Fy,Fz)
		  real(psn), dimension(mx,my,mz) :: div, Fx,Fy,Fz 
          integer :: i,j,k
          div=0 ! Jx must be available to store the divergene data 
#ifdef twoD
          do k=1,1
#else
          do k=2,mz-2
#endif
               do j=2,my-2
                    do i=2,mx-2
#ifdef twoD  
                        div(i,j,k)=Fx(i,j,k)-Fx(i-1,j,k)+Fy(i,j,k)-Fy(i,j-1,k)
#else
                        div(i,j,k)=Fx(i,j,k)-Fx(i-1,j,k)+Fy(i,j,k)-Fy(i,j-1,k)+Fz(i,j,k)-Fz(i,j,k-1)
#endif 
                    end do  
               end do 
          end do 
     end subroutine CalcDivE
	 

     
     
!-----------------------------------------------------------------------------------------------------------------------
!     The following subrotuines computes donwsampled spatial average of various field quantities 
!-----------------------------------------------------------------------------------------------------------------------      
	 
     subroutine CalcPrtlFlux(ch, qnty_id, Fx, Fy, Fz)
		  integer :: qnty_id 
          real(psn), dimension(mx,my,mz) :: Fx, Fy, Fz
		  integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) :: Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) :: invg, qnty_x, qnty_y, qnty_z
		  real(psn) :: pEx,pEy,pEz,pBx,pBy,pBz

          Fx=0.0_psn; Fy=0.0_psn; Fz=0.0_psn

          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1.0_psn/sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
					select case (qnty_id)
						case(1) ! charge flux
							qnty_x = up(n)*invg*qp(n); qnty_y = vp(n)*invg*qp(n); qnty_z = wp(n)*invg*qp(n); 
						case(2)	! velocity^2
	                        qnty_x=qp(n)*(up(n)*invg)**2; qnty_y=qp(n)*(vp(n)*invg)**2; qnty_z=qp(n)*(wp(n)*invg)**2;
						case(3) ! E.V
						    call InterpEMfield(xp(n),yp(n),zp(n),pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)
	                        qnty_x=qp(n)*up(n)*invg*pEx; qnty_y=qp(n)*vp(n)*invg*pEy; qnty_z=qp(n)*wp(n)*invg*pEz;
							
					end select

					
					call AddToGrid(Fx,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,qnty_x)
                    call AddToGrid(Fy,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,qnty_y)
					call AddToGrid(Fz,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,qnty_z)

               end if
          end do
          call NormaliseFldDensity3(Fx,Fy,Fz)
     end subroutine CalcPrtlFlux
	 
     subroutine CalcPrtlDen(ch, qnty_id, Fx)
		  integer :: qnty_id 
          real(psn), dimension(mx,my,mz) :: Fx
		  integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) :: Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) :: qnty
		  real(psn) :: mc2

          Fx=0.0_psn;

		  mc2 = 0 
		  if(FlvrCharge(ch).gt.0) mc2 = (qmi/flvrqm(ch))*massi*c*c
		  if(FlvrCharge(ch).lt.0) mc2 = (qme/flvrqm(ch))*masse*c*c
		  
          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
					select case (qnty_id)
						case(1) ! charge density
							qnty = qp(n)
						case(2)	! energy density
	                        qnty = (sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)-1.0_psn)*qp(n)*mc2
							
					end select

					call AddToGrid(Fx,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,qnty)

               end if
          end do
          call NormaliseFldDensity1(Fx)
     end subroutine CalcPrtlDen
	 	 
	 
	 subroutine AddToGrid(Fx,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp,qnty)
		 real(psn), dimension(mx,my,mz) :: Fx
		 real(psn) :: Wx,Wy,Wz,Wxp,Wyp,Wzp
		 integer   :: i,j,k,ip,jp,kp
		 real(psn) :: qnty
		 
         Fx(i  ,j  ,k  )= Fx(i  ,j  ,k  )+ Wx *Wy *Wz*qnty
         Fx(ip ,j  ,k  )= Fx(ip ,j  ,k  )+ Wxp*Wy *Wz*qnty
         Fx(i  ,jp ,k  )= Fx(i  ,jp ,k  )+ Wx *Wyp*Wz*qnty
         Fx(ip ,jp ,k  )= Fx(ip ,jp ,k  )+ Wxp*Wyp*Wz*qnty
#ifndef twoD
         Fx(i  ,j  ,kp )= Fx(i  ,j  ,kp )+ Wx *Wy *Wzp*qnty
         Fx(ip ,j  ,kp )= Fx(ip ,j  ,kp )+ Wxp*Wy *Wzp*qnty
         Fx(i  ,jp ,kp )= Fx(i  ,jp ,kp )+ Wx *Wyp*Wzp*qnty
         Fx(ip ,jp ,kp )= Fx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*qnty
#endif
	 end subroutine AddToGrid
	

     subroutine DownsampleGridIndex(x,y,z,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
          real(psn), intent(in) :: x,y,z
          integer  , intent(out):: i,j,k,ip,jp,kp
          real(psn), intent(out):: Wx,Wy,Wz,Wxp,Wyp,Wzp
          real                  :: fp

               i=floor2real(real(x-3+xborders(procxind)-fdataxi),fsave_ratio)+fdataxi-xborders(procxind)+3
               j=floor2real(real(y-3+yborders(procyind)-fdatayi),fsave_ratio)+fdatayi-yborders(procyind)+3

               ip=i+fsave_ratio
               jp=j+fsave_ratio
			   
			   
               fp=x-i-binlen+0.5
               Wxp=max(min(fp,1.0),0.0)
               Wx =max(min(1-fp,1.0),0.0)
               fp =y-j-binlen+0.5
               Wyp=max(min(fp,1.0),0.0)
               Wy =max(min(1-fp,1.0),0.0)

               i=max(1,i)
               j=max(1,j)
               ip=min(mx,ip)
               jp=min(my,jp)

#ifndef twoD
               k=floor2real(real(z-3+zborders(proczind)-fdatazi),fsave_ratio)+fdatazi-zborders(proczind)+3
               kp=k+fsave_ratio
               fp =z-k-binlen+0.5
               Wzp=max(min(fp,1.0),0.0)
               Wz =max(min(1.0-fp,1.0),0.0)
               k=max(1,k)
               kp=min(mz,kp)
#else
               Wz=1
               k=1
			   Wzp =0
			   kp =1 
#endif
#ifdef cyl
               call AxisPrtlWt(x,Wx,Wxp)	 
#endif
     end subroutine DownsampleGridIndex
	 
     integer function floor2real(n1,n2)
          real    :: n1
          integer :: n2
          if(n1.ge.0) then 
               floor2real=n2*int(n1/n2)
          else 
               floor2real=n2*int(n1/n2)-n2
          end if
     end function floor2real

	      

#ifndef cyl	 
!subroutines to nomalised the qunaitties to per unit cell      
     subroutine NormaliseFldDensity3(Fx,Fy,Fz)
		 real(psn), dimension(mx,my,mz) :: Fx,Fy,Fz 
#ifdef twoD          
               Fx=Fx/(fsave_ratio**2)
               Fy=Fy/(fsave_ratio**2)
               Fz=Fz/(fsave_ratio**2)
#else
               Fx=Fx/(fsave_ratio**3)
               Fy=Fy/(fsave_ratio**3)
               Fz=Fz/(fsave_ratio**3)
#endif
     end subroutine NormaliseFldDensity3

     subroutine NormaliseFldDensity1(Fx)
		 real(psn), dimension(mx,my,mz) :: Fx
#ifdef twoD
               Fx=Fx/(fsave_ratio**2)
#else
               Fx=Fx/(fsave_ratio**3)
#endif
    end subroutine NormaliseFldDensity1


#endif	

     
end module savedata_routines