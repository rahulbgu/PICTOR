module savedata_routines
     use parameters
     use vars
     use interpolation, only : InterpEMfield
     implicit none 
     integer :: istart,jstart,kstart,iend,jend,kend
     integer :: fdataxi,fdatayi,fdatazi
     integer :: fdataxf,fdatayf,fdatazf
     real(psn) :: binlen 
     integer :: ldep_edge_send,rdep_edge_send,ldep_edge_recv,rdep_edge_recv
     integer :: ldep_proc_send,rdep_proc_send,ldep_procxind_recv,rdep_procxind_recv     
     
contains
     
     subroutine GetSizeofCollectPrtl
          integer :: i
          tosave_prtl_arr_size=0
          do i=1,used_prtl_arr_size
               if(tagp(i).ne.0) then
                    tosave_prtl_arr_size=tosave_prtl_arr_size+1
               end if
          end do
          do i=1,used_test_prtl_arr_size
               if(tagtp(i).ne.0) then
                    tosave_prtl_arr_size=tosave_prtl_arr_size+1
               end if
          end do
		  
          !print*,'tosave_prtl_arr_size is',tosave_prtl_arr_size
     end subroutine GetSizeofCollectPrtl
     subroutine CollectPrtl(vid)
          integer :: vid,i,j
          
          if(vid.le.10) then 
             j=1
             do i=1,used_prtl_arr_size
                    if(tagp(i).ne.0) then 
                          select case (vid)
                          case(1) 
                            pdata_real(j)=real(xp(i))
                          case(2)
                            pdata_real(j)=real(yp(i))
                          case(3)
                            pdata_real(j)=real(zp(i))
                          case(4)
                            pdata_real(j)=real(up(i))
                          case(5)
                            pdata_real(j)=real(vp(i))
                          case(6)
                            pdata_real(j)=real(wp(i))
                          case(7)
                            pdata_real(j)=real(qp(i))
                          case(8)    
                            pdata_int(j)=tagp(i)
                          case(9)
                            pdata_int(j)=flvp(i)
						  case(10)
                            pdata_real(j)=real(var1p(i))
                          end select              
                          j=j+1
                     end if               
             end do
			 !test particles 
             do i=1,used_test_prtl_arr_size
                    if(tagtp(i).ne.0) then 
                          select case (vid)
                          case(1) 
                            pdata_real(j)=real(xtp(i))
                          case(2)
                            pdata_real(j)=real(ytp(i))
                          case(3)
                            pdata_real(j)=real(ztp(i))
                          case(4)
                            pdata_real(j)=real(utp(i))
                          case(5)
                            pdata_real(j)=real(vtp(i))
                          case(6)
                            pdata_real(j)=real(wtp(i))
                          case(7)
                            pdata_real(j)=real(qtp(i))
                          case(8)    
                            pdata_int(j)=tagtp(i)
                          case(9)
                            pdata_int(j)=flvtp(i)
                          case(10)
                            pdata_real(j)=real(var1tp(i))
                          end select              
                          j=j+1
                     end if               
             end do
			 
			 !test particles 
			 
			 
         end if
          
	      if(vid.ge.11) then 
	           do j=1,tosave_prtl_arr_size
	                      select case (vid)          
	                     case(11)
	                         pdata_real(j)=pdata_local_field(j,1)
	                     case(12)
	                         pdata_real(j)=pdata_local_field(j,2)
	                     case(13)
	                         pdata_real(j)=pdata_local_field(j,3)
	                     case(14)
	                         pdata_real(j)=pdata_local_field(j,4)
	                     case(15)
	                         pdata_real(j)=pdata_local_field(j,5)
	                     case(16)
	                         pdata_real(j)=pdata_local_field(j,6)
	                     case(17)
	                         pdata_real(j)=pdata_local_field(j,1) ! now it contains value of the local current 
	                     case(18)
	                         pdata_real(j)=pdata_local_field(j,2)
	                     case(19)
	                         pdata_real(j)=pdata_local_field(j,3)
	                     end select 
	           end do 
	      end if
          
     end subroutine CollectPrtl
     
     subroutine CalcPrtlLocalEMField 
          integer :: i,j
          j=1
          do i=1,used_prtl_arr_size
               if(tagp(i).ne.0) then 
                    call InterpEMField(xp(i),yp(i),zp(i),pdata_local_field(j,1),pdata_local_field(j,2),&
                    pdata_local_field(j,3),pdata_local_field(j,4),pdata_local_field(j,5),pdata_local_field(j,6),Ex,Ey,Ez,Bx,By,Bz)             
                    j=j+1
               end if
          end do 
          do i=1,used_test_prtl_arr_size
               if(tagtp(i).ne.0) then 				   
                    call InterpEMField(xtp(i),ytp(i),ztp(i),pdata_local_field(j,1),pdata_local_field(j,2),&
                    pdata_local_field(j,3),pdata_local_field(j,4),pdata_local_field(j,5),pdata_local_field(j,6),Ex,Ey,Ez,Bx,By,Bz)             
                    j=j+1
               end if
          end do 
     end subroutine CalcPrtlLocalEMField
     subroutine CalcPrtlLocalCurr 
          integer :: i,j
          j=1
          do i=1,used_prtl_arr_size
               if(tagp(i).ne.0) then 				   
                    call InterpEMField(xp(i),yp(i),zp(i),pdata_local_field(j,1),pdata_local_field(j,2),&
                    pdata_local_field(j,3),pdata_local_field(j,4),pdata_local_field(j,5),pdata_local_field(j,6),Jx,Jy,Jz,Bx,By,Bz)
                    j=j+1
               end if
          end do 
          do i=1,used_test_prtl_arr_size
               if(tagtp(i).ne.0) then 				   
                    call InterpEMField(xtp(i),ytp(i),ztp(i),pdata_local_field(j,1),pdata_local_field(j,2),&
                    pdata_local_field(j,3),pdata_local_field(j,4),pdata_local_field(j,5),pdata_local_field(j,6),Jx,Jy,Jz,Bx,By,Bz)
                    j=j+1
               end if
          end do 
     end subroutine CalcPrtlLocalCurr
     

     subroutine GetSizeofCollectFld
          
          if( (fdataxi.lt.xborders(procxind(proc)+1)) .and. (fdataxf.ge.xborders(procxind(proc))) ) then 
               if(fdataxi.ge.xborders(procxind(proc))) then 
                    istart=fdataxi-xborders(procxind(proc))+3
               else 
                    istart=ceil2(xborders(procxind(proc))-fdataxi,fsave_ratio)-(xborders(procxind(proc))-fdataxi)+3
               end if
               iend=fsave_ratio*int((fdataxf-fdataxi)/fsave_ratio)+fdataxi-xborders(procxind(proc))+3
               if(procxind(proc).eq.nSubDomainsX-1) then !this is to include the right edge on right most proc
                     iend=min(mx-2,iend) 
               else 
                     iend=min(mx-3,iend)
               end if
               fdatax=max(0,floor2(iend-istart,fsave_ratio)/fsave_ratio+1)
          else 
               fdatax=0
          end if
          
          if( (fdatayi.lt.yborders(procyind(proc)+1)) .and. (fdatayf.ge.yborders(procyind(proc))) ) then 
               if(fdatayi.ge.yborders(procyind(proc))) then 
                    jstart=fdatayi-yborders(procyind(proc))+3
               else 
                    jstart=ceil2(yborders(procyind(proc))-fdatayi,fsave_ratio)-(yborders(procyind(proc))-fdatayi)+3
               end if
               jend=fsave_ratio*int((fdatayf-fdatayi)/fsave_ratio)+fdatayi-yborders(procyind(proc))+3
               if(procyind(proc).eq.nSubDomainsY-1) then !this is to include the right edge on top most proc
                     jend=min(my-2,jend) 
               else 
                     jend=min(my-3,jend)
               end if
               fdatay=max(0,floor2(jend-jstart,fsave_ratio)/fsave_ratio+1)               
          else 
               fdatay=0
          end if          
   
#ifndef twoD   
          if( (fdatazi.lt.zborders(proczind(proc)+1)) .and. (fdatazf.ge.zborders(proczind(proc))) ) then 
               if(fdatazi.ge.zborders(proczind(proc))) then 
                    kstart=fdatazi-zborders(proczind(proc))+3
               else 
                    kstart=ceil2(zborders(proczind(proc))-fdatazi,fsave_ratio)-(zborders(proczind(proc))-fdatazi)+3
               end if
               kend=fsave_ratio*int((fdatazf-fdatazi)/fsave_ratio)+fdatazi-zborders(proczind(proc))+3
               if(proczind(proc).eq.nSubDomainsZ-1) then !this is to include the right edge on top most proc
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
     end subroutine CollectFld
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
     integer function floor2real(n1,n2)
          real    :: n1
          integer :: n2
          if(n1.ge.0) then 
               floor2real=n2*int(n1/n2)
          else 
               floor2real=n2*int(n1/n2)-n2
          end if
     end function floor2real



     subroutine CalcEnergy(energy_this) !Calcualtes the total energy in EM field and particles
          real(psn), dimension(4) :: energy_this
          integer :: n,i,j,k
          real(psn) :: elc_energy,ion_energy,gamma,Bfld_energy,Efld_energy 
             elc_energy=0
             ion_energy=0
             Bfld_energy=0
             Efld_energy=0
          do n=1,used_prtl_arr_size
               gamma=sqrt(1+up(n)**2+vp(n)**2+wp(n)**2)     
#ifdef mulflvr
               if(flvp(n).eq.1)  ion_energy=ion_energy+(gamma-1)*abs(qp(n))
               if(flvp(n).eq.2)  elc_energy=elc_energy+(gamma-1)*abs(qp(n))
#else                     
               if(qp(n).gt.0)  ion_energy=ion_energy+(gamma-1)*abs(qp(n))
               if(qp(n).lt.0)  elc_energy=elc_energy+(gamma-1)*abs(qp(n))
#endif
          end do
          ion_energy=ion_energy*massi*c*c
          elc_energy=elc_energy*masse*c*c
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
          Bfld_energy=Bfld_energy*0.5
          Efld_energy=Efld_energy*0.5
          energy_this(1)=elc_energy
          energy_this(2)=ion_energy
          energy_this(3)=Bfld_energy
          energy_this(4)=Efld_energy
          !print*,'Total',prtl_energy+Bfld_energy+Efld_energy
     end subroutine CalcEnergy
     subroutine CalcDivE ! currently written to calucate divergence of electric field in 2D 
          integer :: i,j,k
          Jx=0 ! Jx must be available to store the divergene data 
#ifdef twoD
          do k=1,1
#else
          do k=2,mz-2
#endif
               do j=2,my-2
                    do i=2,mx-2
#ifdef twoD  
                        Jx(i,j,k)=Ex(i,j,k)-Ex(i-1,j,k)+Ey(i,j,k)-Ey(i,j-1,k)
#else
                        Jx(i,j,k)=Ex(i,j,k)-Ex(i-1,j,k)+Ey(i,j,k)-Ey(i,j-1,k)+Ez(i,j,k)-Ez(i,j,k-1)
#endif 
                    end do  
               end do 
          end do 
     end subroutine CalcDivE
     
     
!-----------------------------------------------------------------------------------------------------------------------
!     The following subrotuines computes donwsampled spatial average of various field quantities 
!-----------------------------------------------------------------------------------------------------------------------      


     subroutine CalcPrtlDensity(ch)
          integer :: n,i,j,k,ip,jp,kp,ch
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          Jx=0.0_psn
          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call  DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*qp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*qp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*qp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*qp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*qp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*qp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*qp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*qp(n)
#endif
               end if
          end do
        call NormaliseFldDensity1

     end subroutine CalcPrtlDensity
	 
     subroutine CalcTestPrtlDensity(ch)
          integer :: n,i,j,k,ip,jp,kp,ch
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          Jx=0.0_psn
          do n=1,used_test_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call  DownsampleGridIndex(xtp(n),ytp(n),ztp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*qtp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*qtp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*qtp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*qtp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*qtp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*qtp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*qtp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*qtp(n)
#endif
               end if
          end do
        call NormaliseFldDensity1
     end subroutine CalcTestPrtlDensity
     
     subroutine CalcPrtlChargeFlux(ch)
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg

          Jx=0.0_psn
          Jy=0.0_psn
          Jz=0.0_psn

          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
               call  DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1/sqrt(1+up(n)**2+vp(n)**2+wp(n)**2)
                    vx=up(n)*invg
                    vy=vp(n)*invg
                    vz=wp(n)*invg


                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qp(n)
#endif


               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcPrtlChargeFlux
	 
     subroutine CalcTestPrtlChargeFlux(ch)
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg

          Jx=0.0_psn
          Jy=0.0_psn
          Jz=0.0_psn

          do n=1,used_test_prtl_arr_size
               if(flvtp(n).eq.ch) then
               call  DownsampleGridIndex(xtp(n),ytp(n),ztp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1/sqrt(1+utp(n)**2+vtp(n)**2+wtp(n)**2)
                    vx=utp(n)*invg
                    vy=vtp(n)*invg
                    vz=wtp(n)*invg


                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qtp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qtp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qtp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qtp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qtp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qtp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qtp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qtp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qtp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qtp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qtp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qtp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qtp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qtp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qtp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qtp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qtp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qtp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qtp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qtp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qtp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qtp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qtp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qtp(n)
#endif
               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcTestPrtlChargeFlux
     
     subroutine CalcPrtlMeanVelSquare(ch)
          ! this subroutine calcuates average of cacluate average of velocity square 
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg

          Jx=0
          Jy=0
          Jz=0


          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call  DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1/sqrt(1+up(n)**2+vp(n)**2+wp(n)**2)
                    
                    vx=(up(n)*invg)**2
                    vy=(vp(n)*invg)**2
                    vz=(wp(n)*invg)**2


                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qp(n)
#endif

               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcPrtlMeanVelSquare
	 
     subroutine CalcTestPrtlMeanVelSquare(ch)
          ! this subroutine calcuates average of cacluate average of velocity square 
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg

          Jx=0
          Jy=0
          Jz=0


          do n=1,used_test_prtl_arr_size
               if(flvtp(n).eq.ch) then
                    call  DownsampleGridIndex(xtp(n),ytp(n),ztp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1/sqrt(1+utp(n)**2+vtp(n)**2+wtp(n)**2)
                    
                    vx=(utp(n)*invg)**2
                    vy=(vtp(n)*invg)**2
                    vz=(wtp(n)*invg)**2

                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qtp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qtp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qtp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qtp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qtp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qtp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qtp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qtp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qtp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qtp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qtp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qtp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qtp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qtp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qtp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qtp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qtp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qtp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qtp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qtp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qtp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qtp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qtp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qtp(n)
#endif

               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcTestPrtlMeanVelSquare
     
     subroutine CalcPrtlEdotV(ch)
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg,pEx,pEy,pEz
          Jx=0
          Jy=0
          Jz=0

          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call  DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1.0_psn/sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)
                    call InterpEMfield(xp(n),yp(n),zp(n),pEx,pEy,pEz,vx,vy,vz,Ex,Ey,Ez,Bx,By,Bz)
                    vx=up(n)*invg*pEx
                    vy=vp(n)*invg*pEy
                    vz=wp(n)*invg*pEz

                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qp(n)
#endif                   
               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcPrtlEdotV
	 
     subroutine CalcTestPrtlEdotV(ch)
          integer :: n,i,j,k,ch,ip,jp,kp
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn) ::vx,vy,vz,invg,pEx,pEy,pEz
          Jx=0
          Jy=0
          Jz=0

          do n=1,used_test_prtl_arr_size
               if(flvtp(n).eq.ch) then
                    call  DownsampleGridIndex(xtp(n),ytp(n),ztp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
                    invg=1.0_psn/sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)
                    call InterpEMfield(xtp(n),ytp(n),ztp(n),pEx,pEy,pEz,vx,vy,vz,Ex,Ey,Ez,Bx,By,Bz)
                    vx=utp(n)*invg*pEx
                    vy=vtp(n)*invg*pEy
                    vz=wtp(n)*invg*pEz

                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*vx*qtp(n)
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*vx*qtp(n)
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*vx*qtp(n)
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*vx*qtp(n)
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*vx*qtp(n)
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*vx*qtp(n)
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*vx*qtp(n)
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vx*qtp(n)
#endif

                    Jy(i  ,j  ,k  )=Jy(i  ,j  ,k  )+ Wx *Wy *Wz*vy*qtp(n)
                    Jy(ip ,j  ,k  )=Jy(ip ,j  ,k  )+ Wxp*Wy *Wz*vy*qtp(n)
                    Jy(i  ,jp ,k  )=Jy(i  ,jp ,k  )+ Wx *Wyp*Wz*vy*qtp(n)
                    Jy(ip ,jp ,k  )=Jy(ip ,jp ,k  )+ Wxp*Wyp*Wz*vy*qtp(n)
#ifndef twoD
                    Jy(i  ,j  ,kp )=Jy(i  ,j  ,kp )+ Wx *Wy *Wzp*vy*qtp(n)
                    Jy(ip ,j  ,kp )=Jy(ip ,j  ,kp )+ Wxp*Wy *Wzp*vy*qtp(n)
                    Jy(i  ,jp ,kp )=Jy(i  ,jp ,kp )+ Wx *Wyp*Wzp*vy*qtp(n)
                    Jy(ip ,jp ,kp )=Jy(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vy*qtp(n)
#endif
                    Jz(i  ,j  ,k  )=Jz(i  ,j  ,k  )+ Wx *Wy *Wz*vz*qtp(n)
                    Jz(ip ,j  ,k  )=Jz(ip ,j  ,k  )+ Wxp*Wy *Wz*vz*qtp(n)
                    Jz(i  ,jp ,k  )=Jz(i  ,jp ,k  )+ Wx *Wyp*Wz*vz*qtp(n)
                    Jz(ip ,jp ,k  )=Jz(ip ,jp ,k  )+ Wxp*Wyp*Wz*vz*qtp(n)
#ifndef twoD
                    Jz(i  ,j  ,kp )=Jz(i  ,j  ,kp )+ Wx *Wy *Wzp*vz*qtp(n)
                    Jz(ip ,j  ,kp )=Jz(ip ,j  ,kp )+ Wxp*Wy *Wzp*vz*qtp(n)
                    Jz(i  ,jp ,kp )=Jz(i  ,jp ,kp )+ Wx *Wyp*Wzp*vz*qtp(n)
                    Jz(ip ,jp ,kp )=Jz(ip ,jp ,kp )+ Wxp*Wyp*Wzp*vz*qtp(n)
#endif                   
               end if
          end do
          call NormaliseFldDensity3
     end subroutine CalcTestPrtlEdotV
     
     subroutine CalcPrtlEnergySpatial(ch)
        ! this can possibly be optimised using Jx Jy Jz together as in the previous commented out subroutine 
          integer :: n,i,j,k,ip,jp,kp,ch
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp,eng
          Jx=0
                           
          do n=1,used_prtl_arr_size
               if(flvp(n).eq.ch) then
                    call  DownsampleGridIndex(xp(n),yp(n),zp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)                         
                    eng=sqrt(1.0_psn+up(n)**2+vp(n)**2+wp(n)**2)-1.0_psn
                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*eng
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*eng
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*eng
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*eng
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*eng
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*eng
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*eng
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*eng
#endif
               end if
          end do
        call NormaliseFldDensity1
     end subroutine CalcPrtlEnergySpatial
	 
     subroutine CalcTestPrtlEnergySpatial(ch)
          integer :: n,i,j,k,ip,jp,kp,ch
          real(psn) ::Wx,Wy,Wz,Wxp,Wyp,Wzp,eng
          Jx=0
                           
          do n=1,used_test_prtl_arr_size
               if(flvtp(n).eq.ch) then
                    call  DownsampleGridIndex(xtp(n),ytp(n),ztp(n),i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)                         
                    eng=sqrt(1.0_psn+utp(n)**2+vtp(n)**2+wtp(n)**2)-1.0_psn
                    Jx(i  ,j  ,k  )=Jx(i  ,j  ,k  )+ Wx *Wy *Wz*eng
                    Jx(ip ,j  ,k  )=Jx(ip ,j  ,k  )+ Wxp*Wy *Wz*eng
                    Jx(i  ,jp ,k  )=Jx(i  ,jp ,k  )+ Wx *Wyp*Wz*eng
                    Jx(ip ,jp ,k  )=Jx(ip ,jp ,k  )+ Wxp*Wyp*Wz*eng
#ifndef twoD
                    Jx(i  ,j  ,kp )=Jx(i  ,j  ,kp )+ Wx *Wy *Wzp*eng
                    Jx(ip ,j  ,kp )=Jx(ip ,j  ,kp )+ Wxp*Wy *Wzp*eng
                    Jx(i  ,jp ,kp )=Jx(i  ,jp ,kp )+ Wx *Wyp*Wzp*eng
                    Jx(ip ,jp ,kp )=Jx(ip ,jp ,kp )+ Wxp*Wyp*Wzp*eng
#endif
               end if
          end do
        call NormaliseFldDensity1
     end subroutine CalcTestPrtlEnergySpatial

     subroutine DownsampleGridIndex(x,y,z,i,j,k,ip,jp,kp,Wx,Wy,Wz,Wxp,Wyp,Wzp)
          real(xpsn), intent(in) :: x
          real(ypsn), intent(in) :: y
          real(zpsn), intent(in) :: z 
          integer  , intent(out):: i,j,k,ip,jp,kp
          real(psn), intent(out):: Wx,Wy,Wz,Wxp,Wyp,Wzp
          real(psn)             :: fp
          
               !i=fsave_ratio*int((x-3)/fsave_ratio)+3
               !j=fsave_ratio*int((y-3)/fsave_ratio)+3
               !i=fsave_ratio*int((x-3+xborders(procxind(proc))-fdataxi)/fsave_ratio)+fdataxi-xborders(procxind(proc))+3
               !j=fsave_ratio*int((y-3+yborders(procyind(proc))-fdatayi)/fsave_ratio)+fdatayi-yborders(procyind(proc))+3
               i=floor2real(real(x-3+xborders(procxind(proc))-fdataxi),fsave_ratio)+fdataxi-xborders(procxind(proc))+3
               j=floor2real(real(y-3+yborders(procyind(proc))-fdatayi),fsave_ratio)+fdatayi-yborders(procyind(proc))+3
               
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
               k=floor2real(real(z-3+zborders(proczind(proc))-fdatazi),fsave_ratio)+fdatazi-zborders(proczind(proc))+3
               kp=k+fsave_ratio
               fp =z-k-binlen+0.5
               Wzp=max(min(fp,1.0),0.0)
               Wz =max(min(1-fp,1.0),0.0)     
               k=max(1,k)
               kp=min(mz,kp)
#else
               Wz=1
               k=1
#endif
     end subroutine DownsampleGridIndex
     
!subroutines to nomalised the qunaitties to per unit cell      
     subroutine NormaliseFldDensity3 
#ifdef twoD          
               Jx=Jx/(fsave_ratio**2)
               Jy=Jy/(fsave_ratio**2)
               Jz=Jz/(fsave_ratio**2)
#else
               Jx=Jx/(fsave_ratio**3)
               Jy=Jy/(fsave_ratio**3)
               Jz=Jz/(fsave_ratio**3)
#endif
     end subroutine NormaliseFldDensity3
     subroutine NormaliseFldDensity2
#ifdef twoD
               Jx=Jx/(fsave_ratio**2)
               Jy=Jy/(fsave_ratio**2)
#else
               Jx=Jx/(fsave_ratio**3)
               Jy=Jy/(fsave_ratio**3)
#endif
    end subroutine NormaliseFldDensity2 
     subroutine NormaliseFldDensity1
#ifdef twoD
               Jx=Jx/(fsave_ratio**2)
#else
             Jx=Jx/(fsave_ratio**3)
#endif
    end subroutine NormaliseFldDensity1


!-----------------------------------End of downsampled field quantities subroutines ------------- 
     
     

end module savedata_routines