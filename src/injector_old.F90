module injector 
	use parameters
	use vars
	use deposit
	use prob_dist
	use memory
	implicit none 
	procedure(fld_ext), pointer :: ydrift =>null()
	procedure(fld_ext), pointer :: ydrift_upper =>null()
	procedure(fld_ext), pointer :: ydrift_lower =>null()
	procedure(func_ext), pointer :: InjTempIon =>null()
	procedure(func_ext), pointer :: InjTempElc =>null()
	procedure(func_ext), pointer :: InjDen =>null()
	procedure(func_ext), pointer :: InjDenUpper =>null()
	procedure(func_ext), pointer :: InjDenLower =>null()
	procedure(fld_ext),pointer   :: InjMagFld
	procedure(fld_ext),pointer   :: InjElcFldUpper
	
	real(psn)  :: vmax_inj_upper, vmax_inj_lower
	abstract interface
	    function func_ext(x,y,z)
	 	     import :: psn
			 real(psn) :: x,y,z
			 real(psn) :: func_ext
	 	end function
		subroutine fld_ext(x,y,z,fx,fy,fz)
			import :: psn
			real(psn) :: x,y,z,fx,fy,fz
		end subroutine
	end interface 
	
contains 
!------------------------------------------	
! Particles
!------------------------------------------	
	
    subroutine InflowBC_Prtl_Top(yinj)
		real(psn) :: yinj
		!call RefPrtlFluidFrameUpperY(yinj)
		call ClearOldPrtlTop(yinj)
		call InjectNewPrtlY(1,yinj)
	end subroutine InflowBC_Prtl_Top
    subroutine InflowBC_Prtl_Bottom(yinj)
		real(psn) :: yinj
		!call RefPrtlFluidFrameLowerY(yinj)
		call ClearOldPrtlBottom(yinj)
		call InjectNewPrtlY(2,yinj)
	end subroutine InflowBC_Prtl_Bottom
	
	subroutine ClearOldPrtlBottom(yinj)
		real(psn) :: yinj,yinj_local
		integer :: n,count 
		
   	    yinj_local=yinj-yborders(procyind(proc))+3
		if(yinj_local.lt.1) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(yp(n).le.yinj_local) then 
				call DeletePrtl(n)
			    count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlBottom
	subroutine ClearOldPrtlTop(yinj)
		real(psn) :: yinj,yinj_local
		integer :: n,count 
		
   	    yinj_local=yinj-yborders(procyind(proc))+3
		if(yinj_local.gt.my) return
		count=0 
		do n=1,used_prtl_arr_size
			if(qp(n).eq.0) cycle
			if(yp(n).ge.yinj_local) then 
				call DeletePrtl(n)
				count=count+1
			end if  
		end do 
		np=np-count
	end subroutine ClearOldPrtlTop
	
	
	subroutine InjectNewPrtlY(side,yinj)
		integer :: side !1=upper, 2=lower
		real(psn) :: yinj,y1,y2
		real(psn) :: xglobal,yglobal,zglobal
		integer :: nprtl_new,n,tag
		real(psn) :: rnd_acpt,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
		real(psn) :: vx,vy,vz
		real(psn) :: Temp

		 
		 call SetInjDomain(side,yinj,nprtl_new,y1,y2)
		 
		 do n=1,nprtl_new
			 
		    call GetNewPrtlPos(xmin,xmax,y1,y2,zmin,zmax,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
			
			call ydrift(xglobal,yglobal,zglobal,vx,vy,vz)
			
			call random_number(rnd_acpt)
			if(rnd_acpt.le.InjDen(xglobal,yglobal,zglobal)) then	
				
				Temp=InjTempIon(xglobal,yglobal,zglobal)
				call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
				call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
				call GetPrtlTag(tag,1)
				call InsertParticleAt(used_prtl_arr_size+1,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,0.0_psn)
				Temp=InjTempElc(xglobal,yglobal,zglobal)
				call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
				call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
				call GetPrtlTag(tag,2)
				call InsertParticleAt(used_prtl_arr_size+2,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,0.0_psn)
				used_prtl_arr_size=used_prtl_arr_size+2
				np=np+2
				
			 end if 
		 end do 
	end subroutine InjectNewPrtlY
	
	subroutine SetInjDomain(side,yinj,nprtl_new,y1,y2)
		integer :: side !1=upper, 2=lower
		real(psn) :: yinj,y1,y2,dy
		integer   :: nprtl_new
		
		 if(side.eq.1) then 
		     y1=max(yinj,real(yborders(procyind(proc))))
	   	     y2=min(yinj+vmax_inj_upper,real(yborders(procyind(proc)+1)))
			 ydrift=>ydrift_upper
			 InjDen=>InjDenUpper
		 end if
	   	 if(side.eq.2) then 
		     y1=max(yinj-vmax_inj_lower,real(yborders(procyind(proc))))
	   	     y2=min(yinj,real(yborders(procyind(proc)+1)))
			 ydrift=>ydrift_lower
			 InjDen=>InjDenLower
		 end if
	  
 
	   	 y1=y1-yborders(procyind(proc))+3 !switch to the local cordinate 
	   	 y2=y2-yborders(procyind(proc))+3 
		 dy=max(0.0_psn,y2-y1)
		 nprtl_new= epc*dy*(xmax-xmin)*(zmax-zmin)
		 if(nprtl_new.eq.0) return 
		 if(nprtl_new.lt.1000) call GetIntPoissonDist(real(real(nprtl_new),psn),nprtl_new)
		 if(used_prtl_arr_size+2*nprtl_new.gt.prtl_arr_size) call ReshapePrtlArr(prtl_arr_size+3*nprtl_new)
		
	end subroutine SetInjDomain
	
	subroutine GetNewPrtlPos(x1,x2,y1,y2,z1,z2,xlocal,ylocal,zlocal,xglobal,yglobal,zglobal)
		real(psn) :: x1,x2,y1,y2,z1,z2
		real(psn) :: xlocal,ylocal,zlocal,xglobal,yglobal,zglobal
		real(psn) :: r1,r2,r3
	 		call random_number(r1)
	 	    call random_number(r2)
	 		call random_number(r3)
	 		xlocal=x1+r1*(x2-x1)
			ylocal=y1+r2*(y2-y1)
			zlocal=z1+r3*(z2-z1)
			
			xglobal=xlocal+xborders(procxind(proc))-3
			yglobal=ylocal+yborders(procyind(proc))-3
#ifndef twoD			
			zglobal=zlocal+zborders(proczind(proc))-3
#else
            zglobal=zlocal
#endif 	

	end subroutine GetNewPrtlPos
			
	
	subroutine RefPrtlFluidFrameUpperY(yinj) !reflect the backward-moving particles in the fluid frame
		real(psn) :: yinj,y1,vx,vy,vz,x0,y0,z0
		real(psn) :: ox,oy,oz
		integer :: n
		
		if( (yinj-yborders(procyind(proc))+3).gt.my-1) return 
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		oz=zborders(proczind(proc))-3
		
		do n=1,used_prtl_arr_size
			call ydrift_upper(xp(n)+ox,yp(n)+oy,zp(n)+oz,vx,vy,vz)
			y1= yinj+vy*c-yborders(procyind(proc))+3 ! instantaneous location of the drifted boundary layer
			if((qp(n).ne.0).and.(yp(n).gt.y1)) then 
				call ReflectPrtlMovingFrameY(vy,up(n),vp(n),wp(n))
				x0=xp(n)
				y0=yp(n)
				z0=zp(n)
				yp(n)=2.0_psn*y1-y0
				call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
			end if
		end do  
	end subroutine RefPrtlFluidFrameUpperY
	
	subroutine RefPrtlFluidFrameLowerY(yinj) !reflect the backward-moving particles in the fluid frame
		real(psn) :: yinj,y1,vx,vy,vz,x0,y0,z0
		real(psn) :: ox,oy,oz
		integer :: n
		
		if( (yinj-yborders(procyind(proc))+3).lt.2) return 
		
		ox=xborders(procxind(proc))-3
		oy=yborders(procyind(proc))-3
		oz=zborders(proczind(proc))-3
		
		do n=1,used_prtl_arr_size
			call ydrift_lower(xp(n)+ox,yp(n)+oy,zp(n)+oz,vx,vy,vz)
			y1= yinj+vy*c-yborders(procyind(proc))+3 ! instantaneous location of the drifted boundary layer
			if((qp(n).ne.0).and.(yp(n).lt.y1)) then 
				call ReflectPrtlMovingFrameY(vy,up(n),vp(n),wp(n))
				x0=xp(n)
				y0=yp(n)
				z0=zp(n)
				yp(n)=2.0_psn*y1-y0
				call DepositCurrentPIC(x0,y0,z0,xp(n),yp(n),zp(n),qp(n))
			end if
		end do  
	end subroutine RefPrtlFluidFrameLowerY

	
	
	subroutine ReflectPrtlMovingFrameY(Drift,ugamma,vgamma,wgamma)
		real(psn) :: Drift,DriftGamma,DriftBeta
		real(psn) :: ugamma,vgamma,wgamma
		real(psn) :: gamma

		if(abs(Drift).lt.1.0) then
	         DriftBeta=Drift
	         DriftGamma=1.0_psn/sqrt((1.0_psn-Drift)*(1.0_psn+Drift))
	    else
	         DriftGamma=abs(Drift)
	         DriftBeta=sqrt((Drift-1.0_psn)*(Drift+1.0_psn))/Drift
	    end if

		gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
		DriftBeta=-DriftBeta
		vgamma=DriftGamma*vgamma + DriftGamma*DriftBeta*gamma!in the drifting frame

		vgamma=-vgamma
		gamma=sqrt(1.0_psn+ugamma**2+vgamma**2+wgamma**2)
		DriftBeta=-DriftBeta
		vgamma=DriftGamma*vgamma + DriftGamma*DriftBeta*gamma! back to the simualtion frame
	end subroutine ReflectPrtlMovingFrameY
	
	!------------------------------------------------------------
	!subroutine to genrate unique tag for particles, warning:: tagging start from 1 if simulation is restarted
	!may not work in some cases when large number of particles are tagged 
	!------------------------------------------------------------
	subroutine GetPrtlTag(tag,FlvID)
		implicit none
		integer :: tag,FlvID
	
		if(TagCounter(FlvID).gt.psave_ratio) TagCounter(FlvID)=1
		if(TagCounter(FlvID).eq.1) then  
			tag=CurrentTagID(FlvID)
		    if(mod(CurrentTagID(FlvID),NtagProcLen).eq.0) CurrentTagID(FlvID)=CurrentTagID(FlvID)+NtagProcLen*(nproc-1)
	        CurrentTagID(FlvID)=CurrentTagID(FlvID)+1
	    else
			tag=0
		end if	
		TagCounter(FlvID)=TagCounter(FlvID)+1
	end subroutine GetPrtlTag

!------------------------------------------	
! Fields
!------------------------------------------	

	subroutine InflowBC_Fld(yinj,side)
		integer :: side ! 1=upper, 2=lower
		real(psn) :: yinj,yinj_local
		real(psn):: b_x,b_y,b_z,vx,vy,vz,x0,y0,z0
	    integer :: i,j,k, ji_loop1, ji_loop2, jf_loop1, jf_loop2 
		
		yinj_local= yinj-yborders(procyind(proc))+3 
		if((side.eq.1).and.(yinj_local.gt.my)) return
		if((side.eq.2).and.(yinj_local.lt.1)) return
		if(side.eq.1) then
			ji_loop1=ceiling(yinj_local-0.5_psn); jf_loop1=my
			ji_loop2=ceiling(yinj_local); jf_loop2=my
			ydrift=>ydrift_upper
		end if
		if(side.eq.2) then
			ji_loop1=1; jf_loop1=floor(yinj_local-0.5_psn)
			ji_loop2=1; jf_loop2=floor(yinj_local)
			ydrift=>ydrift_lower
		end if
		
		x0=xborders(procxind(proc))-3
		y0=yborders(procyind(proc))-3
		z0=zborders(proczind(proc))-3
#ifdef twoD
        z0=0
#endif		
		 do k=1,mz
	          do j=ji_loop1,jf_loop1
				  do i=1,mx
					  call InjMagFld(real(i,psn)+x0,real(j+0.5_psn,psn)+y0,real(k,psn)+z0,b_x,b_y,b_z)
					  call ydrift(real(i,psn)+x0,real(j+0.5_psn,psn)+y0,real(k,psn)+z0,vx,vy,vz)
					  Ey(i,j,k)= vx * b_z - vz * b_x
					  
					  call InjMagFld(real(i,psn)+x0,real(j+0.5_psn,psn)+y0,real(k+0.5_psn,psn)+z0,b_x,b_y,b_z)
					  Bx(i,j,k)=b_x
					  
					  call InjMagFld(real(i+0.5_psn,psn)+x0,real(j+0.5_psn,psn)+y0,real(k,psn)+z0,b_x,b_y,b_z)
					  Bz(i,j,k)=b_z
				  end do
			  end do 
			  do j=ji_loop2,jf_loop2
				  do i=1,mx
					  
					  call InjMagFld(real(i+0.5_psn,psn)+x0,real(j,psn)+y0,real(k,psn)+z0,b_x,b_y,b_z)
					  call ydrift(real(i+0.5_psn,psn)+x0,real(j,psn)+y0,real(k,psn)+z0,vx,vy,vz)
					  Ex(i,j,k)= vz * b_y - vy * b_z 
					  
					  call InjMagFld(real(i,psn)+x0,real(j,psn)+y0,real(k+0.5+psn,psn)+z0,b_x,b_y,b_z)
					  call ydrift(real(i,psn)+x0,real(j,psn)+y0,real(k+0.5_psn,psn)+z0,vx,vy,vz)
					  Ez(i,j,k)= vy * b_x - vx * b_y
					  
					  call InjMagFld(real(i+0.5_psn,psn)+x0,real(j,psn)+y0,real(k+0.5_psn,psn)+z0,b_x,b_y,b_z)
					  By(i,j,k)=b_y
				  end do 
			  end do  
	     end do
	end subroutine InflowBC_Fld
	
	
! older injection algorithm: particles are reflected in the fluid frame	

! 	subroutine InjectNewPrtlY(side,yinj)
! 		integer :: side !1=upper, 2=lower
! 		real(psn) :: yinj,y1,y2,dy,dely,yinj_local
! 		real(psn) :: xglobal,yglobal,zglobal
! 		integer :: nprtl_new,n,tag
! 		real(psn) :: r1,r2,r3,rnd_acpt,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma
! 		real(psn) :: vx,vy,vz
! 		real(psn) :: Temp
!
! 	   	 !new particles will be injected between y1 and y2, but only in the region emptied out by the older drifting particles
! 	   	 yinj_local=yinj-yborders(procyind(proc))+3
! 		 if(side.eq.1) then
! 		     y1=max(yinj-vmax_inj_upper,real(yborders(procyind(proc))))
! 	   	     y2=min(yinj,real(yborders(procyind(proc)+1)))
! 			 ydrift=>ydrift_upper
! 			 InjDen=>InjDenUpper
! 		 end if
! 	   	 if(side.eq.2) then
! 		     y1=max(yinj,real(yborders(procyind(proc))))
! 	   	     y2=min(yinj+vmax_inj_lower,real(yborders(procyind(proc)+1)))
! 			 ydrift=>ydrift_lower
! 			 InjDen=>InjDenLower
! 		 end if
!
!
! 	   	 y1=y1-yborders(procyind(proc))+3 !switch to the local cordinate
! 	   	 y2=y2-yborders(procyind(proc))+3
! 		 dy=max(0.0_psn,y2-y1)
! 		 nprtl_new= epc*dy*(xmax-xmin)*(zmax-zmin)
! 		 if(nprtl_new.eq.0) return
! 		 if(nprtl_new.lt.1000) call GetIntPoissonDist(real(real(nprtl_new),psn),nprtl_new)
! 		 if(used_prtl_arr_size+2*nprtl_new.gt.prtl_arr_size) call ReshapePrtlArr(prtl_arr_size+3*nprtl_new)
!
! 		 do n=1,nprtl_new
! 	 		call random_number(r1)
! 	 	    call random_number(r2)
! 	 		call random_number(r3)
! 	 		xlocal=xmin+r1*(xmax-xmin)
! 			ylocal=y1+r2*(y2-y1)
! 			zlocal=zmin+r3*(zmax-zmin)
!
! 			xglobal=xlocal+xborders(procxind(proc))-3
! 			yglobal=ylocal+yborders(procyind(proc))-3
! #ifndef twoD
! 			zglobal=zlocal+zborders(proczind(proc))-3
! #else
!             zglobal=zlocal
! #endif
! 			call ydrift(xglobal,yglobal,zglobal,vx,vy,vz)
! 			if(side.eq.1) dely=ylocal-(yinj_local + vy*c)
! 			if(side.eq.2) dely=(yinj_local + vy*c)-ylocal
!
! 			if(dely.ge.0) then
! 				call random_number(rnd_acpt)
! 				if(rnd_acpt.le.InjDen(xglobal,yglobal,zglobal)) then
! 					Temp=InjTempIon(xglobal,yglobal,zglobal)
! 					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
! 					call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
! 					call GetPrtlTag(tag,1)
! 					call InsertParticleAt(used_prtl_arr_size+1,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,1.0_psn,tag,1,0.0_psn)
! 					Temp=InjTempElc(xglobal,yglobal,zglobal)
! 					call GetVelGamma_MaxBolt(Temp,ugamma,vgamma,wgamma)
! 					call AddDriftVel(ugamma,vgamma,wgamma,vx,vy,vz)
! 					call GetPrtlTag(tag,2)
! 					call InsertParticleAt(used_prtl_arr_size+2,xlocal,ylocal,zlocal,ugamma,vgamma,wgamma,-1.0_psn,tag,2,0.0_psn)
! 					used_prtl_arr_size=used_prtl_arr_size+2
! 					np=np+2
! 			    end if
! 			end if
! 		 end do
! 	end subroutine InjectNewPrtlY
	
end module injector