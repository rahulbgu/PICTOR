module interpolation 
     use parameters
     use vars
     implicit none 
contains 

     
subroutine InterpEMfield(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,E_x,E_y,E_z,B_x,B_y,B_z)
     real(psn), intent(IN):: x
     real(psn), intent(IN):: y
     real(psn), intent(IN):: z
     real(psn), intent(INOUT) :: pEx,pEy,pEz,pBx,pBy,pBz
     real(psn), dimension(mx,my,mz), intent(IN) :: E_x,E_y,E_z,B_x,B_y,B_z
     
     
     !call InterpFields_CellCenters_2D(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
     call InterpFields_GridEdges(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,E_x,E_y,E_z,B_x,B_y,B_z)
 

end subroutine InterpEMfield     
     
!-----------------------------------------------------------------------------------------     
! The following scheme interpolates EM flds from grid points
! Note: This shcheme gives 0 self-force
!-----------------------------------------------------------------------------------------     
#ifdef twoD
subroutine InterpFields_GridPoints(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
     real(psn), intent(IN) :: x
     real(psn), intent(IN) :: y
     real(psn), intent(IN) :: z
     real(psn),intent(INOUT):: pEx,pEy,pEz,pBx,pBy,pBz
     real(psn):: dx,dy
     integer :: i,j
     real(psn) :: w1,w2,w3,w4
     
     i=floor(x)
     j=floor(y)
     dx=x-i
     dy=y-j
     w1=(1-dx)*(1-dy)
     w2=dx    *(1-dy)
     w3=(1-dx)*dy    
     w4=dx    *dy    
     
     
     pEx=w1*(Ex(i-1,j  ,1) + Ex(i  ,j  ,1))+&
         w2*(Ex(i  ,j  ,1) + Ex(i+1,j  ,1))+&
          w3*(Ex(i-1,j+1,1) + Ex(i  ,j+1,1))+&
          w4*(Ex(i  ,j+1,1) + Ex(i+1,j+1,1))
          
     pEy=w1*(Ey(i  ,j-1,1) + Ey(i  ,j  ,1))+&
         w2*(Ey(i+1,j-1,1) + Ey(i+1,j  ,1))+&
          w3*(Ey(i  ,j  ,1) + Ey(i  ,j+1,1))+&
          w4*(Ey(i+1,j  ,1) + Ey(i+1,j+1,1))
          
     pEz=w1* Ez(i  ,j  ,1) +&
         w2* Ez(i+1,j  ,1) +&
          w3* Ez(i  ,j+1,1) +&
          w4* Ez(i+1,j+1,1)
          
     pBx=w1*(Bx(i  ,j-1,1) + Bx(i  ,j  ,1))+&
         w2*(Bx(i+1,j-1,1) + Bx(i+1,j  ,1))+&
          w3*(Bx(i  ,j  ,1) + Bx(i  ,j+1,1))+&
          w4*(Bx(i+1,j  ,1) + Bx(i+1,j+1,1))
     
     pBy=w1*(By(i-1,j  ,1) + By(i  ,j  ,1))+&
         w2*(By(i  ,j  ,1) + By(i+1,j  ,1))+&
          w3*(By(i-1,j+1,1) + By(i  ,j+1,1))+&
          w4*(By(i  ,j+1,1) + By(i+1,j+1,1))
          
     pBz=w1*(Bz(i-1,j-1,1) + Bz(i  ,j-1,1) + Bz(i-1,j  ,1) + Bz(i  ,j  ,1))+&
         w2*(Bz(i  ,j-1,1) + Bz(i+1,j-1,1) + Bz(i  ,j  ,1) + Bz(i+1,j  ,1))+&
          w3*(Bz(i-1,j  ,1) + Bz(i  ,j  ,1) + Bz(i-1,j+1,1) + Bz(i  ,j+1,1))+&
          w4*(Bz(i  ,j  ,1) + Bz(i+1,j  ,1) + Bz(i  ,j+1,1) + Bz(i+1,j+1,1))
          
     
                    pEx=pEx*0.5_psn
                    pEy=pEy*0.5_psn
                    pEz=pEz
#ifndef Bext0
                    pBx=pBx*0.5_psn
                    pBy=pBy*0.5_psn
                    pBz=pBz*0.25_psn
#else
                    pBx=pBx*0.5_psn +Bx_ext0
                    pBy=pBy*0.5_psn +By_ext0
                    pBz=pBz*0.25_psn+Bz_ext0
#endif     
     
end subroutine InterpFields_GridPoints
#endif

!-----------------------------------------------------------------------------------------     
! The following method interpolates EM fields directly from the points EM flds are recorded on Yee-Mesh 
!-----------------------------------------------------------------------------------------     
#ifdef twoD
     subroutine InterpFields_GridEdges(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,E_x,E_y,E_z,B_x,B_y,B_z) !interpolation in 2D 
          real(psn), intent(IN) :: x
          real(psn), intent(IN) :: y
          real(psn), intent(IN) :: z
          real(psn), intent(INOUT):: pEx,pEy,pEz,pBx,pBy,pBz
          real(psn), dimension(mx,my,mz), intent(IN) :: E_x,E_y,E_z,B_x,B_y,B_z
          integer :: i,j
          real(psn):: dx,dy

               i=floor(x+0.5_psn)
               j=floor(y)
               dx=x-i+0.5_psn
               dy=y-j
               call Interp2D(i-1,j,E_x,dx,dy,pEx)
               
               i=floor(x)
               j=floor(y+0.5_psn)
               dx=x-i
               dy=y-j+0.5_psn
               call Interp2D(i,j-1,E_y,dx,dy,pEy)
               
               i=floor(x)
               j=floor(y)
               dx=x-i
               dy=y-j
               call Interp2D(i,j,E_z,dx,dy,pEz)

! now interpolate the magnetic field                
               i=floor(x)
               j=floor(y-0.5_psn)
               dx=x-i
               dy=y-j-0.5_psn
               call Interp2D(i,j,B_x,dx,dy,pBx)
               
               i=floor(x-0.5_psn)
               j=floor(y)
               dx=x-i-0.5_psn
               dy=y-j
               call Interp2D(i,j,B_y,dx,dy,pBy)
               
               i=floor(x-0.5_psn)
               j=floor(y-0.5_psn)
               dx=x-i-0.5_psn
               dy=y-j-0.5_psn
               call Interp2D(i,j,B_z,dx,dy,pBz)

                    pBx=pBx+Bx_ext0
                    pBy=pBy+By_ext0
                    pBz=pBz+Bz_ext0
          
     end subroutine InterpFields_GridEdges

#else           
     subroutine InterpFields_GridEdges(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,E_x,E_y,E_z,B_x,B_y,B_z) ! interpolation in 3D
          real(psn), intent(IN) :: x
          real(psn), intent(IN) :: y
          real(psn), intent(IN) :: z
          real(psn), INTENT(INOUT)::pEx,pEy,pEz,pBx,pBy,pBz
          real(psn), dimension(mx,my,mz), intent(IN) :: E_x,E_y,E_z,B_x,B_y,B_z
          integer :: i,j,k
          real(psn):: dx,dy,dz
          
               i=floor(x+0.5_psn)
               j=floor(y)
               k=floor(z)
               dx=x-i+0.5_psn
               dy=y-j
               dz=z-k
               call Interp3D(i-1,j,k,E_x,dx,dy,dz,pEx)
               
               i=floor(x)
               j=floor(y+0.5_psn)
               k=floor(z)
               dx=x-i
               dy=y-j+0.5_psn
               dz=z-k
               call Interp3D(i,j-1,k,E_y,dx,dy,dz,pEy)
               
               i=floor(x)
               j=floor(y)
               k=floor(z+0.5_psn)
               dx=x-i
               dy=y-j
               dz=z-k+0.5_psn
               call Interp3D(i,j,k-1,E_z,dx,dy,dz,pEz)

! now interpolate the magnetic field                
               i=floor(x)
               j=floor(y-0.5_psn)
               k=floor(z-0.5_psn)
               dx=x-i
               dy=y-j-0.5_psn
               dz=z-k-0.5_psn
               call Interp3D(i,j,k,B_x,dx,dy,dz,pBx)
               
               i=floor(x-0.5_psn)
               j=floor(y)
               k=floor(z-0.5_psn)
               dx=x-i-0.5_psn
               dy=y-j
               dz=z-k-0.5_psn
               call Interp3D(i,j,k,B_y,dx,dy,dz,pBy)
               
               i=floor(x-0.5_psn)
               j=floor(y-0.5_psn)
               k=floor(z)
               dx=x-i-0.5_psn
               dy=y-j-0.5_psn
               dz=z-k
               call Interp3D(i,j,k,B_z,dx,dy,dz,pBz)

                    pBx=pBx+Bx_ext0
                    pBy=pBy+By_ext0
                    pBz=pBz+Bz_ext0
          
     end subroutine InterpFields_GridEdges
#endif     
     


!subroutines used in the Grid-point interpolation method      
     subroutine Interp3D(i,j,k,F,dx,dy,dz,res)
          integer :: i,j,k
          real(psn) :: dx,dy,dz,res
          real(psn),dimension(mx,my,mz) :: F
          res=(1.0_psn-dx) *(1.0_psn-dy)*(1.0_psn-dz)*F(i,  j  ,k)+&
               dx          *(1.0_psn-dy)*(1.0_psn-dz)*F(i+1,j  ,k)+&
               (1.0_psn-dx)*dy          *(1.0_psn-dz)*F(i  ,j+1,k)+&
               (1.0_psn-dx)*(1.0_psn-dy)*dz          *F(i,  j,  k+1)+&
               (1.0_psn-dx)*dy          *dz          *F(i,  j+1,k+1)+&
               dx          *(1.0_psn-dy)*dz          *F(i+1,j,  k+1)+&
               dx          *dy          *(1.0_psn-dz)*F(i+1,j+1,k)+&
               dx          *dy          *dz          *F(i+1,j+1,k+1)
     end subroutine Interp3D
     subroutine Interp2D(i,j,F,dx,dy,res)
          integer :: i,j
          real(psn) :: dx,dy,res
          real(psn),dimension(mx,my,1) :: F
          res=(1.0_psn-dx)*(1.0_psn-dy)*F(i,  j  ,1)+&
              dx          *(1.0_psn-dy)*F(i+1,j  ,1)+&
              (1.0_psn-dx)*dy          *F(i  ,j+1,1)+&
               dx         *dy          *F(i+1,j+1,1)
     end subroutine Interp2D  
!-------------------------------------------------------------------------------------------------------
! Interpolation from grid points where the (average) field values are obtained by averaging staggered E,B field components on the Yee-Mesh 
!------------------------------------------------------------------------------------------------------- 
subroutine InterpEMFldGridPointYeeMesh(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz,Ex,Ey,Ez,Bx,By,Bz)
	     real(psn), intent(IN)    :: x,y,z
	     real(psn), intent(INOUT) :: pEx,pEy,pEz,pBx,pBy,pBz
		 real(psn), dimension(mx,my,mz), intent(IN) :: Ex,Ey,Ez,Bx,By,Bz
         integer :: i,j,k
         real(psn):: dx,dy,dz
#ifdef twoD
         real(psn), dimension(4) :: FldVec,wt
#else
         real(psn), dimension(8) :: FldVec,wt 
#endif
		 
		  i=x
		  j=y
		  dx=x-i
		  dy=y-j 
#ifndef twoD 
          k=z
          dz=z-k 		  
	      wt(1)=(1.0_psn-dx)*(1.0_psn-dy)*(1.0_psn-dz)
	      wt(2)=dx          *(1.0_psn-dy)*(1.0_psn-dz)
	      wt(3)=(1._psn-dx) *dy          *(1.0_psn-dz)
	      wt(4)=dx          *dy          *(1.0_psn-dz)
	      wt(5)=(1.0_psn-dx)*(1.0_psn-dy)*dz
	      wt(6)=dx          *(1.0_psn-dy)*dz
	      wt(7)=(1.0_psn-dx)*dy          *dz
	      wt(8)=dx          *dy          *dz
#else
          k=1 
		  wt(1)=(1.0_psn-dx)*(1.0_psn-dy)
		  wt(2)=dx          *(1.0_psn-dy)
		  wt(3)=(1._psn-dx) *dy          
		  wt(4)=dx          *dy          
#endif  
		  FldVec(1)=Ex(i-1,j,k)+Ex(i,j,k)
		  FldVec(2)=Ex(i,j,k)+Ex(i+1,j,k)
		  FldVec(3)=Ex(i-1,j+1,k)+Ex(i,j+1,k)
		  FldVec(4)=Ex(i,j+1,k)+Ex(i+1,j+1,k)
#ifndef twoD 		  
		  FldVec(5)=Ex(i-1,j,k+1)+Ex(i,j,k+1)
		  FldVec(6)=Ex(i,j,k+1)+Ex(i+1,j,k+1)
		  FldVec(7)=Ex(i-1,j+1,k+1)+Ex(i,j+1,k+1)
		  FldVec(8)=Ex(i,j+1,k+1)+Ex(i+1,j+1,k+1)
#endif 
		  pEx=sum(0.5_psn*wt*FldVec) 
		  
		  FldVec(1)=Ey(i,j-1,k)+Ey(i,j,k)
		  FldVec(2)=Ey(i+1,j-1,k)+Ey(i+1,j,k)
		  FldVec(3)=Ey(i,j,k)+Ey(i,j+1,k)
		  FldVec(4)=Ey(i+1,j,k)+Ey(i+1,j+1,k)
#ifndef twoD 		  
		  FldVec(5)=Ey(i,j-1,k+1)+Ey(i,j,k+1)
		  FldVec(6)=Ey(i+1,j-1,k+1)+Ey(i+1,j,k+1)
		  FldVec(7)=Ey(i,j,k+1)+Ey(i,j+1,k+1)
		  FldVec(8)=Ey(i+1,j,k+1)+Ey(i+1,j+1,k+1)
#endif
		  pEy=sum(0.5_psn*wt*FldVec)

#ifndef twoD 		  
		  FldVec(1)=Ez(i,j,k-1)+Ez(i,j,k)
		  FldVec(2)=Ez(i+1,j,k-1)+Ez(i+1,j,k)
		  FldVec(3)=Ez(i,j+1,k-1)+Ez(i,j+1,k)
		  FldVec(4)=Ez(i+1,j+1,k-1)+Ez(i+1,j+1,k)		  
		  FldVec(5)=Ez(i,j,k)+Ez(i,j,k+1)
		  FldVec(6)=Ez(i+1,j,k)+Ez(i+1,j,k+1)
		  FldVec(7)=Ez(i,j+1,k)+Ez(i,j+1,k+1)
		  FldVec(8)=Ez(i+1,j+1,k)+Ez(i+1,j+1,k+1)
          pEz=sum(0.5_psn*wt*FldVec)
#else
		  FldVec(1)=Ez(i,j,k)
		  FldVec(2)=Ez(i+1,j,k)
		  FldVec(3)=Ez(i,j+1,k)
		  FldVec(4)=Ez(i+1,j+1,k)
          pEz=sum(wt*FldVec)
#endif 
 

#ifndef twoD 		  
		  FldVec(1)=Bx(i,j-1,k-1)+Bx(i,j,k-1)+Bx(i,j-1,k)+Bx(i,j,k)
		  FldVec(2)=Bx(i+1,j-1,k-1)+Bx(i+1,j,k-1)+Bx(i+1,j-1,k)+Bx(i+1,j,k)  	  
		  FldVec(3)=Bx(i,j,k-1)+Bx(i,j+1,k-1)+Bx(i,j,k)+Bx(i,j+1,k)	
		  FldVec(4)=Bx(i+1,j,k-1)+Bx(i+1,j+1,k-1)+Bx(i+1,j,k)+Bx(i+1,j+1,k)	
		  FldVec(5)=Bx(i,j-1,k)+Bx(i,j,k)+Bx(i,j-1,k+1)+Bx(i,j,k+1)	
		  FldVec(6)=Bx(i+1,j-1,k)+Bx(i+1,j,k)+Bx(i+1,j-1,k+1)+Bx(i+1,j,k+1)  	  
		  FldVec(7)=Bx(i,j,k)+Bx(i,j+1,k)+Bx(i,j,k+1)+Bx(i,j+1,k+1)	
		  FldVec(8)=Bx(i+1,j,k)+Bx(i+1,j+1,k)+Bx(i+1,j,k+1)+Bx(i+1,j+1,k+1)		
          pBx=sum(0.25_psn*wt*FldVec) +Bx_ext0
#else
		  FldVec(1)=Bx(i,j-1,k)+Bx(i,j,k)	
		  FldVec(2)=Bx(i+1,j-1,k)+Bx(i+1,j,k)  	  
		  FldVec(3)=Bx(i,j,k)+Bx(i,j+1,k)	
		  FldVec(4)=Bx(i+1,j,k)+Bx(i+1,j+1,k)		
		  pBx=sum(0.5_psn*wt*FldVec) +Bx_ext0
#endif 
		  
#ifndef twoD 		  
		  FldVec(1)=By(i-1,j,k-1)+By(i-1,j,k)+By(i,j,k-1)+By(i,j,k) 
		  FldVec(2)=By(i,j,k-1)+By(i,j,k)+By(i+1,j,k-1)+By(i+1,j,k) 
		  FldVec(3)=By(i-1,j+1,k-1)+By(i-1,j+1,k)+By(i,j+1,k-1)+By(i,j+1,k)
		  FldVec(4)=By(i,j+1,k-1)+By(i,j+1,k)+By(i+1,j+1,k-1)+By(i+1,j+1,k)
		  FldVec(5)=By(i-1,j,k)+By(i-1,j,k+1)+By(i,j,k)+By(i,j,k+1)
		  FldVec(6)=By(i,j,k)+By(i,j,k+1)+By(i+1,j,k)+By(i+1,j,k+1)
		  FldVec(7)=By(i-1,j+1,k)+By(i-1,j+1,k+1)+By(i,j+1,k)+By(i,j+1,k+1)
		  FldVec(8)=By(i,j+1,k)+By(i,j+1,k+1)+By(i+1,j+1,k)+By(i+1,j+1,k+1)
          pBy=sum(0.25_psn*wt*FldVec) +By_ext0
#else
		  FldVec(1)=By(i-1,j,k)+By(i,j,k)
		  FldVec(2)=By(i,j,k)+By(i+1,j,k)
		  FldVec(3)=By(i-1,j+1,k)+By(i,j+1,k)
		  FldVec(4)=By(i,j+1,k)+By(i+1,j+1,k)
		  pBy=sum(0.5_psn*wt*FldVec) +By_ext0
#endif 
		  
		  FldVec(1)=Bz(i-1,j-1,k)+Bz(i-1,j,k)+Bz(i,j-1,k)+Bz(i,j,k)
		  FldVec(2)=Bz(i,j-1,k)+Bz(i,j,k)+Bz(i+1,j-1,k)+Bz(i+1,j,k)
		  FldVec(3)=Bz(i-1,j,k)+Bz(i-1,j+1,k)+Bz(i,j,k)+Bz(i,j+1,k)
		  FldVec(4)=Bz(i,j,k)+Bz(i,j+1,k)+Bz(i+1,j,k)+Bz(i+1,j+1,k)
#ifndef twoD 		
		  FldVec(5)=Bz(i-1,j-1,k+1)+Bz(i-1,j,k+1)+Bz(i,j-1,k+1)+Bz(i,j,k+1)
		  FldVec(6)=Bz(i,j-1,k+1)+Bz(i,j,k+1)+Bz(i+1,j-1,k+1)+Bz(i+1,j,k+1)
		  FldVec(7)=Bz(i-1,j,k+1)+Bz(i-1,j+1,k+1)+Bz(i,j,k+1)+Bz(i,j+1,k+1)
		  FldVec(8)=Bz(i,j,k+1)+Bz(i,j+1,k+1)+Bz(i+1,j,k+1)+Bz(i+1,j+1,k+1)
#endif 
		  pBz=sum(0.25*wt*FldVec) +Bz_ext0
	
end subroutine InterpEMFldGridPointYeeMesh
     
     
!-------------------------------------------------------------------------------------------------------
! The following interpolation method assumes that grid qnuatites are evalauated at the center of the cell, 
! the linear interpolation is done from the cell center to a specified location  
!------------------------------------------------------------------------------------------------------- 

subroutine InterpVectorField_CellCenters(x,y,z,pFx,pFy,pFz,F_x,F_y,F_z)
     real(psn):: x,y,z
     real(psn):: pFx,pFy,pFz
     real(psn), dimension(mx,my,mz):: F_x,F_y,F_z 
     real(psn):: dx,dy
     integer :: i,j,k
     real(psn) :: w1,w2,w3,w4
#ifndef twoD
    real(psn) ::w5,w6,w7,w8
     real(psn) ::dz
#endif

#ifdef twoD          
     i=floor(x-0.5_psn)
     j=floor(y-0.5_psn)
     dx=x-i-0.5_psn
     dy=y-j-0.5_psn
     w1=(1.0_psn-dx)*(1.0_psn-dy)
     w2=dx          *(1.0_psn-dy)
     w3=(1.0_psn-dx)*dy    
     w4=dx          *dy    
     
     pFx=w1*F_x(i,j,1) +w2*F_x(i+1,j,1) +w3*F_x(i,j+1,1) +w4*F_x(i+1,j+1,1)
     pFy=w1*F_y(i,j,1) +w2*F_y(i+1,j,1) +w3*F_y(i,j+1,1) +w4*F_y(i+1,j+1,1)
     pFz=w1*F_z(i,j,1) +w2*F_z(i+1,j,1) +w3*F_z(i,j+1,1) +w4*F_z(i+1,j+1,1)

#else      
     i=floor(x-0.5)
     j=floor(y-0.5)
     k=floor(z-0.5)
     dx=x-i-0.5
     dy=y-j-0.5
     dz=z-k-0.5
     w1=(1-dx)*(1-dy)*(1-dz)
     w2=dx    *(1-dy)*(1-dz)
     w3=(1-dx)*dy    *(1-dz)
     w4=dx    *dy    *(1-dz)
     w5=(1-dx)*(1-dy)*dz
     w6=dx    *(1-dy)*dz
     w7=(1-dx)*dy    *dz
     w8=dx    *dy    *dz

     pFx=w1*F_x(i,j,k)   +w2*F_x(i+1,j,k)   +w3*F_x(i,j+1,k)   +w4*F_x(i+1,j+1,k)+&
         w5*F_x(i,j,K+1) +w6*F_x(i+1,j,K+1) +w7*F_x(i,j+1,k+1) +w8*F_x(i+1,j+1,k+1)
     pFy=w1*F_y(i,j,k)   +w2*F_y(i+1,j,k)   +w3*F_y(i,j+1,k)   +w4*F_y(i+1,j+1,k)+&
         w5*F_y(i,j,K+1) +w6*F_y(i+1,j,K+1) +w7*F_y(i,j+1,k+1) +w8*F_y(i+1,j+1,k+1)     
     pFz=w1*F_z(i,j,k)   +w2*F_z(i+1,j,k)   +w3*F_z(i,j+1,k)   +w4*F_z(i+1,j+1,k)+&
         w5*F_z(i,j,K+1) +w6*F_z(i+1,j,K+1) +w7*F_z(i,j+1,k+1) +w8*F_z(i+1,j+1,k+1)          
#endif
          
end subroutine InterpVectorField_CellCenters


!-------------------------------------------------------------------------------------------------------
! The following interpolation method used vecE,B matrices for interpolation, similar to how the interpolation
! is done in movdep.F90
!------------------------------------------------------------------------------------------------------- 

subroutine Interp_VecEM_GridPoints(x,y,z,pEx,pEy,pEz,pBx,pBy,pBz)
	real(psn), intent(IN)    :: x,y,z
	real(psn), intent(INOUT) :: pEx,pEy,pEz,pBx,pBy,pBz
	real(psn) :: dx,dy,dz
	integer   :: ip0, jp0, kp0, short_arr_ind, nn
#ifdef twoD 
    real(psn), dimension(4) :: wt
#else
    real(psn), dimension(8) :: wt 
#endif
	
	ip0=x
    jp0=y
			  
	dx=x-ip0
	dy=y-jp0 

#ifndef twoD 
    kp0=z     
    dz=z-kp0 	
!$OMP SIMD
    do nn=1,8
		wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
    end do 	  

#else
!$OMP SIMD
	do nn=1,4 
		wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
	end do     
#endif   

#ifdef twoD 
    short_arr_ind= (jp0-1)*mx + ip0 
#else
    short_arr_ind= (kp0-1)*mx*my  + (jp0-1)*mx + ip0 
#endif
	
	pEx=0.0_psn; pEy=0.0_psn; pEz=0.0_psn; pBx=0.0_psn; pBy=0.0_psn; pBz=0.0_psn;	  
		  
!$OMP SIMD
#ifdef twoD
	do nn=1,4
#else
    do nn=1,8
#endif
		  pEx=pEx+wt(nn)*VecEx(nn,short_arr_ind)
		  pEy=pEy+wt(nn)*VecEy(nn,short_arr_ind)
		  pEz=pEz+wt(nn)*VecEz(nn,short_arr_ind)
		  pBx=pBx+wt(nn)*VecBx(nn,short_arr_ind)
		  pBy=pBy+wt(nn)*VecBy(nn,short_arr_ind)
		  pBz=pBz+wt(nn)*VecBz(nn,short_arr_ind)

	end do
	
	
end subroutine Interp_VecEM_GridPoints

!-------------------------------------------------------------------------------------------------------
! The following interpolation method  is used in "savedata_routines" 
! it uses existing VecEx, e.t.c .. to interpolated fld quantities at the particle's location 
!------------------------------------------------------------------------------------------------------- 

subroutine Interp_VecFld_GridPoints(x,y,z,pFx,pFy,pFz,vecFx,vecFy,vecFz)
	real, intent(IN)          :: x,y,z
	real, intent(INOUT)       :: pFx,pFy,pFz
	real(psn), dimension(:,:) :: vecFx,vecFy,vecFz
	real(psn) :: dx,dy,dz
	integer   :: ip0, jp0, kp0, short_arr_ind, nn
#ifdef twoD 
    real(psn), dimension(4) :: wt
#else
    real(psn), dimension(8) :: wt 
#endif
	
	ip0=x
    jp0=y

	dx=x-ip0
	dy=y-jp0 

#ifndef twoD  
	kp0=z    
    dz=z-kp0 	
!$OMP SIMD
    do nn=1,8
		wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)*(wtz1(nn)+wtz2(nn)*dz)
    end do 	  

#else
!$OMP SIMD
	do nn=1,4 
		wt(nn)=(wtx1(nn)+wtx2(nn)*dx)*(wty1(nn)+wty2(nn)*dy)
	end do     
#endif   

#ifdef twoD 
    short_arr_ind= (jp0-1)*mx + ip0 
#else
    short_arr_ind= (kp0-1)*mx*my  + (jp0-1)*mx + ip0 
#endif
	
	pFx=0.0_psn; pFy=0.0_psn; pFz=0.0_psn; 
		  
!$OMP SIMD
#ifdef twoD
	do nn=1,4
#else
    do nn=1,8
#endif
		  pFx=pFx+wt(nn)*VecFx(nn,short_arr_ind)
		  pFy=pFy+wt(nn)*VecFy(nn,short_arr_ind)
		  pFz=pFz+wt(nn)*VecFz(nn,short_arr_ind)

	end do
	
	
end subroutine Interp_VecFld_GridPoints
     

end module interpolation 