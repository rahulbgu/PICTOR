module bc_domain_size
	use parameters
	use vars
	use mem_fld
	use bc_exterior_field
	use subdomains

    implicit none
	integer, parameter :: domain_size_change_period = 120 ! domain size is adjusted every this many steps 
	integer, parameter :: expanding_domain_buff_cells = int(domain_size_change_period*c*1.1)+8 !mutiple of 4	
contains
	
	subroutine AutoResizeDomain

	    if(modulo(t,domain_size_change_period).ne.0) return
		
		if(indepLBaxis.eq.0) then
			if(bc_face(1)%speed.ne.0) call ResizeDomainLeft(0,nSubDomainsX,xborders,procxind)
			if(bc_face(2)%speed.ne.0) call ResizeDomainRight(0,nSubDomainsX,xborders,procxind)
		end if
		if(indepLBaxis.eq.1) then
			if(bc_face(3)%speed.ne.0) call ResizeDomainLeft(1,nSubDomainsY,yborders,procyind)
			if(bc_face(4)%speed.ne.0) call ResizeDomainRight(1,nSubDomainsY,yborders,procyind)
		end if
		if(indepLBaxis.eq.2) then
			if(bc_face(5)%speed.ne.0) call ResizeDomainLeft(2,nSubDomainsZ,zborders,proczind)
			if(bc_face(6)%speed.ne.0) call ResizeDomainRight(2,nSubDomainsZ,zborders,proczind)
		end if
		
	end subroutine AutoResizeDomain
	
	subroutine InitialResizeDomain
		if(indepLBaxis.eq.0) then
			call ResizeDomainLeft(0,nSubDomainsX,xborders,procxind)
			call ResizeDomainRight(0,nSubDomainsX,xborders,procxind)
		end if 
		if(indepLBaxis.eq.1) then
			call ResizeDomainLeft(1,nSubDomainsY,yborders,procyind)
			call ResizeDomainRight(1,nSubDomainsY,yborders,procyind)
		end if
		if(indepLBaxis.eq.2) then
			call ResizeDomainLeft(2,nSubDomainsZ,zborders,proczind)
			call ResizeDomainRight(2,nSubDomainsZ,zborders,proczind)
		end if
	end subroutine InitialResizeDomain
	
	
	subroutine ResizeDomainRight(axis,nSubDomains,borders,procind)
		 integer :: axis, nSubDomains, procind
		 integer, dimension(0:nSubDomains) :: borders
		 real(dbpsn) :: xb !max. position of the boundary effects 
		 integer :: inc, x1, fid
		 
		 fid = 2*axis+2
		 
		 xb = bc_face(fid)%pos_fld + bc_face(fid)%attn_thickness
		 if(bc_face(fid)%type_prtl.ne.'prdc') xb = max(xb,bc_face(fid)%pos_prtl)
		 		 
		 inc = 0 
		 if(bc_face(fid)%speed.gt.0) then 
			 call CalcExpansionIncRight(inc, xb, borders(nSubDomains) )
		 else 
			 !Note: currently subdomains at the boundary are not allowed to shrink independently 
			 ! box remiains a cuboid; otherwise at least "savedata" and this routine need to become more general 
			 x1 = borders(nSubDomains-1)
			 call MPI_AllREDUCE(MPI_IN_PLACE,x1,1,MPI_INTEGER,mpi_max, MPI_COMM_WORLD)
			 xb = min(real(borders(nSubDomains),dbpsn),xb)
			 xb = max(real(x1,dbpsn),xb)
			 if(xb.lt.borders(nSubDomains)) call CalcContractionInc(inc, xb , borders(nSubDomains) )
			 
		 end if		 
		 
		 borders(nSubDomains) = borders(nSubDomains) + inc
		 if(indepLBaxis.eq.axis) ProcGrid(iproc,jproc)%borders(nSubDomains) = borders(nSubDomains) 
		 
		 if( inc.ne.0 .and. procind.eq.nSubDomains-1 ) then 
			 
#ifdef gpu
	         	call RecvFullDomainEMFldFromGPU
#endif		
			 	if(axis.eq.0) call ResizeFld(mx,my,mz,mx+inc,my,mz,0,0,0)
				if(axis.eq.1) call ResizeFld(mx,my,mz,mx,my+inc,mz,0,0,0)
				if(axis.eq.2) call ResizeFld(mx,my,mz,mx,my,mz+inc,0,0,0)
				
				call SetSubDomainGridSize
				call SetPrtlBoundaries
				call InitAuxFld
				
				!Initialise EM Fld in the newly created domain
				if(bc_face(fid)%type_fld.eq.'iflw') call Set_BC_EM_FlowField(fid)			 
#ifdef gpu
			 	call ResetGridGPU
#endif 	
				
		 end if
		 
	     call MPI_Barrier(MPI_COMM_WORLD) 
		 
		 !In general, all subdomains need to reset their neighbors since size at a ngbr might have changed independently
		 call ResetCommSubdomains 
		
		 box_bounds(fid)=borders(nSubDomains) !Note : in general, may need to find global min,max when needed , e.g. before writing data
		 		 
	end subroutine ResizeDomainRight
	
	
	subroutine ResizeDomainLeft(axis,nSubDomains,borders,procind)
		 integer :: axis, nSubDomains, procind
		 integer, dimension(0:nSubDomains) :: borders 
		 real(dbpsn) :: xb !min. position of the boundary effects 
		 integer :: inc, x1, fid
		 
		 fid = 2*axis+1
				 
		 xb = bc_face(fid)%pos_fld - bc_face(fid)%attn_thickness
		 if(bc_face(fid)%type_prtl.ne.'prdc') xb = min(xb, bc_face(fid)%pos_prtl)
		 
		 inc=0
		 if(bc_face(fid)%speed.lt.0) then 
			 call CalcExpansionIncLeft(inc, xb, borders(0) )
		 else 
			 !Note: currently subdomains at the boundary are not allowed to shrink independently 
			 ! box remiains a cuboid; otherwise at least "savedata" and this routine need to become more general 
			 x1 = borders(1)
			 call MPI_AllREDUCE( MPI_IN_PLACE, x1, 1, MPI_INTEGER, mpi_min, MPI_COMM_WORLD)
			 xb = max(real(borders(0),dbpsn),xb)
			 xb = min(real(x1,dbpsn),xb)
			 if(xb.gt.borders(0)) call CalcContractionInc(inc, min(real(x1,dbpsn),xb), borders(0) )
		 end if
		 
		 		 		 
		 borders(0) = borders(0) - inc
		 if(indepLBaxis.eq.axis) ProcGrid(iproc,jproc)%borders(0) = borders(0) 
		 	 
		 if( inc.ne.0 .and. procind.eq.0 ) then 
			 
#ifdef gpu
	         	call RecvFullDomainEMFldFromGPU
#endif		
			 	if(axis.eq.0) call ResizeFld(mx,my,mz,mx+inc,my,mz,inc,0,0)
				if(axis.eq.1) call ResizeFld(mx,my,mz,mx,my+inc,mz,0,inc,0)
				if(axis.eq.2) call ResizeFld(mx,my,mz,mx,my,mz+inc,0,0,inc)
				
				call SetSubDomainGridSize
				call SetPrtlBoundaries
				call InitAuxFld
				
				if(axis.eq.0) call ShiftPrtlPos(real(inc,psn),xmin,xp,flvp)
				if(axis.eq.1) call ShiftPrtlPos(real(inc,psn),ymin,yp,flvp)
				if(axis.eq.2) call ShiftPrtlPos(real(inc,psn),zmin,zp,flvp)
				
				!Initialise EM Fld in the newly created domain
				if(bc_face(fid)%type_fld.eq.'iflw') call Set_BC_EM_FlowField(fid)			 
#ifdef gpu
			 	call ResetGridGPU
#endif 				
				
		 end if 
		 
		 call MPI_Barrier(MPI_COMM_WORLD) 
		 
		 call ResetCommSubdomains
		 
		 box_bounds(fid)=borders(0) !Note : in general, may need to find global min,max when needed , e.g. before writing data
	end subroutine ResizeDomainLeft
	
	subroutine ResetCommSubdomains
		 call SetNgbr
	     call SetSendRecvFldDomain
	     call SetSendRecvPrtl
		 call SetBoundaryNgbr
	end subroutine ResetCommSubdomains
	
    subroutine ShiftPrtlPos(shift,xmin,x,flv)
		real(psn) , dimension(:) :: x 
		integer, dimension(:) :: flv
		real(psn) :: shift, xmin
		integer  :: n
		do n=1,used_prtl_arr_size
		    x(n) = x(n)+shift
			if(flv(n).eq.0) x(n) = xmin ! x-pos of empty prtl slots need to be reset 
	    end do
    end subroutine ShiftPrtlPos
	
		
	subroutine CalcExpansionIncRight(inc,xb,x2)
		integer :: inc, x2
		real(dbpsn) :: xb
		
		if( x2 .lt. xb + expanding_domain_buff_cells) then
			inc = expanding_domain_buff_cells 
		end if 

	end subroutine CalcExpansionIncRight
	
	subroutine CalcExpansionIncLeft(inc,xb,x2)
		integer :: inc, x2
		real(dbpsn) :: xb
		
		if( xb  .lt. x2 + expanding_domain_buff_cells) then
			inc = expanding_domain_buff_cells 
		end if 

	end subroutine CalcExpansionIncLeft
	
	subroutine CalcContractionInc(inc,xb,x2)
		integer :: inc, x2
		real(dbpsn) :: xb
		
		if( abs(x2-xb) .gt. 16) then
			inc = -4*floor(abs(x2-xb)/4.0) +4 !reduce the subdomain size, in mutiple of 4
		end if 
		
	end subroutine CalcContractionInc
	
end module bc_domain_size