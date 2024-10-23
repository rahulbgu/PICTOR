module vars
     use parameters
#ifdef OPEN_MP	 
	 USE OMP_LIB
#endif  
     use hdf5
	 use mpi_f08 
     implicit none
     
     integer, parameter   :: dbpsn=kind(1.0d0)
	 logical              :: restart=.false. !is changed to .true. if a saved restart data is found
     integer              :: mx,my,mz
	 real(psn), parameter :: me=1.0_psn ! dimensionless electron mass (normalised to the electron mass) of electrons = 1 by choice
	 real(psn), parameter :: mi=mi_me ! dimensionless mass of ions (normalised to the electron mass) = mi_me
	 real(psn), parameter :: grid_inv_dx = 1.0_psn/grid_dx
	 real(psn), parameter :: grid_inv_dy = 1.0_psn/grid_dy
	 real(psn), parameter :: grid_inv_dz = 1.0_psn/grid_dz
	 real(psn), parameter :: cell_volume = grid_dx*grid_dy*grid_dz
     
     !----------------------------------------------------
     ! Grid variables 
     !----------------------------------------------------
	 real(psn),dimension(:,:,:), allocatable :: Ex,Ey,Ez,Bx,By,Bz
     real(psn),dimension(:,:,:), allocatable :: Jx,Jy,Jz
     real(psn) ::Bx_ext0,By_ext0,Bz_ext0! constant external magnetic field 
     real(psn),dimension(:,:,:), allocatable :: FilteredEx,FilteredEy,FilteredEz

	 type :: vec_fld
	 	real(psn), dimension(:,:,:), allocatable  :: x
		real(psn), dimension(:,:,:), allocatable  :: y
		real(psn), dimension(:,:,:), allocatable  :: z
	 end type vec_fld
     !----------------------------------------------------
     ! Short Arrays used to process ordered particles in a vectorized way
     !----------------------------------------------------
	 integer, parameter :: VecBlockSize=32
	 real(psn), dimension(:,:), allocatable :: VecEx,VecEy,VecEz,VecBx,VecBy,VecBz
	 real(psn), dimension(:,:,:), allocatable :: VecJ
	 
     !----------------------------------------------------
     ! Array used for particles and test particles 
     !----------------------------------------------------
	 integer, dimension(:), allocatable:: flvp,tagp, procp
	 real(psn), dimension(:), allocatable:: qp,xp,yp,zp,up,vp,wp,var1p
	 real(psn), dimension(:), allocatable:: qmp, wtp 

	 integer, dimension(:), allocatable:: flvp_temp,tagp_temp, procp_temp
	 real(psn), dimension(:), allocatable:: qp_temp,xp_temp,yp_temp,zp_temp,up_temp,vp_temp,wp_temp,var1p_temp
	 real(psn), dimension(:), allocatable:: qmp_temp, wtp_temp 
     
	 !----------------------------------------------------
     !particle datatype : used for transferring particles 
     !----------------------------------------------------
     type :: particle
		  sequence

          real(psn) :: q		  
          real(psn) :: x
          real(psn) :: y
          real(psn) :: z
          real(psn) :: u
          real(psn) :: v
          real(psn) :: w
		  real(psn) :: var1 !multipurpose variable, problem dependent 
		  
		  !cases where each particle has distinct q/m and weight :: make optinal?
		  real(psn) :: qm !charge to mass ratio
		  real(psn) :: wt !particle weight

		  integer   :: flv

		  !particle tracking :: make optional
		  integer   :: tag
		  integer   :: proc

		  !dummy : to ensure data alignments in double precision. 
		  integer   :: dum
		  
	 end type particle
     !in case the above partile datatype changes, the following must be modified : 
     !a) mpi datatype for particle transfer b) AppendParticle 
     
	 !---------------------------------------------------------------------------
     ! charge, mass of the computational electrons
     !---------------------------------------------------------------------------
	 real(psn), parameter :: ompe = c/c_ompe
	 real(psn), parameter :: qi=(g0*ompe**2)/((real(epc,dbpsn)/cell_volume)*(1.0_psn+me/mi))
	 real(psn), parameter :: qe=-qi
	 real(psn), parameter :: qme=-1.0_psn
	 real(psn), parameter :: qmi=1.0_psn/mi
	 real(psn), parameter :: masse=abs(qe)*me
	 real(psn), parameter :: massi=qi*mi

	 
	 !---------------------------------------------------------------------------
     !variables used to characterize particles in the main array
     !---------------------------------------------------------------------------
     !real(psn) :: qmi,qme,qe,qi,ompe,masse,massi
     real(psn), dimension(:),allocatable :: flvrqm
	 real(psn), dimension(:),allocatable :: FlvrCharge !initial charge
     integer  , dimension(:),allocatable :: FlvrSaveFldData
     integer  , dimension(:),allocatable :: FlvrType! By default, 0: simulation particle 1: Test Particle
     integer  , dimension(:),allocatable :: FlvrSaveRatio 
     integer        :: Nflvr ! Number of plasma species 

	 
	 type :: FlvrProperty
	 	integer    :: ionization = 0 ! 0 : no ionization
		integer    :: ionization_elc_flv = 0  
		integer    :: Z_nucleus = 0
	 end type FlvrProperty

	 type(FlvrProperty), dimension(:), allocatable :: flvr_prpt
     
     !----------------------------------------------------
     !variables used to save data 
     !----------------------------------------------------

	 integer, dimension(:), allocatable :: CurrentTagID ! tag counter
	 integer, dimension(:), allocatable :: CurrentTagProcID ! tag proc ID
	 integer, dimension(:) , allocatable :: prtl_arr_size_all
	 integer :: nprtl_save_this
	 real(psn) :: psave_all_gmin = huge(1.0_psn)
     
	 real,dimension(:,:,:),allocatable :: fdata
     integer :: fdatax,fdatay,fdataz
	 real(psn) :: binlen 
     integer :: fdataxi,fdatayi,fdatazi
     integer :: fdataxf,fdatayf,fdatazf
	 
     !----------------------------------------------------
     !indices and variables
     !----------------------------------------------------
     integer :: t !counter for time steps  
     integer :: tstart ! starting time step number  
	 integer :: restart_time! time step to load the restart data
     integer :: prtl_arr_size=0! size of particle array
     integer :: used_prtl_arr_size=0 !maximum index of an active particle in the particle array
     integer :: np=0 ! total number of particles on this processor
	 integer, parameter :: outp_arr_block_size=10240
	 integer            :: max_nprtl_out = 0
     !------------------------------------------------------
     !varaibles to set the bounds 
     !------------------------------------------------------- 
     real(psn) :: xmin,xmax,ymax,ymin,zmax,zmin,xlen,ylen,zlen ! physical boundaries of particles at this processor 
	 integer, dimension(6) :: box_bounds
     !------------------------------------------------------
     !varaibles to set boundary conditions for field and particles
     !------------------------------------------------------- 
	 logical :: inflowBC = .false.
	 logical :: pmlBC = .false.
	 	 
 	 type FlowFldType
 		 real(psn) :: Attenuate  
 		 procedure(vector_global), pointer, nopass  :: MagFld =>null()
 		 procedure(vector_global), pointer, nopass :: Drift =>null()
 		 real(psn) :: vmax	
 	 end type FlowFldType
	
	 type BC_Plane
		 real(dbpsn) :: pos_fld, pos_prtl
		 character (len=4) :: type_fld, type_prtl
		 real(psn) :: speed = 0.0_psn! speed of the boundary, if moving
		 real(psn) :: attn_thickness = 0.0_psn! attenuate EM fld
		 real(psn) :: attn_scale = 1.0_psn ! scale length for the attenuation of EM fld
		 type(FlowFldType) :: flw
		 !Note: inj or their indices should be an present here for better customisation and performance
	 end type BC_Plane
	 
	 type(BC_Plane), dimension(6) :: bc_face 
		  
     !----------------------------------------------------------
     !variable used for communication 
     !----------------------------------------------------------
     integer :: proc,nproc! proc = MPI rank; nproc = total number of MPI tasks/subdomains
     !---------------------------------------------------------
     !variables used for read and write (HDF5)
     !---------------------------------------------------------
     INTEGER(HID_T) :: fid, dset_id,dspace_id,memspace,plist
     INTEGER(HSIZE_T), dimension(1) :: data_dim1,max_dim1,offset1,chunk_dim1,&
                                       new_size1,dcount1
     INTEGER(HSIZE_T), dimension(3) :: data_dim3,max_dim3,offset3,chunk_dim3,new_size3
     integer :: err,rank
     character (len=1024) :: fname,cmd,str1,str2

     !---------------------------------------------------------
     !variables used for particle exchange among subdomains
     !---------------------------------------------------------

	 integer :: np_recv
     		 
	 type neighbor
		 integer :: proc
		 integer, dimension(6) :: ind
		 real(psn) :: xshift, yshift, zshift ! shift in/out particles location 
		 
		 type(particle), dimension(:), allocatable :: p 
		 integer :: pcount
		 type(MPI_Request) :: mpi_req_prtl, mpi_req_prtl_size
					  
		 real(psn), dimension(:,:,:), allocatable :: Fldx, Fldy,Fldz	 
		 type(MPI_Request) :: mpi_req_fldx, mpi_req_fldy, mpi_req_fldz 
		 type(MPI_Status)  :: mpi_stat_fldx, mpi_stat_fldy, mpi_stat_fldz
		 integer :: mpi_tag_offset = 0	 	
	 end type neighbor
	 
	 type ngbr_list
		 integer :: num_ngbrs = 0
		 type(neighbor), dimension(:), allocatable :: ngbr 
	     integer, dimension(:), allocatable :: edges 
	 end type ngbr_list
	 
	 type(ngbr_list), dimension(:), allocatable :: ngbr_send, ngbr_recv
	 type(ngbr_list), dimension(:), allocatable :: lb_send, lb_recv
		        
     !----------------------------------------------------------------
     !variables used for smoothening 
     !----------------------------------------------------------------
     real(psn),dimension(:,:,:), allocatable :: F0
     real(psn) :: wtm1=0.25_psn
     real(psn) :: wt0=0.5_psn
     real(psn) :: wtp1=0.25_psn
     
     !----------------------------------------------------------------
     !variable used to monitor performance
     !----------------------------------------------------------------
     real(dbpsn), dimension(50) :: exec_time 
     
     !----------------------------------------------------------------
     !Constants
     !----------------------------------------------------------------
     real(psn), parameter :: fldc=c*EMspeedBoost
     real(psn), parameter :: fld_halfc=c*0.5_psn*EMspeedBoost
     real(psn), parameter :: cinv=1.0_psn/c
     real(psn), parameter :: sqc=c**2
     real(psn), parameter :: pi=3.14159265358979323846264_psn
     real(psn), parameter :: SplitQ_MIN=1.0_psn/8.0_psn

	 !----------------------------------------------------------------
     ! Scaling factors to convert quantities to Physical units
     !----------------------------------------------------------------
	 real(dbpsn), parameter  :: skin_depth_cm = (5.31e+5)*(sqrt(ne0))
	 real(dbpsn), parameter  :: elc_per_prtl = (ne0*(skin_depth_cm**3))/(epc*(c_ompe**3))
	 real(dbpsn), parameter  :: grid_spacing_cm = skin_depth_cm / c_ompe  
	 real(dbpsn), parameter  :: E_SI = sqrt(ne0*9.2466e+03)/sqrt(epc*masse*c*c) !V/m

	 
     !-------------------------------------------------------------------
     ! to keep track of physical domain on each proc
     !-------------------------------------------------------------------	 
     integer , dimension(:), allocatable  :: xborders, yborders, zborders
	 integer :: procxind, procyind, proczind
	 
	 type proc_list
		 integer :: count
		 integer, dimension(:), allocatable :: borders, procs 
	 end type proc_list
	 type(proc_list), dimension(:,:), allocatable :: ProcGrid
	 integer :: isizeProcGrid, jsizeProcGrid ! size of 2D ProcGrid, the third direction (proc_list) has variable size
	 integer :: iproc, jproc, kproc ! index location of the proc in ProcGrid
	  
	 integer :: indepLBaxis = 0! 0=x,1=y,2=z; independent load balancing along this axis
     type(MPI_Comm) :: comm_indepLBaxis, comm_kproc0

	 
     !-------------------------------------------------------------------
     !OpenMP related varaibles 
     !-------------------------------------------------------------------
	 integer :: Nthreads,ThreadID
	 
     !-------------------------------------------------------------------
     !Curved BC related varaibles 
     !-------------------------------------------------------------------
	 logical:: curved_bc=.false.! becomes .true. if a curved (in the Cart. grid) boundary condition is included   
	 real(psn), dimension(:,:,:), allocatable :: e_lx, e_ly , e_lz, b_arx, b_ary, b_arz ! used in conformal FDTD schemes 
	 real(psn), dimension(:,:), allocatable :: usc_wtr_bx, usc_wtr_by, usc_wtr_bz ! wts used in the USC method, r (c)= row (column) of V matrix
	 real(psn), dimension(:,:), allocatable :: usc_wtc_bx, usc_wtc_by, usc_wtc_bz
	 integer, dimension(:,:,:), allocatable :: usc_fdtd_bx, usc_fdtd_by, usc_fdtd_bz ! 0: regular fdtd; positive integer, use the USC method, integer -> index in the wt. list 
	 real(psn), dimension(:,:,:), allocatable :: usc_norm_bx, usc_norm_by, usc_norm_bz ! matrix to store normalisations for the wts
	 real(psn), dimension(:,:,:), allocatable :: usc_db1, usc_db2
	 
     !-------------------------------------------------------------------
     !Interface defining structure of the user defined fucntions and subroutines
     !-------------------------------------------------------------------
 	 abstract interface
 	     function scalar_global(x,y,z)
 	 	     import :: psn, dbpsn
 			 real(dbpsn) :: x,y,z
 			 real(psn) :: scalar_global
 	 	 end function
 	     function scalar_global_dp(x,y,z)
 	 	     import :: dbpsn
 			 real(dbpsn) :: x,y,z
 			 real(dbpsn) :: scalar_global_dp
 	 	 end function
 	     function scalar_local(x,y,z)
 	 	     import :: psn
 			 real(psn) :: x,y,z
 			 real(psn) :: scalar_local
 	 	 end function
 		 subroutine vector_global(x,y,z,fx,fy,fz)
 			import :: psn, dbpsn
 			real(dbpsn) :: x,y,z
 			real(psn)   :: fx,fy,fz
 		 end subroutine
 	     function func1D(x)
 	 	     import :: psn
 			 real(psn) :: x
 			 real(psn) :: func1D
 	 	 end function
 	 end interface 
     !-------------------------------------------------------------------
     ! Data Type to define properties of Phase-Space for a given species 
     !-------------------------------------------------------------------
	 type :: PhaseSpaceProperty
		 integer :: Flvr
		 integer :: multiplicity=1
 		 procedure(scalar_global), pointer, nopass :: Density =>null()
 		 procedure(vector_global), pointer, nopass :: DriftVelocity =>null()
 		 procedure(scalar_global), pointer, nopass :: Temperature =>null()
 		 procedure(func1D), pointer, nopass :: SpeedDist =>null()
		 real(psn) :: Vmax
		 real(psn) :: fraction = 1.0_psn
 		 real(psn), dimension(:), allocatable :: Table, PDF_Table
	 end type PhaseSpaceProperty
	 integer :: PDF_TableSize=10000
	 
	 type(PhaseSpaceProperty), dimension(100) :: PSP_list
	 integer :: used_PSP_list_size = 0
	 
     !-------------------------------------------------------------------
     ! constants used in interpolation and deposition
     !-------------------------------------------------------------------
#ifdef twoD
     real(psn), dimension(4), parameter :: wtx1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
     real(psn), dimension(4), parameter :: wtx2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
     real(psn), dimension(4), parameter :: wty1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
     real(psn), dimension(4), parameter :: wty2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
  
     real(psn), dimension(8), parameter :: wtFx= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(8), parameter :: wtFy= (/0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(8), parameter :: wtFz= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)

  	 real(psn), dimension(8), parameter :: Jwtx1= (/1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
  	 real(psn), dimension(8), parameter :: Jwtx2= (/0.0_psn, 0.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
  	 real(psn), dimension(8), parameter :: Jwty1= (/1.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(8), parameter :: Jwty2= (/-1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
  	 real(psn), dimension(8), parameter :: Jwtz1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
  	 real(psn), dimension(8), parameter :: Jwtz2= (/-1.0_psn, -1.0_psn, -1.0_psn, -1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
#else 
  	 real(psn), dimension(8), parameter :: wtx1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
  	 real(psn), dimension(8), parameter :: wtx2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
  	 real(psn), dimension(8), parameter :: wty1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn  /)
  	 real(psn), dimension(8), parameter :: wty2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn /)
  	 real(psn), dimension(8), parameter :: wtz1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn  /)
     real(psn), dimension(8), parameter :: wtz2= (/-1.0_psn, -1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn  /)
  
  	 real(psn), dimension(12), parameter :: wtFx= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(12), parameter :: wtFy= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(12), parameter :: wtFz= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)

  	 real(psn), dimension(12), parameter :: Jwtx1= (/1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn /)
  	 real(psn), dimension(12), parameter :: Jwtx2= (/0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn /)
  	 real(psn), dimension(12), parameter :: Jwty1= (/1.0_psn, 0.0_psn, 1.0_psn, 0.0_psn,  1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn /)
  	 real(psn), dimension(12), parameter :: Jwty2= (/-1.0_psn, 1.0_psn, -1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn/)
  	 real(psn), dimension(12), parameter :: Jwtz1= (/1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 1.0_psn, 1.0_psn, 1.0_psn, 1.0_psn /)
  	 real(psn), dimension(12), parameter :: Jwtz2= (/-1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, -1.0_psn, -1.0_psn, 1.0_psn, 1.0_psn, 0.0_psn, 0.0_psn, 0.0_psn, 0.0_psn /) 
#endif 

     !-------------------------------------------------------------------
     ! variables used in the cylinderical version
     !-------------------------------------------------------------------
	 real(psn) :: dtheta,inv_dtheta ! angular resolution
	 real(psn) :: ax_perm_area
	 real(psn) :: rshift ! radial positions are obtained by shifting the x-positions by 'rshift'

	 !-------------------------------------------------------------------
     ! Ionization
     !-------------------------------------------------------------------
	 logical :: ionize_prtl = .false.

	 !-------------------------------------------------------------------
     ! Convolution PML 
     !-------------------------------------------------------------------
	 type(vec_fld) :: pml_E_psi1, pml_E_psi2
	 type(vec_fld) :: pml_B_psi1, pml_B_psi2

  
end module vars
