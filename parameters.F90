module parameters
!----------------------------------
!all user defined input parameters are listed here
!---------------------------------- 
integer, parameter :: psn=kind(1.0d0)	!computation precision; 1.0e0: single precision, 1.0d0: double precision  
integer, parameter :: nx=2048		!no. of grid cells along x  
integer, parameter :: ny=256	    !no. of grid cells along y 
integer, parameter :: nz=4000		!no. of grid cells along z (used only if THREE_DIMS=y (3D simulation) in the Makefile)
integer            :: nSubDomainsX=2	!*16!*2!32 ! no. of sub-domains along x
integer            :: nSubDomainsY=4	!8!*2! no. of sub-domains along y
integer            :: nSubDomainsZ=1	! no. of sub-domains along z (relevant only for 3D)
integer, parameter :: tfinish=200000	! maximum number of time steps
real(psn),parameter:: saveAndExitAfterHours = 2.0 !in hours (>0) :  stop execution, save restart data, and exit 

! physical parameters 
real(psn), parameter :: g0=1.0_psn	!if a relativistic plasma, g0 is the Lorentz factor of particles (determines the electron skin-depth) 
integer, parameter   :: epc=1		!number of electrons per unit cell of volume = 1 (dx=dy=dz=1), determines the electron charge
real(psn), parameter :: mi_me=3671.0_psn!mass ratio = ion mass/electron mass = mi/me 
real(psn), parameter :: c_ompe=30.0_psn! no. of cells per electron skin-depth 
real(psn), parameter :: c=0.45_psn!0.98*0.707_psn! the speed of light (i.e., EM waves propagate about a distance c in one time step)
real(psn), parameter :: EMspeedBoost=1.0_psn! EM waves are propagated faster than c by this factor 
integer , parameter  :: curr_filter=0	!3!no. of times the current filter (moving average) is applied  
integer , parameter  :: nMoverEMfilter=0!8!no. of times the EM fld averaging is performed, only relavant for the test particles 
real(psn),parameter  :: ne0=1.20e19 ! (optional) electron density in cm^-3, sets physical units 

!------------------------------------------------------------------------------------------------
! Grid parameters
!------------------------------------------------------------------------------------------------
!grid spacing
real(psn), parameter :: grid_dx = 1.0_psn
real(psn), parameter :: grid_dy = 1.0_psn
real(psn), parameter :: grid_dz = 1.0_psn

!1D electrons per cell, (defines the 'even placement' of particles in a cell; epc = epc_x*epc_y*epc_z) 
integer, parameter   :: epc_x = 1
integer, parameter   :: epc_y = 1
integer, parameter   :: epc_z = 2

!------------------------------------------------------------------------------------------------
! Saving the output
!------------------------------------------------------------------------------------------------
! Set the period to "-1" to disable saving a particular data, e.g. "fld_save_period=-1" will disble saving the field data   
character(len=256) :: data_folder=trim("data") 	!name of the folder to save the output files
integer            :: fld_save_period=1	!in time steps, saves Ex,Ey,Ez,Bx,By,Bz by default 
integer            :: prtl_save_period=100	!00000!default: spatial cordinates: x,y,z;three-velocity: u,v,w; and q(charge),flv(species type)
integer            :: spec_save_period=100	!energy spectrum of the particles
integer, parameter :: psave_ratio=1	!(default ratio) 1 in every "psave_ratio" particles are tagged for tracking
integer            :: fsave_ratio=1		!(default downsampling) field qunatities are saved in a reduced resolution by this factor
logical            :: save_prtl_local_fld=.true.! to save the E abd B field at the particle's location 
logical            :: save_prtl_local_curr=.false.! save the value of local current at the particle's location 
logical            :: save_density=.true. 	! density of particles
logical            :: save_ch_flux_fld=.true. 	! flux of the charge(q), i.e., q<v_i>, i=x,y,z
logical            :: save_velsq_fld=.false. 	! q<v^2_i>, i=x,y,z
logical            :: save_tot_curr=.true. 	! the net current, the exact quantiites used in the Maxwell's equations
logical            :: save_divE=.false. 	! divergence of electric field
logical            :: save_EdotV_fld=.false. 	! <E.V>, averaging is done over cells of size fsave_ratio
logical            :: save_eng_spat=.false. 	! total energy of particles in cells of size fsave_ratio
!The following parameters are useful to save various useful qunatities 
integer            :: prtl_mean_save_period=-1 	!several mean quantitities(e.q., E.V), but now averaged over all the particles
integer, parameter :: restart_save_period=1000	! period to save state of simulation (checkpoint) to restart  

!------------------------------------------------------------------------------------------------
!Cylinderical grid (coordiante system must be changed in "Makefile" as well, default is Cartesian)  
!------------------------------------------------------------------------------------------------
logical, parameter :: inc_axis=.true. 	! Incude the axis in the domain, no inner BC 
logical, parameter :: RZtwoD=.false. 	!no derivative along the theta-direction; RZ -plane only; see cyl_em_update* (GPU version only)  

!------------------------------------------------------------------------------------------------
!Advanced settings  :: Performance
!------------------------------------------------------------------------------------------------
integer, parameter :: prtl_reorder_period=20		!sort particles in order to improve cache performance 
integer, parameter :: performance_save_period=100000 	!period to save performance of each core
integer, parameter :: load_balance_period=120	!(optional) time period to change the size of sub-domains to balance the load 
integer            :: load_balance_yz_period = 120 !(optional) change the number of proc. along x-axis to balance the load across yz axes

!------------------------------------------------------------------------------------------------
!Advanced settings  :: Particle spectrum
!------------------------------------------------------------------------------------------------
real   ,parameter  :: Gamma_spec_binwidth=0.1/10	!dlog(Gamma)
real   ,parameter  :: Speed_spec_binwidth=0.001 	!d speed/c;  

!------------------------------------------------------------------------------------------------
!Advanced settings  :: GPU
!------------------------------------------------------------------------------------------------
integer, parameter :: gpu_prtl_arr_size=512000000 	! max. no. of particles on GPU
logical            :: curr_filter_gpu=.false.		!.true.= filtering on GPU;.false. use the host filtering routines



end module parameters
