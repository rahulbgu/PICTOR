module parameters
!----------------------------------
!all user defined input parameters are listed here
!----------------------------------
!integer, parameter :: psn=selected_real_kind(15,307)  
integer, parameter :: psn=kind(1.0e0)!precision of the simulation 
integer, parameter :: nx=1024!32*8!*16*8!no. of grid cells along x  
integer, parameter :: ny=512!32*8!no. of grid cells along y 
integer, parameter :: nz=64!*16!no. of grid cells along z (relevant only for 3D simulation)
integer, parameter :: nSubDomainsX=2!*16!*2!32 ! no. of sub-domains along x
integer, parameter :: nSubDomainsY=2!8!*2! no. of sub-domains along y
integer, parameter :: nSubDomainsZ=1! no. of sub-domains along z (relevant only for 3D)
integer, parameter :: tfinish=5000! total number of time steps

! physical parameters 
real(psn), parameter :: g0=1.0_psn!if a relativistic plasma, g0 is the Lorentz factor of particles (it only sets the electron skin-depth) 
integer, parameter   :: epc=4!number of the electrons per cell 
real(psn), parameter :: mi=32.0_psn ! mass of ions
real(psn), parameter :: me=1.0_psn ! mass of electrons
real(psn), parameter :: c_ompe=8.0_psn ! no. of cells per electron skin-depth 
real(psn), parameter :: c=0.45_psn ! speed of light (photons propagate about a distance c in one time step)
real(psn), parameter :: EMspeedBoost=1.0_psn! EM waves are propagated faster than c by this factor 
integer , parameter  :: curr_filter=2!3!no. of times the urrent filter (moving average) is applied  
integer , parameter  :: nMoverEMfilter=0!8!no. of times the EM fld averaging is done, only relavant for the test particles 

!parameters for saving the output
!b) Set the period to "-1" to disable saving a particular data, e.g. "fld_save_period=-1" will disble saving the field data   
character(len=256) :: data_folder=trim("data") !name of the folder to save the output files
integer            :: fld_save_period=100!in time steps, saves Ex,Ey,Ez,Bx,By,Bz by default 
integer            :: prtl_save_period=100!default: spatial cordinates: x,y,z;three-velocity: u,v,w; and q(charge),flv(species type)
integer            :: spec_save_period=100!energy spectrum of the particles
integer, parameter :: psave_ratio=64!(default ratio) 1 in every "psave_ratio" particles are tagged for tracking
integer            :: fsave_ratio=4!(default downsampling) field qunatities are saved in a reduced resolution by this factor
logical            :: save_prtl_local_fld=.true.! to save the E abd B field at the particle's location 
logical            :: save_prtl_local_curr=.true.! save the value of local current at the particle's location 
logical            :: save_density=.true. ! density of particles
logical            :: save_ch_flux_fld=.true. ! flux of the charge(q), i.e., q<v_i>, i=x,y,z
logical            :: save_velsq_fld=.false. ! q<v^2_i>, i=x,y,z
logical            :: save_tot_curr=.true. ! the net current, the exact quantiites used in the Maxwell's equations
logical            :: save_divE=.false. ! divergence of the Electric field
logical            :: save_EdotV_fld=.false. ! <E.V>, averaging is done over a bigger cells of size fsave_ratio
logical            :: save_eng_spat=.false. ! total energy of particle in a cell of size fsave_ratio
!The following parameters are useful to save various useful qunatities 
integer            :: prtl_mean_save_period=-1 !serval mean quantitities(e.q., E.V), but now averaged over all the particles

!------------------------------------------------------------------------------------------------
!Advanced settings
!------------------------------------------------------------------------------------------------

!parameters for restarting the simulation from a saved state
integer, parameter :: restart_save_period=100000000 ! period to save state of the simulation  
logical, parameter :: restart=.false. !set to true if the simulation is to be restarted from a saved data, false otherwise
integer, parameter :: restart_time=3000 ! time step no. to load the restart data  

! parameters for performance optimisation
integer, parameter :: prtl_reorder_period=20!sort particles in order to improve cache performance 
integer, parameter :: performance_save_period=300000 !period to save performance of each core
integer, parameter :: domain_change_period=120!(optional) time period to change the size of sub-domains to balance the load  

!parameters for computing particle spectrum
integer,parameter  :: prtl_spec_type=3 !1: only gamma spectrum, 2: only speed spectrum, 3 :both, gamma and speed   
real   ,parameter  :: Gamma_spec_binwidth=0.1/10 !dlog(Gamma)
real   ,parameter  :: Speed_spec_binwidth=0.001 !d speed 

!parameters for GPU 
integer, parameter :: gpu_prtl_arr_size=110000000 ! max. no. of particles on GPU
integer, parameter :: gpu_test_prtl_arr_size=1000000 ! max. no. of test particles on GPU
integer            :: buff_size_prtl_gpu=10000000 ! size of the arrays that temporarily hold prtls coming/leaving the sub-domain
integer            :: buff_size_test_prtl_gpu=100000 ! same as above but for the test particles

end module parameters
