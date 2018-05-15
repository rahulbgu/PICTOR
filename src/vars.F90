module vars
     use parameters
#ifdef OPEN_MP	 
	 USE OMP_LIB
#endif  
     use hdf5
     implicit none
     
     integer, parameter :: dbpsn=kind(1.0d0)
     integer :: mx,my,mz
     
     real(psn),dimension(:,:,:), allocatable :: Ex,Ey,Ez,Bx,By,Bz
     real(psn),dimension(:,:,:), allocatable :: Jx,Jy,Jz
     real(psn) ::Bx_ext0,By_ext0,Bz_ext0! constant external magnetic field 
     real(psn),dimension(:,:,:), allocatable :: FilteredEx,FilteredEy,FilteredEz
     !----------------------------------------------------
     ! Short Arrays used to process ordered particles in a vectorized way
     !----------------------------------------------------
	 integer :: ShortFldArrSizeY,ShortFldArrSizeZ
	 !real(psn), dimension(:,:), allocatable ::VecJ, VecEM! VecEx,VecEy,VecEz,VecBx,VecBy,VecBz,
	 real(psn), dimension(:,:), allocatable :: VecEx,VecEy,VecEz,VecBx,VecBy,VecBz
	 real(psn), dimension(:,:,:), allocatable :: VecJ
	 
     !----------------------------------------------------
     ! Array used for particles and test particles 
     !----------------------------------------------------
	 integer, dimension(:), allocatable:: flvp,tagp
	 real(psn), dimension(:), allocatable:: qp,xp,yp,zp,up,vp,wp,var1p 
	 integer, dimension(:), allocatable:: flvtp,tagtp
	 real(psn), dimension(:), allocatable:: qtp,xtp,ytp,ztp,utp,vtp,wtp,var1tp 
	 integer, dimension(:), allocatable:: flvp_temp,tagp_temp
	 real(psn), dimension(:), allocatable:: qp_temp,xp_temp,yp_temp,zp_temp,up_temp,vp_temp,wp_temp,var1p_temp
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
		  integer   :: flv
		  integer   :: tag
     end type particle
     !in case the above partile datatype changes, the following must be modified : 
     !a) mpi datatype for particle transfer b) CreateFreeSlot c) AppendParticle 
     !---------------------------------------------------------------------------
     !variables used in particle array declarations
    !---------------------------------------------------------------------------
     type(particle), dimension(:),allocatable :: p,ptemp     
     !some parameters
     real(psn) :: qmi,qme,qe,qi,ompe,masse,massi
     real(psn), dimension(:),allocatable :: flvrqm
     integer  , dimension(:),allocatable :: FlvrSaveFldData
     integer  , dimension(:),allocatable :: FlvrType! By default, 0: simulation particle 1: Test Particle
     integer  , dimension(:),allocatable :: FlvrSplitPrtl 
     integer        :: Nflvr ! Number of plasma species 
     integer        :: Nelc  !Initial number of electrons, defined in initialise.F90  
#ifdef longX
    integer, parameter :: xpsn=dbpsn
#else
    integer, parameter :: xpsn=psn
#endif 
#ifdef longY
    integer, parameter :: ypsn=dbpsn
#else
    integer, parameter :: ypsn=psn
#endif
#ifdef longZ
    integer, parameter :: zpsn=dbpsn
#else
    integer, parameter :: zpsn=psn
#endif
     
     !----------------------------------------------------
     !variables used to save data 
     !----------------------------------------------------
     integer, dimension(:), allocatable :: CurrentTagID,TagCounter !UsedTagID0 is not needed
     integer                            :: NtagProcLen=100
     real     ,dimension(:),allocatable :: pdata_real
     integer  ,dimension(:),allocatable :: pdata_int
	 real, dimension(:), allocatable    :: outq,outx,outy,outz,outu,outv,outw,outvar1
	 real, dimension(:), allocatable    :: outtag,outflv
	 real, dimension(:), allocatable    :: outpEx,outpEy,outpEz,outpBx,outpBy,outpBz
     real(psn),dimension(:,:),allocatable :: pdata_local_field
#ifdef twoD     
     integer,dimension(0:nSubDomainsX*nSubDomainsY-1) :: prtl_arr_size_all
     integer,dimension(0:nSubDomainsX*nSubDomainsY-1,1:3) :: fld_size_all
#else 
    integer,dimension(0:nSubDomainsX*nSubDomainsY*nSubDomainsZ-1) :: prtl_arr_size_all
    integer,dimension(0:nSubDomainsX*nSubDomainsY*nSubDomainsZ-1,1:3) :: fld_size_all
#endif     
     integer :: tosave_prtl_arr_size
     real,dimension(:,:,:),allocatable :: fdata
	 real,dimension(:,:,:), allocatable :: outEx,outEy,outEz,outBx,outBy,outBz,outJx,outJy,outJz,outJx1,outJy1,outJz1,outJx2,outJy2,outJz2   
     integer :: fdatax,fdatay,fdataz
     real(psn), dimension(4) :: energy
     integer :: fdataxi_box,fdataxf_box,fdatayi_box,fdatayf_box,fdatazi_box,fdatazf_box !limits for saving fld data
     !----------------------------------------------------
     !indices and variables
     !----------------------------------------------------
     integer :: t !counter for time steps  
     integer :: tstart ! starting time step number  
     !integer :: i,j,k,l,m,n
     integer :: prtl_arr_size! size of particle array
     integer :: used_prtl_arr_size !maximum index of an active particle in the particle array
	 integer :: prtl_random_insert_index ! index to insert an random particle that is out of spatial order 
	 integer, dimension(:,:), allocatable :: SortedPrtlCountYZ ! it stores number of particles in columns paralles to x-axis  
     integer :: np ! total number of active test particles on this processor
     integer :: test_prtl_arr_size! size of test particle array
     integer :: used_test_prtl_arr_size !maximum index of an active test particle in the particle array
	 integer :: test_prtl_random_insert_index ! index to insert an random test particle that is out of spatial order 
     integer :: ntp ! total number of active test particles on this processor  
     !------------------------------------------------------
    !varaibles to set the boundaries
     !------------------------------------------------------- 
     real(psn) :: xmin,xmax,ymax,ymin,zmax,zmin,xlen,ylen,zlen ! physical boundaries of particles at this processor 
     !----------------------------------------------------------
     !variable used for communication 
     !----------------------------------------------------------
     integer :: ierr,proc,nproc,lproc,rproc,tproc,bproc,uproc,dproc!l=left,r=right,t=top,b=bottom,3D: u=up(+z),d=down(-z)
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
     !varaibles used for particle exchange between machines 
     !---------------------------------------------------------
     integer :: lcross,rcross,tcross,bcross,ucross,dcross
     integer :: lpcross,rpcross,tpcross,bpcross,upcross,dpcross
     integer :: ltpcross,rtpcross,ttpcross,btpcross,utpcross,dtpcross
     integer :: linp_count,rinp_count,tinp_count,binp_count,uinp_count,dinp_count ! no. of incoming particle 
     integer :: lintp_count,rintp_count,tintp_count,bintp_count,uintp_count,dintp_count ! np. of incoming test particles 
     type(particle), dimension(:), allocatable :: linp,rinp,tinp,binp,uinp,dinp ! same array is used for both partilces and test particles
     type(particle), dimension(:), allocatable :: loutp,routp,toutp,boutp,uoutp,doutp
	 integer :: loutp_size,routp_size,toutp_size,boutp_size,uoutp_size,doutp_size     
     integer :: linp_size,rinp_size,tinp_size,binp_size,uinp_size,dinp_size     
     !---------------------------------------------------------------
     !variables used for Current transfer 
     !---------------------------------------------------------------
     real(psn),dimension(:,:,:), allocatable :: buff_lJx,buff_rJx,buff_lJy,buff_rJy,buff_lJz,buff_rJz
     real(psn),dimension(:,:,:), allocatable :: buff_bJx,buff_tJx,buff_bJy,buff_tJy,buff_bJz,buff_tJz
     real(psn),dimension(:,:,:), allocatable :: buff_dJx,buff_uJx,buff_dJy,buff_uJy,buff_dJz,buff_uJz 
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
     real(kind=8),dimension(50) :: exec_time 
     
     !----------------------------------------------------------------
     !some auxiliary variables 
     !----------------------------------------------------------------
     real(psn), parameter :: fldc=c*EMspeedBoost
     real(psn), parameter :: fld_halfc=c*0.5_psn*EMspeedBoost
     real(psn), parameter :: cinv=1.0_psn/c
     real(psn), parameter :: sqc=c**2
     real(psn),parameter :: pi=3.14159265358979323846264_psn
     real(psn), parameter :: SplitQ_MIN=1.0_psn/8.0_psn
     !-------------------------------------------------------------------
     !variables used to regualte memory usage 
     !-------------------------------------------------------------------
     integer, parameter :: size_hist_bin_size=150
     integer, dimension(size_hist_bin_size) :: np_history=0
     integer ::  hist_ind=1
     !-------------------------------------------------------------------
     ! to keep track of physical domain on each proc
     !-------------------------------------------------------------------
     integer , dimension (0:nSubDomainsX) :: xborders,xborders_new
     integer , dimension (0:nSubDomainsY) :: yborders,yborders_new
#ifdef twoD     
    integer , dimension (0:1)      :: zborders
     integer , dimension (0:nSubDomainsX*nSubDomainsY-1) ::procxind,procyind,proczind
#else 
    integer , dimension (0:nSubDomainsZ) :: zborders
    integer , dimension (0:nSubDomainsX*nSubDomainsY*nSubDomainsZ-1) ::procxind,procyind,proczind 
#endif     
     !-------------------------------------------------------------------
     !varaibles used for the purpose of load balancing 
     !-------------------------------------------------------------------
#ifdef twoD     
    integer,dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:0) :: proc_grid !indices of these matrices correspond to actual physical location 
#else 
    integer,dimension(0:nSubDomainsX-1,0:nSubDomainsY-1,0:nSubDomainsZ-1) :: proc_grid !indices of these matrices correspond to actual physical location 
#endif        
     integer :: loadbalance_profiling_period=1
     real, dimension(100) :: hist_time_movdep=0,hist_time_comm=0
     real                 :: mean_time_movdep  ,mean_time_comm 
	 
     !-------------------------------------------------------------------
     !OpenMP related varaibles 
     !-------------------------------------------------------------------
	 integer :: Nthreads,ThreadID
  
end module vars
