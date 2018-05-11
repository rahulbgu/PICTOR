#MakeFile for the PIC simulation 
SHELL=/bin/sh
#---------------------------------------------------------------------------------------
#  Change the followings to set up the PIC simulation
#---------------------------------------------------------------------------------------

SETUP=setup_alfventurb
SETUP=setup_weibel
SETUP=setup_shock
#SETUP=setup_shear_flow
#SETUP=setup_Bell
#SETUP=setup_friction
#SETUP=setup_harris_sheet


#THREE_DIMS=y# y = three spatial dimensions
#OPEN_MP=y #to include openmp support NOTE: include -qopenmp (for intel) or -fopenmp (for gcc) in 
#GPU=y #option A: (GPU+CPU) GPUs reduce workload from CPU by diving a subdomain  
#GPU_EXCLUSIVE=y #option B (preffered) : subdomains lie entirely on GPUs
#GPU_USE_INTRINSICS=y #for a more optimized gpu version
 
#----------------------------------------------------------------------------------------  
# settings for the compiler 
#----------------------------------------------------------------------------------------
cc= gcc#use "icc" for intel compilers 
FC= h5pfc
FLAGS= -O3 -cpp -fno-range-check -ffree-line-length-none #use to run the simulation using GCC compiler
#FLAGS= -O3 -cpp -ipo -xHost#for intel compiler 
#FLAGS= -O3 -g -fbacktrace -fcheck=all -cpp -fno-range-check -ffree-line-length-none -fbounds-check -fimplicit-none#-fimplicit-none # use for debugging 
#FLAGS= -O3 -Mcuda=7.5 # for GPU runs (using PGI fortran) 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAKEFILE SETTINGS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#--------------------------------------------------------------------- 
#    DO NOT CHANGE THE FOLLOWING UNLESS THE MAIN CODE IS MODIFIED  
#-----------------------------------------------------------------------------------------
ifeq ($(THREE_DIMS),y)
else
	CUSTOM :=$(CUSTOM) -DtwoD
endif

TYPE=PIC
CUSTOM :=$(CUSTOM) -Dmulflvr#to have muplitle species of different q/m 

ifeq ($(TYPE),PIC)
	 MOVDEP=movdepVEC.o#for standard PIC simualtion
	 CUSTOM :=$(CUSTOM) -DPIC
else ifeq ($(TYPE),GyroKinetic)
	MOVDEP=movdepGyroKinetic.o#for gyrokinetic PIC simulation   
	CUSTOM :=$(CUSTOM) -DGyroKinetic 
else ifeq ($(TYPE),Hybrid1)
	MOVDEP=movdepPIC.o
	 CUSTOM :=$(CUSTOM) -DHybrid1
else ifeq ($(TYPE),Hybrid2)
	MOVDEP=movdepPIC.o
	 CUSTOM :=$(CUSTOM) -DHybrid2 
else ifeq ($(TYPE),Hybrid3)
	MOVDEP=movdepPIC.o
	CUSTOM :=$(CUSTOM) -DHybrid3
endif



CPU=y
GPU_INCLUDE=n
ifeq ($(GPU_EXCLUSIVE),y)
	CPU=n
endif

#Add GPU flag for preprocessing 
ifeq ($(CPU),y)
	CUSTOM :=$(CUSTOM) -DCPU
endif
#Add GPU flag for preprocessing 
ifeq ($(GPU),y)
	CUSTOM :=$(CUSTOM) -Dgpu
	GPU_INCLUDE=y
endif
#add openmp flag for preprocessing code 
ifeq ($(OPEN_MP),y)
	CUSTOM :=$(CUSTOM) -DOPEN_MP
endif
#Add Exclusive GPU flag 
ifeq ($(GPU_EXCLUSIVE),y)
	CUSTOM :=$(CUSTOM) -DGPU_EXCLUSIVE
	GPU_INCLUDE=y
endif
#Use some more optimised GPU subroutines
ifeq ($(GPU_USE_INTRINSICS),y)
	CUSTOM :=$(CUSTOM) -DGPU_USE_INTRINSICS
endif


 
FFLAGS= $(CUSTOM) $(FLAGS)
EXECUTABLE= PICTOR

ifeq ($(GPU_INCLUDE),y)
OBJS= parameters.o vars.o var_gpu.o communication.o memory.o interpolation.o prtl_stats.o deposit.o $(MOVDEP) comm_fldprtl.o loadprtlout_gpu.o particles.o fields.o \
	 HDF5write.o savedata_routines.o comm_savedata.o comm_fld_gpu.o savedata_gpu.o savedata.o help_setup.o DerivedQnty.o comm_loadbalance.o loadbalance.o $(SETUP).o reload.o initialise.o \
	 comm_prtl_gpu.o inc_scan.o thrust_module.o memory_gpu.o fields_gpu.o movdep_gpu.o initialise_gpu.o MAIN.o
else
OBJS= parameters.o vars.o communication.o memory.o interpolation.o prtl_stats.o deposit.o $(MOVDEP) comm_fldprtl.o loadprtlout.o particles.o fields.o HDF5write.o\
	 savedata_routines.o comm_savedata.o savedata.o help_setup.o DerivedQnty.o comm_loadbalance.o loadbalance.o $(SETUP).o reload.o initialise.o MAIN.o
endif 

all: $(EXECUTABLE)

%.o:%.F90
	$(FC) $(FFLAGS) $(INCPATH) -c $^
%.o:%.cu
	nvcc -c $^	

$(EXECUTABLE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(EXECUTABLE)
