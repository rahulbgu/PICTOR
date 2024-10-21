#MakeFile for the PIC simulation 
SHELL=/bin/sh
#---------------------------------------------------------------------------------------
#  Configure the PIC simulation
#---------------------------------------------------------------------------------------
SETUP=setup_ionize
SETUP=setup_test_pml


#THREE_DIMS=y# y = three spatial dimensions
#CORD=cyl# cordinate system (default is Cartesian; 'cyl' is cylinderical)
#OPEN_MP=y#to enable openmp support NOTE: include -qopenmp(intel) or -fopenmp (gcc) in "FLAGS"
#GPU=y#Use GPUs

# ----- Algorithm Selection  -----
#FDTD=Lehe# options :  Lehe ; ( the standard Yee scheme is the default)



 
#----------------------------------------------------------------------------------------  
# settings for the compiler 
#----------------------------------------------------------------------------------------
cc= gcc#use "icc" for intel compilers 
FC=h5pfc
FLAGS= -O3 -cpp -fno-range-check -ffree-line-length-none -fallow-argument-mismatch#use to run the simulation using GCC compiler
#FLAGS= -O3 -cpp -ipo -xHost#for intel compiler 
FLAGS= -O0 -g -finit-real=nan -Wmaybe-uninitialized -fbacktrace -fcheck=all -cpp -fno-range-check -ffree-line-length-none -fbounds-check -fimplicit-none -fallow-argument-mismatch#-fimplicit-none # use for debugging 
#FLAGS= -O3 -ta=tesla:cc70 -Mcuda# for GPU runs (using PGI fortran) 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF MAKEFILE SETTINGS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#---------------------------------------------------------------------------------------- 
#   DO NOT CHANGE ANYTHING BELOW THIS UNLESS THE STRUCTURE OF SOURCE CODE IS MODIFIED  
#-----------------------------------------------------------------------------------------
ifeq ($(THREE_DIMS),y)
else
	CUSTOM :=$(CUSTOM) -DtwoD
endif

ifeq ($(CORD),cyl)
	CUSTOM :=$(CUSTOM) -Dcyl
endif

CPU=y#Default option: everything on CPU

#Switch to GPU
ifeq ($(GPU),y)
	CUSTOM :=$(CUSTOM) -Dgpu
	CPU=n
endif
#Add CPU flag for preprocessing 
ifeq ($(CPU),y)
	CUSTOM :=$(CUSTOM) -DCPU
endif
#add openmp flag for preprocessing code 
ifeq ($(OPEN_MP),y)
	CUSTOM :=$(CUSTOM) -DOPEN_MP
endif
#Use some more optimised GPU subroutines
ifeq ($(GPU_USE_INTRINSICS),y)
	CUSTOM :=$(CUSTOM) -DGPU_USE_INTRINSICS
endif



# -------- Algorithms -----------
FDTD_SRC_EXT =# default
ifeq ($(FDTD),Lehe)
#	FDTD_SRC = em_update_lehe
FDTD_SRC_EXT =_lehe# default
endif





 
FFLAGS= $(CUSTOM) $(FLAGS)
EXECUTABLE= PICTOR
VPATH=src:src/cyl:src/bc:src/io:src/maxwell:src/movdep:src/mem:src/init:src/comm:src/comm/subdomains:src/ionization
OBJDIR=build

ifeq ($(GPU),y)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o var_gpu.o communication.o comm_fld.o subdomains.o prtl_tag.o prob_dist.o mem_prtl.o mem_fld.o interpolation.o prtl_stats.o reload.o initialise.o deposit.o deposit_gpu.o loadprtlout.o movdep.o usc_CFDTD_gpu.o usc_CFDTD.o bc_curved_surf.o comm_prtl.o em_update.o fields.o \
	 HDF5write.o particles_gpu.o savedata_gpu.o savedata_routines.o  comm_fld_gpu.o savedata.o help_setup.o DerivedQnty.o \
	 interp_gpu.o em_update_gpu.o fields_gpu.o movdep_gpu.o initialise_gpu.o bc_exterior_field.o bc_open.o bc_attenuate.o bc_inflow.o bc_domain_size.o loadbalance_gpu.o loadbalance_indep.o bc_gpu.o bc.o $(SETUP).o MAIN.o)

ifeq ($(CORD),cyl)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o var_gpu.o cyl_vars.o communication.o subdomains.o prob_dist.o prtl_tag.o mem_prtl.o mem_fld.o interpolation.o prtl_stats.o cyl_em_update.o deposit.o deposit_gpu.o movdep.o cyl_comm_fld.o comm_fld.o cyl_bc_gpu.o cyl_bc.o loadprtlout_gpu.o comm_prtl.o cyl_filter.o fields.o\
	 HDF5write.o cyl_savedata_routines.o savedata_gpu.o savedata_routines.o  cyl_comm_fld_gpu.o comm_fld_gpu.o savedata.o help_setup.o cyl_help_setup.o DerivedQnty.o reload.o initialise.o \
	 particles_gpu.o interp_gpu.o cyl_em_update_gpu.o cyl_filter_gpu.o fields_gpu.o cyl_movdep_gpu.o initialise_gpu.o cyl_initialise.o bc_domain_size.o loadbalance_gpu.o loadbalance_indep.o cyl_loadbalance.o $(SETUP).o MAIN.o)
endif


else


OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o mem_fld.o mem_prtl.o communication.o prtl_tag.o prob_dist.o comm_fld.o loadprtlout.o comm_prtl.o subdomains.o interpolation.o prtl_stats.o reload.o init_prtl.o initialise.o movdep_routines.o deposit.o movdep.o bc_exterior_field.o bc_reflect.o bc_open.o bc_adjust_speed.o bc_attenuate.o bc_pml.o bc_inflow.o usc_CFDTD.o bc_curved_surf.o bc.o em_update$(FDTD_SRC_EXT).o current.o HDF5write.o\
	 read_data.o savedata_routines.o  savedata.o help_setup.o DerivedQnty.o bc_domain_size.o loadbalance_indep.o loadbalance_procgrid.o $(SETUP).o MAIN.o)
ifeq ($(CORD),cyl)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o mem_fld.o mem_prtl.o cyl_common.o communication.o prtl_tag.o prob_dist.o comm_fld.o loadprtlout.o comm_prtl.o subdomains.o interpolation.o prtl_stats.o reload.o init_prtl.o initialise.o movdep_routines.o deposit.o cyl_movdep_strang.o bc_exterior_field.o bc_reflect.o bc_open.o bc_adjust_speed.o bc_attenuate.o bc_inflow.o usc_CFDTD.o bc_curved_surf.o bc.o cyl_em_update.o current.o HDF5write.o\
	 read_data.o cyl_savedata_routines.o savedata_routines.o  savedata.o help_setup.o DerivedQnty.o bc_domain_size.o loadbalance_indep.o loadbalance_procgrid.o $(SETUP).o MAIN.o)
endif

endif 

all: $(EXECUTABLE)

$(OBJDIR)/%.o:%.F90
	$(FC) $(FFLAGS) $(INCPATH) -c $^ -o $@
$(OBJDIR)/%.o:%.cu
	nvcc -c $^ -o $@

$(EXECUTABLE): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS)
	mv *.mod $(OBJDIR)/

clean: 
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/*.mod
	rm -f $(EXECUTABLE)
