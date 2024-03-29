#MakeFile for the PIC simulation 
SHELL=/bin/sh
#---------------------------------------------------------------------------------------
#  Configure the PIC simulation
#---------------------------------------------------------------------------------------

SETUP=setup_alfventurb
SETUP=setup_shock
#SETUP=setup_shear_flow
SETUP=setup_Bell
#SETUP=setup_cond_cyl_test
SETUP=setup_bennett_pinch_cart
#SETUP=setup_weibel_rel
#SETUP=setup_harris_sheet
#SETUP=setup_reconnection
SETUP=setup_PUI_shock
SETUP=setup_weibel


#THREE_DIMS=y# y = three spatial dimensions
#CORD=cyl# cordinate system (default is Cartesian; 'cyl' is cylinderical)
#OPEN_MP=y#to enable openmp support NOTE: include -qopenmp(intel) or -fopenmp (gcc) in "FLAGS"
#GPU=y#Use GPUs  
 
#----------------------------------------------------------------------------------------  
# settings for the compiler 
#----------------------------------------------------------------------------------------
cc= gcc#use "icc" for intel compilers 
FC= h5pfc
FLAGS= -O3 -cpp -fno-range-check -ffree-line-length-none -fallow-argument-mismatch#use to run the simulation using GCC compiler
#FLAGS= -O3 -cpp -ipo -xHost#for intel compiler 
#FLAGS= -O3 -g -fbacktrace -fcheck=all -cpp -fno-range-check -ffree-line-length-none -fbounds-check -fimplicit-none -fallow-argument-mismatch#-fimplicit-none # use for debugging 
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


 
FFLAGS= $(CUSTOM) $(FLAGS)
EXECUTABLE= PICTOR
VPATH=src:src/cyl
OBJDIR=build

ifeq ($(GPU),y)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o var_gpu.o communication.o prtl_tag.o prob_dist.o memory.o interpolation.o prtl_stats.o reload.o initialise.o deposit.o deposit_gpu.o movdep.o comm_fldprtl.o usc_CFDTD_gpu.o usc_CFDTD.o bc_curved_surf.o loadprtlout_gpu.o particles.o em_update.o fields.o \
	 HDF5write.o savedata_gpu.o savedata_routines.o comm_savedata.o comm_fld_gpu.o savedata.o help_setup.o DerivedQnty.o comm_prtl_gpu.o particles_gpu.o \
	 interp_gpu.o em_update_gpu.o fields_gpu.o movdep_gpu.o initialise_gpu.o comm_loadbalance.o loadbalance_gpu.o loadbalance.o injector_gpu.o injector.o bc_gpu.o bc.o $(SETUP).o MAIN.o)

ifeq ($(CORD),cyl)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o var_gpu.o cyl_vars.o communication.o prob_dist.o prtl_tag.o memory.o interpolation.o prtl_stats.o cyl_em_update.o deposit.o deposit_gpu.o movdep.o cyl_comm_fldprtl.o comm_fldprtl.o cyl_bc_gpu.o cyl_bc.o loadprtlout_gpu.o particles.o cyl_filter.o fields.o\
	 HDF5write.o cyl_savedata_routines.o savedata_gpu.o savedata_routines.o comm_savedata.o cyl_comm_fld_gpu.o comm_fld_gpu.o savedata.o help_setup.o cyl_help_setup.o DerivedQnty.o reload.o initialise.o \
	 comm_prtl_gpu.o particles_gpu.o interp_gpu.o cyl_em_update_gpu.o cyl_filter_gpu.o fields_gpu.o cyl_movdep_gpu.o initialise_gpu.o cyl_initialise.o comm_loadbalance.o loadbalance_gpu.o loadbalance.o cyl_loadbalance.o $(SETUP).o MAIN.o)
endif


else


OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o communication.o prtl_tag.o prob_dist.o memory.o interpolation.o prtl_stats.o reload.o initialise.o deposit.o movdep.o comm_fldprtl.o injector.o usc_CFDTD.o bc_curved_surf.o bc.o loadprtlout.o particles.o em_update.o fields.o HDF5write.o\
	 savedata_routines.o comm_savedata.o savedata.o help_setup.o DerivedQnty.o comm_loadbalance.o loadbalance.o $(SETUP).o MAIN.o)
ifeq ($(CORD),cyl)
OBJS= $(addprefix $(OBJDIR)/, parameters.o vars.o cyl_vars.o communication.o prtl_tag.o prob_dist.o memory.o interpolation.o prtl_stats.o cyl_em_update.o deposit.o cyl_comm_fldprtl.o comm_fldprtl.o cyl_bc_cpu.o cyl_bc.o cyl_movdep.o loadprtlout.o particles.o \
      cyl_filter.o fields.o HDF5write.o cyl_savedata_routines.o savedata_routines.o comm_savedata.o savedata.o help_setup.o cyl_help_setup.o DerivedQnty.o comm_loadbalance.o loadbalance.o cyl_loadbalance.o $(SETUP).o reload.o initialise.o cyl_initialise.o MAIN.o)
endif

endif 

all: $(EXECUTABLE)

$(OBJDIR)/%.o:%.F90
	$(FC) $(FFLAGS) $(INCPATH) -c $^ -o $@
$(OBJDIR)/%.o:%.cu
	nvcc -c $^ -o $@

$(EXECUTABLE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)
	mv *.mod $(OBJDIR)/

clean: 
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/*.mod
	rm -f $(EXECUTABLE)
