include $(PETSC_DIR)/conf/variables

CPPFLAGS = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
OPTIM_FLAGS = -O2 -xHost -ipo -cpp 
AM_FCFLAGS = $(OPTIM_FLAGS) -cpp -traceback -zero -unroll0 -dynamic $(CPPFLAGS)

bin_PROGRAMS = pfm
pfm_SOURCES = shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 initialize_geometry.f90 calculate_sulfidation_rates.f90 solve_pf.f90 solve_mu.f90 solve_pH.f90 solve_ang.f90 solve_pot.f90 write_parameters.f90 thermo.f90 kMC_filmdissolve.f90 kMC_h2form.f90 kMC_vacdebonding.f90 allocate.f90 diffusivities.f90
pfm_LDFLAGS = -Wl,-rpath,$(SPPARKS_LIB_DIR) -L$(SPPARKS_LIB_DIR) -lspparks
pfm_LDADD = $(PETSC_LIB) -L$(SPPARKS_LIB_DIR) -lspparks

purge:
	make clean ; rm -f *.exe *.out *.mod *.o *__genmod.f90 input.metfilm *.spparksscript *.vts
