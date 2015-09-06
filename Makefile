include ${PETSC_DIR}/conf/variables

# Define fortran compiler and linker
FC=mpif90

# Define compiler flags
CFLAG= -CB -fp-model precise -cpp

# Define optimization flags
OFLAG= -O0 -xhost

# Include files
INC= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

# SPPARKS library
SPPARKS_LIB= -L${SPPARKS_LIB_DIR} -lspparks

all: spparks_filmenv.o spparks_metfilm.o
	$(FC) -o model.exe $(CFLAG) $(OFLAG) $(INC) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 swap_pf.f90 pfsolve.f90 pll_pfsolve.f90 nonlin_feval.f90 nonlin_jacob.f90 pll_musolve.f90 musolve.f90 pll_pHsolve.f90 pll_angsolve.f90 orsolve.f90 elpotsolve.f90 electro_feval.f90 electro_jacob.f90 swap_electro.f90 dissolve_film.f90 voids_create.f90 write_parameters.f90 read_geometry.f90 gather_pf.f90 gather_mu.f90 gather_opyr.f90 gather_electro.f90 write_fields.f90 calc_grad_pf.f90 thermo.f90 spparks_filmenv.o spparks_metfilm.o swap_mu.f90 swap_or.f90 allocate.f90 distrib_pf.f90 diffusivities.f90 estimate_timestep.f90 ${PETSC_LIB} ${SPPARKS_LIB}

spparks_filmenv.o: spparks_filmenv.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_filmenv.f90

spparks_metfilm.o: spparks_metfilm.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_metfilm.f90

clean:
	rm -f *.exe *.out *.mod *.o kgrid pgrid

