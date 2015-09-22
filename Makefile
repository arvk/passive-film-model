include ${PETSC_DIR}/conf/variables

# Define fortran compiler and linker
FC=mpif90

# Define compiler flags
CFLAG= -cpp -traceback -gen-interfaces -warn interfaces -u -fpe3

# Define optimization flags
OFLAG= -O3 -xHost -ipo

# Include files
INC= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

# SPPARKS library
SPPARKS_LIB= -L${SPPARKS_LIB_DIR} -lspparks

all: spparks_filmenv.o spparks_metfilm.o spparks_vacmet.o
	$(FC) -o model.exe $(CFLAG) $(OFLAG) $(INC) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 calculate_sulfidation_rates.f90 pll_pfsolve.f90 pll_musolve.f90 pll_pHsolve.f90 pll_angsolve.f90 pll_potsolve.f90 write_parameters.f90 read_geometry.f90 thermo.f90 spparks_filmenv.o spparks_metfilm.o spparks_vacmet.f90 allocate.f90 distrib_pf.f90 diffusivities.f90 estimate_timestep.f90 ${PETSC_LIB} ${SPPARKS_LIB}

spparks_filmenv.o: spparks_filmenv.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_filmenv.f90

spparks_metfilm.o: spparks_metfilm.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_metfilm.f90

spparks_vacmet.o: spparks_vacmet.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_vacmet.f90

clean:
	rm -f *.exe *.out *.mod *.o *__genmod.f90

