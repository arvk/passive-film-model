include ${PETSC_DIR}/conf/variables

# Define fortran compiler and linker
FC=mpif90

# Define compiler flags
CFLAG= -cpp -traceback -zero -unroll0

# Define optimization flags
OFLAG= -O2 -xHost -ipo

# Include files
INC= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

# SPPARKS library
SPPARKS_LIB= -L${SPPARKS_LIB_DIR} -lspparks

all: spparks_filmenv.o spparks_metfilm.o spparks_vacmet.o
	$(FC) -o model.exe $(CFLAG) $(OFLAG) $(INC) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 calculate_sulfidation_rates.f90 solve_pf.f90 solve_mu.f90 solve_pH.f90 solve_ang.f90 solve_pot.f90 write_parameters.f90 thermo.f90 spparks_filmenv.o spparks_metfilm.o spparks_vacmet.f90 allocate.f90 distrib_pf.f90 diffusivities.f90 estimate_timestep.f90 ${PETSC_LIB} ${SPPARKS_LIB}

spparks_filmenv.o: spparks_filmenv.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_filmenv.f90

spparks_metfilm.o: spparks_metfilm.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_metfilm.f90

spparks_vacmet.o: spparks_vacmet.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 spparks_vacmet.f90

clean:
	rm -f *.exe *.out *.mod *.o *__genmod.f90

