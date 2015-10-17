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

all: kMC_filmdissolve.o kMC_h2form.o kMC_vacdebonding.o
	$(FC) -o model.exe $(CFLAG) $(OFLAG) $(INC) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 calculate_sulfidation_rates.f90 solve_pf.f90 solve_mu.f90 solve_pH.f90 solve_ang.f90 solve_pot.f90 write_parameters.f90 thermo.f90 kMC_filmdissolve.o kMC_h2form.o kMC_vacdebonding.o allocate.f90 distrib_pf.f90 diffusivities.f90 ${PETSC_LIB} ${SPPARKS_LIB}

kMC_filmdissolve.o: kMC_filmdissolve.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 kMC_filmdissolve.f90

kMC_h2form.o: kMC_h2form.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 kMC_h2form.f90

kMC_vacdebonding.o: kMC_vacdebonding.f90
	$(FC) $(CFLAG) -c $(OFLAG) $(INC) shared_modules.f90 kMC_vacdebonding.f90

clean:
	rm -f *.exe *.out *.mod *.o *__genmod.f90

purge:
	rm -f *.exe *.out *.mod *.o *__genmod.f90 input.metfilm *.spparksscript

