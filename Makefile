include ${PETSC_DIR}/conf/variables

# Define fortran compiler and linker
FC=mpif90

# Define compiler flags
CFLAG= -CB -fp-model precise -cpp -dM

# Define optimization flags
OFLAG= -O0 #3 -xhost -ipo

# Include files
INC= -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include

all: 
	$(FC) -o model.exe $(CFLAG) $(OFLAG) $(INC) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 swap_pf.f90 pfsolve.f90 musolve.f90 orsolve.f90 write_parameters.f90 read_geometry.f90 gather_pf.f90 gather_mu.f90 gather_opyr.f90 write_fields.f90 calc_grad_pf.f90 calc_lap_pf.f90 calc_lap_mu.f90 thermo.f90 initialize_kmc.f90 swap_kmc.f90 swap_mu.f90 swap_or.f90 kmcsolve.f90 gather_kmc.f90 allocate.f90 couple_kmc_pf.f90 distrib_kmc.f90 distrib_pf.f90 diffusivities.f90 estimate_timestep.f90 feval_env_met.f90  feval_env_pyr.f90  feval_met_pht.f90  feval_pht_env.f90  feval_pht_pyr.f90  feval_pyr_met.f90 feval_env_pht.f90  feval_met_env.f90  feval_met_pyr.f90  feval_pht_met.f90  feval_pyr_env.f90  feval_pyr_pht.f90 jacob_env_met.f90  jacob_env_pyr.f90  jacob_met_pht.f90  jacob_pht_env.f90  jacob_pht_pyr.f90  jacob_pyr_met.f90 jacob_env_pht.f90  jacob_met_env.f90  jacob_met_pyr.f90  jacob_pht_met.f90  jacob_pyr_env.f90  jacob_pyr_pht.f90  ${PETSC_LIB}

clean: 
	rm -f *.exe *.out *.mod kgrid pgrid

