# Define fortran compiler and linker
FC=mpif90

# Define compiler flags
CFLAG= -CB

# Define optimization flags
OFLAG= -O3 -xhost -ipo

all: 
	$(FC) $(CFLAG) $(OFLAG) shared_modules.f90 wrapper.f90 read_parameters.f90 distrib_params.f90 seed_prng.f90 initialize_geometry.f90 swap_pf.f90 pfsolve.f90 musolve.f90 orsolve.f90 iccg.f90 write_parameters.f90 read_geometry.f90 gather_pf.f90 gather_mu.f90 gather_opyr.f90 write_fields.f90 calc_grad_pf.f90 calc_lap_pf.f90 calc_lap_mu.f90 thermo.f90 initialize_kmc.f90 swap_kmc.f90 swap_mu.f90 swap_or.f90 kmcsolve.f90 gather_kmc.f90 allocate.f90 couple_kmc_pf.f90 distrib_kmc.f90 distrib_pf.f90 diffusivities.f90 estimate_timestep.f90 -o model.exe

clean: 
	rm -f *.exe *.out *.mod kgrid pgrid

