program passive_film_model
  use commondata
  use fields
  use kmc_data
  implicit none
#include <finclude/petscsys.h>


  integer :: x, y, z   ! Loop variables
  integer :: iter      ! Loop variable for current iteration 
  integer :: rank_loop
  integer :: deallocatestatus
  integer :: ierr,status(MPI_STATUS_SIZE)   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize Parallelization
  call mpi_init(ierr)
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

  call mpi_comm_size(MPI_COMM_WORLD,procs,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

  !! Define Root process (rank 0)
  if(rank .eq. 0)then
     isroot = .TRUE.
  else 
     isroot = .FALSE.
  end if


  !=======================
  ! SYSTEM INITIALIZATION
  !=======================


  !! Read in input parameters -- Root only
  if (isroot) then
     call read_parameters()
  end if

  !! Distribute input parameters -- MPI Broadcast by Root
  call distrib_params()

  !! Import the thermodynamic module -- All processors
  call thermo()
  call diffusivities()
  call estimate_timestep()
  call seed_prng()

  !! Allocate all matrices -- All processors
  call allocate_matrices()

  !! Construct/read the system -- All processors
  if (isroot) then
     if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then
        call read_geometry()
     else 
        call system("rm -rf *.out")   ! Cleanup before running
        call initialize_geometry()
     end if
  end if

  call distrib_pf()

  if(isroot) then
     call write_parameters()
  end if

  !! Barrier before beginning time loop
  call mpi_barrier(MPI_COMM_WORLD,ierr)

  ! =========
  ! TIME LOOP
  ! =========

  do iter = 1,nomc

     !! Solve PF equations
     call pfsolve(iter)
     call musolve(iter)
     call orsolve(iter)

     if (include_dissolve) then
        call dissolve_film()
     end if


     if (iter.eq.nomc/10) then
        call voids_create()
     end if
  
  
     if (mod(iter,freq_scale).eq.0) then

        if (isroot) then
           call initialize_kmc()
        end if

        call distrib_kmc()

        call kmcsolve(iter)

!!!! Couple KMC surface to PF

        call gather_pf()
        call gather_mu()
        call gather_opyr()
        call gather_kmc()

        if (isroot) then           
           call couple_kmc_pf()
        end if
        call distrib_pf()

     end if


!!!         Write output files at fixed times
     if (mod(iter,(nomc/noimg)-1).eq.0) then

        call gather_pf()
        call gather_mu()
        call gather_opyr()

        if(isroot)then
           call write_fields(iter)
        end if

     end if

  end do

  call PetscFinalize(ierr)
  call mpi_finalize(ierr)


end program passive_film_model
