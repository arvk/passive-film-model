program passive_film_model
  use commondata
  use fields
  use kmc_data
  implicit none
#include <finclude/petscsys.h>

  integer :: iter  ! Current iteration number (Loop)
  integer :: ierr,status(MPI_STATUS_SIZE)  ! MPI variables

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  !! Initialize Parallelization
  call mpi_init(ierr)
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,procs,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)


  if (rank .eq. 0) then ! If I am the parent process (rank == 0)
     isroot = .TRUE.
  else                  ! If I am NOT the parent process (rank =/= 0)
     isroot = .FALSE.
  end if


  if (isroot) call read_parameters()  ! IF I am ROOT, read in input parameters
  if (isroot) call write_parameters() ! IF I am ROOT, write out simulation parameters

  call distrib_params()    ! MPI-Distribute input parameters to non-parent processors
  call thermo()            ! Calculate phase stabilities
  call diffusivities()     ! Calculate phase diffusivities
  call estimate_timestep() ! Estimate timestep for all field evolution equations
  call seed_prng()         ! Seed the Pseudo-random-number-generator
  call allocate_matrices() ! Allocate all field matrices

  if (isroot) then
     if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then   ! If it is a restarted run, read PF matrices
        call read_geometry()
     else                                                       ! If it is a fresh calculation
        call initialize_geometry()
     end if
  end if

  call distrib_pf()    ! Distribute all PF-MU-OR-ELPOT matrices to non-parent processors

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning time loop

  do iter = 1,nomc  ! TIME LOOP

     call pfsolve(iter)  ! Solve PF evolution equations
     call musolve(iter)  ! Solve MU evolution equations
     call orsolve(iter)  ! Solve OR evolution equations

     if (include_dissolve) call dissolve_film() ! If dissolve tag is set, dissolve film
     if (include_electro) call elpotsolve(iter) ! If electro tag is set, solve the electric potential field evolution
     ! if (iter.eq.nomc/10) call voids_create() ! Create voids every nomc/10 steps


!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

     if (mod(iter,freq_scale).eq.0) then  ! Start kMC routines every freq_scale steps

        call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS
        call spparks_filmenv()

     end if

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

     if (mod(iter,max(floor(real(nomc/noimg)),1)).eq.0) then   ! Write output files

        call gather_pf()      ! Collect phase field to the parent process
        call gather_mu()      ! Collect mu field to the parent process
        call gather_opyr()    ! Collect orientation field to the parent process
        call gather_electro() ! Collect electric potential field to the parent process

        if (isroot) write(6,'(A,I5.5,A,I6.6,A,I6.6)') " INFO: Writing snapshot ",iter/max(floor(real(nomc/noimg)),1)," at iteration ",iter,"/",nomc
        if (isroot) call write_fields(iter) ! IF I am ROOT, write out simulation parameters

     end if


  end do

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  !! Finalize Parallelization
  call PetscFinalize(ierr)
  call mpi_finalize(ierr)


end program passive_film_model
