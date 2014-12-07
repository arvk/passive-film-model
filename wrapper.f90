program passive_film_model
  use commondata
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: iter      ! Loop variable for current iteration 
  integer :: rank_loop
  integer :: deallocatestatus
  integer :: ierr,status(MPI_STATUS_SIZE)   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize Parallelization
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,procs,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)

  !! Define parent process (rank 0)
  if(rank .eq. 0)then
     isparent = .TRUE.
  else
     isparent = .FALSE.
  end if

  !! Read in input parameters -- Parent Processor only
  if (isparent) then
     call read_parameters()
  end if

  !! Distribute input parameters -- MPI Broadcast by root
  call distrib_params()

  !! Import the thermodynamic module -- All processors
  call thermo()
  call diffusivities()

  !! Allocate all matrices -- All processors
  call allocate_matrices()

  !! Construct/read the system -- All processors
  if (isparent) then
     if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then
        call read_geometry()
     else 
        call system("rm -rf *.out")   ! Cleanup before running
        call initialize_geometry()
     end if
  end if

  call distrib_pf()

  if(isparent) then
     call write_parameters()
  end if

  !! Equilibration of chemical potentials
  !#####! call pre_equilibrate()

  call mpi_barrier(MPI_COMM_WORLD,ierr)

  do iter = 1,nomc

     !! Solve PF equations
     call pfsolve(iter)

     ! if (iter.lt. (nomc/100)) then
     !    call pre_equilibrate()
     ! end if

!        read(*,*)


     if (mod(iter,freq_scale).eq.0) then

        if (isparent) then
           call initialize_kmc()
        end if

        call distrib_kmc()

        call kmcsolve(iter)

!!!! Couple KMC surface to PF

        call gather_pf()
        call gather_kmc()

        if (isparent) then           
           call couple_kmc_pf()
        end if
        call distrib_pf()

     end if


     ! !!         Write output files at fixed times
     if (mod(iter,(nomc/noimg)-1).eq.0) then

        call gather_pf()

        if(isparent)then
           call write_fields(iter)
        end if

     end if


     !! Draw progress bar
     !     call draw_progress_bar(iter,nomc,dt)

     !  write(6,*)   ! Write out empty line after the progress bar

  end do

  call mpi_finalize(ierr)


end program passive_film_model
