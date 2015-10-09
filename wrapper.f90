!================================================================================
! Passive Film Model
! Aravind Krishnamoorthy, aravindk@mit.edu, Massachusetts Institute of Technology
!
! Copyright 2015 Aravind Krishnamoorthy. This software is distributed under
! the GNU General Public License. See doc/LICENSE.GPL for the license text.
!================================================================================
program passive_film_model
  use commondata
  use fields
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  !!##This program simulates the growth and breakdown of iron sulfide films formed in sour corrosive conditions using a combined phase-field and kinetic Monte Carlo algorithm.
  !----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  KSP ksp_mu     !! Linear chemical potential field solver
  KSP ksp_pH     !! Linear pH field solver
  KSP ksp_ang    !! Linear pyrite crystal shape solver
  SNES snes_pf   !! Non-linear phase-field solver
  SNES snes_pot  !! Non-linear electrical potential field solver
  PetscInt :: iter                               !! Current iteration number in the time-stepping loop
  PetscErrorCode ierr !! MPI error flag
  PetscInt :: status(MPI_STATUS_SIZE)       !! MPI status variables
  PetscInt :: x,y,z                              !! Coordinates inside the simulation cell
  PetscReal :: phase_volume_in_simcell(0:(nphases-1))  !! Number of gridpoints in the simulation cell occupied by a given phase
  PetscScalar :: metal_content_in_simcell     !! Amount of metal phase in the simulation cell
  PetscInt :: fesphase                        !! Index to identify the type of FeS phase
  PetscViewer :: snapshot_writer              !! Pointer to write-out simcell snapshot to file
  character*5 :: image_ID                       !! Index of current image (derived from current iteration number)
  type(context) simstate                        !! Field variables stored in PETSc vectors and DMDA objects

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Initialize Parallelization
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

  call DMDACreate3D(MPI_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
       & DMDA_STENCIL_STAR,psx_g,psy_g,psz_g,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,nfields, &
       & 1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,simstate%lattval,ierr)

  call DMDAGetCorners(simstate%lattval,simstate%startx,simstate%starty,simstate%startz, &
       & simstate%widthx,simstate%widthy,simstate%widthz,ierr)

  call DMCreateGlobalVector(simstate%lattval, simstate%slice, ierr)
  call DMCreateGlobalVector(simstate%lattval, simstate%exslice, ierr)



  call allocate_matrices() ! Allocate all field matrices
  call diffusivities()     ! Calculate phase diffusivities
  call seed_prng()         ! Seed the Pseudo-random-number-generator
  call thermo()            ! Calculate phase stabilities
  call estimate_timestep() ! Estimate timestep for all field evolution equations

  call initialize_geometry(simstate)

  call calculate_sulfidation_rates()    ! Calculate sulfidation rates on different FeS phases

  call VecStrideNorm(simstate%slice,nmet,NORM_1,metal_content_in_simcell,ierr)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning time loop

  do iter = 1,nomc

     call VecCopy(simstate%slice,simstate%exslice,ierr)

     call SNESCreate(MPI_COMM_WORLD,snes_pf,ierr); call solve_pf(iter,snes_pf,simstate); call SNESDestroy(snes_pf,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_mu,ierr); call solve_mu(iter,ksp_mu,simstate); call KSPDestroy(ksp_mu,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_pH,ierr); call solve_pH(iter,ksp_pH,simstate); call KSPDestroy(ksp_pH,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_ang,ierr); call solve_ang(iter,ksp_ang,simstate); call KSPDestroy(ksp_ang,ierr)
     call SNESCreate(MPI_COMM_WORLD,snes_pot,ierr); call solve_pot(iter,snes_pot,simstate); call SNESDestroy(snes_pot,ierr)

     if (mod(iter,kmc_freq).eq.0) call kMC_vacdebonding(iter,simstate,metal_content_in_simcell)
     if (mod(iter,kmc_freq).eq.0) call kMC_h2form(iter,simstate)
     if (mod(iter,kmc_freq).eq.0) call kMC_filmdissolve(iter,simstate,metal_content_in_simcell)

     if (mod(iter,stat_freq).eq.0) then
        do fesphase = 0,(nphases-1)
           call VecStrideNorm(simstate%slice,fesphase,NORM_1,phase_volume_in_simcell(fesphase),ierr)
        end do
           if (isroot) write(6,'(A,F9.0,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2)') " INFO: TIME= ", iter*dt, " s. MET= ", phase_volume_in_simcell(nmet), " MKW= ", phase_volume_in_simcell(nmkw), " PHT= ", phase_volume_in_simcell(npht), " PYR= ", phase_volume_in_simcell(npyr), " ENV= ", phase_volume_in_simcell(nenv)
        end if

     if (mod(iter,(nomc/num_images)).eq.0) then
        write(image_ID,'(I5.5)') iter/max(floor(real(nomc/num_images)),1)
        call PetscViewerASCIIOpen(MPI_COMM_WORLD,"SIMCELL_"//image_ID//".out",snapshot_writer,ierr)
        call VecView(simstate%slice,snapshot_writer,ierr)
        call PetscViewerDestroy(snapshot_writer,ierr)
     end if

  end do


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Destroy Unused objects
  call VecDestroy(simstate%slice,ierr)
  call VecDestroy(simstate%exslice,ierr)
  call DMDestroy(simstate%lattval,ierr)

  ! Finalize Parallelization
  call PetscFinalize(ierr)
end program passive_film_model
