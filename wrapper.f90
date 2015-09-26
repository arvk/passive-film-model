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

  !!#This program simulates the growth and breakdown of iron sulfide films formed in sour corrosive conditions using a combined phase-field and kinetic Monte Carlo algorithm.
  !---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  KSP ksp_mu, ksp_pH, ksp_ang
  SNES snes_pf, snes_pot
  PetscScalar, pointer :: statepointer(:,:,:,:)
  integer :: iter  !! Current iteration number in the time-stepping loop
  integer :: ierr,status(MPI_STATUS_SIZE)  !! MPI error and status variables
  integer :: x,y,z  !! Coordinates inside the simulation cell
  type(context) simstate !! Field variables stored in PETSc vectors and DMDA objects

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
  call estimate_timestep() ! Estimate timestep for all field evolution equations
  call seed_prng()         ! Seed the Pseudo-random-number-generator
  call thermo()            ! Calculate phase stabilities

  if (isroot) then
     if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then   ! If it is a restarted run, read PF matrices
        call read_geometry()
     else                                                       ! If it is a fresh calculation
        call initialize_geometry()
     end if
  end if

  call calculate_sulfidation_rates()    ! Calculate sulfidation rates on different FeS phases

  call distrib_pf()    ! Distribute all PF-MU-OR-ELPOT matrices to non-parent processors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  do x = simstate%startx , simstate%startx + simstate%widthx-1
     do y = simstate%starty , simstate%starty + simstate%widthy-1
        do z = simstate%startz , simstate%startz + simstate%widthz-1
           statepointer(nmet,x,y,z) = met_g(x+1,y+1,z+1)
           statepointer(nmkw,x,y,z) = mkw_g(x+1,y+1,z+1)
           statepointer(npht,x,y,z) = pht_g(x+1,y+1,z+1)
           statepointer(npyr,x,y,z) = pyr_g(x+1,y+1,z+1)
           statepointer(nenv,x,y,z) = env_g(x+1,y+1,z+1)
           statepointer(nmus,x,y,z) = mu_g(x+1,y+1,z+1)
           statepointer(npH,x,y,z) = 10**(0.0d0-pH_in)
           statepointer(nang,x,y,z) = opyr_g(x+1,y+1,z+1)
           statepointer(npot,x,y,z) = elpot_g(x+1,y+1,z+1)
           statepointer(nvoi,x,y,z) = 1.0d0
        end do
     end do
  end do
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning time loop

  do iter = 1,nomc

     call VecCopy(simstate%slice,simstate%exslice,ierr)

     write(6,*) 'In iteration',iter
     call SNESCreate(MPI_COMM_WORLD,snes_pf,ierr); call para_pfsolve(iter,snes_pf,simstate); call SNESDestroy(snes_pf,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_mu,ierr); call para_musolve(iter,ksp_mu,simstate); call KSPDestroy(ksp_mu,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_pH,ierr); call para_pHsolve(iter,ksp_pH,simstate); call KSPDestroy(ksp_pH,ierr)
     call KSPCreate(MPI_COMM_WORLD,ksp_ang,ierr); call para_angsolve(iter,ksp_ang,simstate); call KSPDestroy(ksp_ang,ierr)
     call SNESCreate(MPI_COMM_WORLD,snes_pot,ierr); call para_potsolve(iter,snes_pot,simstate); call SNESDestroy(snes_pot,ierr)

     call spparks_filmenv(iter,simstate)
     call spparks_metfilm(iter,simstate)
     call spparks_vacmet(iter,simstate)

  end do

  call VecView(simstate%slice,PETSC_NULL_OBJECT,ierr)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Destroy Unused objects
  call VecDestroy(simstate%slice,ierr)
  call VecDestroy(simstate%exslice,ierr)
  call DMDestroy(simstate%lattval,ierr)

  ! Finalize Parallelization
  call PetscFinalize(ierr)
end program passive_film_model
