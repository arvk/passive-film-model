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

  DM da
  KSP ksp_mu, ksp_pH, ksp_ang
  SNES snes_pf, snes_pot
  Vec state, onlymus
  integer :: iter  ! Current iteration number (Loop)
  integer :: ierr,status(MPI_STATUS_SIZE)  ! MPI variables

  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscViewer :: viewer

  integer :: startx,starty,startz,widthx,widthy,widthz
  integer :: x,y,z
  real*8 :: mysum

  type(context) simstate

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

  call DMDACreate3D(MPI_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,psx_g,psy_g,psz_g,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,nfields,1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,simstate%lattval,ierr)

  call DMDAGetCorners(simstate%lattval,simstate%startx,simstate%starty,simstate%startz,simstate%widthx,simstate%widthy,simstate%widthz,ierr)

  call KSPCreate(MPI_COMM_WORLD,ksp_mu,ierr)
  call KSPCreate(MPI_COMM_WORLD,ksp_pH,ierr)
  call KSPCreate(MPI_COMM_WORLD,ksp_ang,ierr)
  call SNESCreate(MPI_COMM_WORLD,snes_pf,ierr)
  call SNESCreate(MPI_COMM_WORLD,snes_pot,ierr)

  call KSPSetDM(ksp_mu,simstate%lattval,ierr)
  call KSPSetDM(ksp_pH,simstate%lattval,ierr)
  call KSPSetDM(ksp_ang,simstate%lattval,ierr)
  call SNESSetDM(snes_pot,simstate%lattval,ierr)

  call KSPSetFromOptions(ksp_mu,ierr)
  call KSPSetFromOptions(ksp_pH,ierr)
  call KSPSetFromOptions(ksp_ang,ierr)
  call SNESSetFromOptions(snes_pot,ierr)

  call allocate_matrices() ! Allocate all field matrices
  call thermo()            ! Calculate phase stabilities
  call diffusivities()     ! Calculate phase diffusivities
  call estimate_timestep() ! Estimate timestep for all field evolution equations
  call seed_prng()         ! Seed the Pseudo-random-number-generator

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

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,state,statepointer,ierr)

  do x = simstate%startx , simstate%startx + simstate%widthx-1
     do y = simstate%starty , simstate%starty + simstate%widthy-1
        do z = simstate%startz , simstate%startz + simstate%widthz-1
           statepointer(nmet,x,y,z) = met_g(x+1,y+1,z+1)
           statepointer(nmkw,x,y,z) = mkw_g(x+1,y+1,z+1)
           statepointer(npht,x,y,z) = pht_g(x+1,y+1,z+1)
           statepointer(npyr,x,y,z) = pyr_g(x+1,y+1,z+1)
           statepointer(nenv,x,y,z) = env_g(x+1,y+1,z+1)
           statepointer(nmus,x,y,z) = mu_g(x+1,y+1,z+1)
           statepointer(npH,x,y,z) = pH_in
           statepointer(nang,x,y,z) = opyr_g(x+1,y+1,z+1)
           statepointer(npot,x,y,z) = elpot_g(x+1,y+1,z+1)
        end do
     end do
  end do
  call DMDAVecRestoreArrayF90(simstate%lattval,state,statepointer,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)


  do iter = 1,nomc
!     write(6,*) 'In iteration',iter
!     call para_musolve(iter,ksp_mu)
!     call para_pHsolve(iter,ksp_pH)
     call para_pfsolve(iter,snes_pf,simstate)
!     call para_angsolve(iter,ksp_ang)
!     call para_potsolve(iter,snes_pot,state)


!  call VecStrideNorm(state,nmet,NORM_1,mysum,ierr)
!  write(6,*) 'SUM', mysum, rank

  end do


  call VecCreate(MPI_COMM_WORLD,onlymus,ierr)
  call VecSetSizes(onlymus,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(onlymus,ierr)

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecStrideGather(state,nmet,onlymus,INSERT_VALUES,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

!  call VecView(onlymus,PETSC_NULL_OBJECT,ierr)

!   do iter = 1,nomc  ! TIME LOOP

!      call pfsolve(iter)  ! Solve PF evolution equations
!      call musolve(iter)  ! Solve MU evolution equations
!      call orsolve(iter)  ! Solve OR evolution equations

!      if (include_dissolve) call dissolve_film() ! If dissolve tag is set, dissolve film
!      if (include_electro) call elpotsolve(iter) ! If electro tag is set, solve the electric potential field evolution
!      ! if (iter.eq.nomc/10) call voids_create() ! Create voids every nomc/10 steps


! !!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

!      if (mod(iter,freq_scale).eq.0) then  ! Start kMC routines every freq_scale steps

!         call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS
!         call spparks_filmenv()
!         call spparks_metfilm()
!      end if

! !!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

!      if (mod(iter,max(floor(real(nomc/noimg)),1)).eq.0) then   ! Write output files

!         call gather_pf()      ! Collect phase field to the parent process
!         call gather_mu()      ! Collect mu field to the parent process
!         call gather_opyr()    ! Collect orientation field to the parent process
!         call gather_electro() ! Collect electric potential field to the parent process

!         if (isroot) write(6,'(A,I5.5,A,I6.6,A,I6.6)') " INFO: Writing snapshot ",iter/max(floor(real(nomc/noimg)),1)," at iteration ",iter,"/",nomc
!         if (isroot) call write_fields(iter) ! IF I am ROOT, write out simulation parameters

!      end if


!   end do

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  !! Destroy Unused objects
  call KSPDestroy(ksp_mu,ierr)
  call KSPDestroy(ksp_pH,ierr)
  call KSPDestroy(ksp_ang,ierr)
  call SNESDestroy(snes_pf,ierr)
  call DMDestroy(da,ierr)


  !! Finalize Parallelization
  call PetscFinalize(ierr)
  call mpi_finalize(ierr)


end program passive_film_model
