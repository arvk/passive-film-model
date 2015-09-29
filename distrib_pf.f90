subroutine distrib_pf()
  use commondata
  use fields
  implicit none
#include <finclude/petscsys.h>
  !! **Distribute field variables -- different \(FeS\) phase fractions, chemical potential etc. to child processors**

  PetscInt :: x, y, z                      !! Coordinates inside the simulation system
  PetscInt :: ierr,status(MPI_STATUS_SIZE) !! MPI error and status flags

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Broadcast global phase-fraction fields
  call mpi_bcast(met_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(mkw_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pht_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pyr_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(env_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! Broadcast global chemical-potential fields
  call mpi_bcast(mu_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! Broadcast global electrical potential
  call mpi_bcast(elpot_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! Broadcast global pyrite orientation field
  call mpi_bcast(opyr_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  ! Broadcast the chemical potential of the environment
  call mpi_bcast(avg_mu_env,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine distrib_pf
