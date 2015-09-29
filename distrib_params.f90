subroutine distrib_params()
  use commondata
  use fields
  implicit none
#include <finclude/petscsys.h>
  !!####Distribute simulation parameters to all child processors

  PetscErrorCode ierr !! MPI error flag

  call mpi_bcast(psx_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psy_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psz_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(T,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pH_in,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nomc,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(noimg,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(include_dissolve,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(include_electro,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(metal_potential,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine distrib_params
