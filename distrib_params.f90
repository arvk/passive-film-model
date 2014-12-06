subroutine distrib_params()
  use commondata
  implicit none
  include 'mpif.h'

  integer :: ierr,status(MPI_STATUS_SIZE)   

  call mpi_bcast(psx_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psy_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psz_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psx,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psy,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(psz,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ksx,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ksy,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ksx_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ksy_g,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(T,1,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pH_in,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nomc,1,MPI_INT,0,MPI_COMM_WORLD,ierr)


end subroutine distrib_params
