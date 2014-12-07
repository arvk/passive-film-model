subroutine distrib_pf()
  use commondata
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: ierr,status(MPI_STATUS_SIZE)   

  call mpi_bcast(met_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pht_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(env_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(pyr_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(mu_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ph_g(1,1,1),psx_g*psy_g*psz_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           met(x,y,z) = met_g(x,y,z-1+(rank*psz_g/procs))
           pht(x,y,z) = pht_g(x,y,z-1+(rank*psz_g/procs))
           env(x,y,z) = env_g(x,y,z-1+(rank*psz_g/procs))
           pyr(x,y,z) = pyr_g(x,y,z-1+(rank*psz_g/procs))
           mu(x,y,z) = mu_g(x,y,z-1+(rank*psz_g/procs))
           ph(x,y,z) = ph_g(x,y,z-1+(rank*psz_g/procs))
        end do
     end do
  end do

  call mpi_bcast(min_mu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(max_mu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(avg_mu_met,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(avg_mu_pht,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(avg_mu_env,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

end subroutine distrib_pf
