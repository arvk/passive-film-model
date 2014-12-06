subroutine gather_pf()
  use commondata
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: rank_loop
  integer :: ierr,status(MPI_STATUS_SIZE)   
  
  if (rank.eq.0) then

     do x = 1,psx
        do y = 1,psy
           do z = 1,psz
              met_g(x,y,z) = met(x,y,z+1)
              pht_g(x,y,z) = pht(x,y,z+1)
              env_g(x,y,z) = env(x,y,z+1)
              mu_g(x,y,z) = mu(x,y,z+1)
              ph_g(x,y,z) = ph(x,y,z+1)
           end do
        end do
     end do

     do rank_loop = 1,procs-1
        call mpi_recv(met_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,1,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(pht_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,3,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(env_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,5,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(mu_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,7,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(ph_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,9,MPI_COMM_WORLD,status,ierr)
     end do

  else

     call mpi_send(met(1,1,2),psx*psy*psz,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
     call mpi_send(pht(1,1,2),psx*psy*psz,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
     call mpi_send(env(1,1,2),psx*psy*psz,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
     call mpi_send(mu(1,1,2),psx*psy*psz,MPI_DOUBLE_PRECISION,0,7,MPI_COMM_WORLD,ierr)
     call mpi_send(ph(1,1,2),psx*psy*psz,MPI_DOUBLE_PRECISION,0,9,MPI_COMM_WORLD,ierr)

  end if
  
end subroutine gather_pf




