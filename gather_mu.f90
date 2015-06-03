subroutine gather_mu
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: rank_loop
  integer :: ierr,status(MPI_STATUS_SIZE)

  if (rank.eq.0) then

     do x = 1,psx
        do y = 1,psy
           do z = 1,psz
              mu_g(x,y,z) = mu(x,y,z+ghost_width)
           end do
        end do
     end do


     do rank_loop = 1,procs-1
        call mpi_recv(mu_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,1,MPI_COMM_WORLD,status,ierr)
     end do

  else

     call mpi_send(mu(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)

  end if

end subroutine gather_mu




