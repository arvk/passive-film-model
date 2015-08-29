subroutine gather_pf
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y, z                      ! Index for x-, y-, and z-direction (Loop)
  integer :: rank_loop                    ! Integer equal to processor rank (Loop)
  integer :: ierr,status(MPI_STATUS_SIZE) ! MPI variables

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  if (rank.eq.0) then

     do x = 1,psx
        do y = 1,psy
           do z = 1,psz
              met_g(x,y,z) = met(x,y,z+ghost_width)
              mkw_g(x,y,z) = mkw(x,y,z+ghost_width)
              pht_g(x,y,z) = pht(x,y,z+ghost_width)
              pyr_g(x,y,z) = pyr(x,y,z+ghost_width)
              env_g(x,y,z) = env(x,y,z+ghost_width)
              pH_g(x,y,z) = pH(x,y,z+ghost_width)
           end do
        end do
     end do


     do rank_loop = 1,procs-1
        call mpi_recv(met_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,1,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(pht_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,2,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(pyr_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,3,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(env_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,4,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(mkw_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,5,MPI_COMM_WORLD,status,ierr)
        call mpi_recv(pH_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,6,MPI_COMM_WORLD,status,ierr)
     end do

  else

     call mpi_send(met(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)
     call mpi_send(pht(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,ierr)
     call mpi_send(pyr(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,ierr)
     call mpi_send(env(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,4,MPI_COMM_WORLD,ierr)
     call mpi_send(mkw(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,5,MPI_COMM_WORLD,ierr)
     call mpi_send(pH(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,6,MPI_COMM_WORLD,ierr)

  end if

end subroutine gather_pf




