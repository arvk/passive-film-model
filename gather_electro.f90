subroutine gather_electro
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
              elpot_g(x,y,z) = elpot(x,y,z+ghost_width)
           end do
        end do
     end do


     do rank_loop = 1,procs-1
        call mpi_recv(elpot_g(1,1,(rank_loop*psz)+1),psx*psy*psz,MPI_DOUBLE_PRECISION,rank_loop,1,MPI_COMM_WORLD,status,ierr)
     end do

  else

     call mpi_send(elpot(1,1,1+ghost_width),psx*psy*psz,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,ierr)

  end if

end subroutine gather_electro




