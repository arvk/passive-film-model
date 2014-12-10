subroutine gather_kmc()
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y   ! Loop variables
  integer :: rank_loop
  integer :: ierr,status(MPI_STATUS_SIZE)   
  
  if (rank.eq.0) then

     do x = 1,ksx+2
        do y = 1,ksy
              kg_recv(x,y) = kg(x,y)
        end do
     end do

     do rank_loop = 1,procs-1
        call mpi_recv(kg_recv(1,(rank_loop*ksy_g/procs)+1),(ksx+2)*ksy,MPI_INT,rank_loop,21,MPI_COMM_WORLD,status,ierr)
     end do

     do x = 1,ksx_g
        do y = 1,ksy_g
           kg_g(x,y) = kg_recv(x+1,y)
        end do
     end do

  else

     call mpi_send(kg(1,2),(ksx+2)*ksy,MPI_INT,0,21,MPI_COMM_WORLD,ierr)

  end if
  
end subroutine gather_kmc




