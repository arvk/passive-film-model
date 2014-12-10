subroutine swap_kmc()
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y   ! Loop variables
  integer :: ierr
  integer :: stat1(MPI_STATUS_SIZE),stat2(MPI_STATUS_SIZE),stat3(MPI_STATUS_SIZE),stat4(MPI_STATUS_SIZE)
  integer :: req1,req2,req3,req4

  do y = 2,ksy+1
     kg(1,y) = kg(ksx+1,y)
     kg(ksx+2,y) = kg(2,y)
  end do

  if ((rank.gt.0).and.(rank.lt.procs-1)) then


     call mpi_isend(kg(1,ksy+1),ksx+2,MPI_INT,rank+1,11,MPI_COMM_WORLD,req1,ierr)
     call mpi_isend(kg(1,2),ksx+2,MPI_INT,rank-1,12,MPI_COMM_WORLD,req2,ierr)

     call mpi_irecv(kg(1,1),ksx+2,MPI_INT,rank-1,11,MPI_COMM_WORLD,req3,ierr)
     call mpi_irecv(kg(1,ksy+2),ksx+2,MPI_INT,rank+1,12,MPI_COMM_WORLD,req4,ierr)

  elseif (rank.eq.0) then

     call mpi_isend(kg(1,ksy+1),ksx+2,MPI_INT,rank+1,11,MPI_COMM_WORLD,req1,ierr)
     call mpi_isend(kg(1,2),ksx+2,MPI_INT,procs-1,12,MPI_COMM_WORLD,req2,ierr)

     call mpi_irecv(kg(1,1),ksx+2,MPI_INT,procs-1,11,MPI_COMM_WORLD,req3,ierr)
     call mpi_irecv(kg(1,ksy+2),ksx+2,MPI_INT,rank+1,12,MPI_COMM_WORLD,req4,ierr)

  else

     call mpi_isend(kg(1,ksy+1),ksx+2,MPI_INT,0,11,MPI_COMM_WORLD,req1,ierr)
     call mpi_isend(kg(1,2),ksx+2,MPI_INT,rank-1,12,MPI_COMM_WORLD,req2,ierr)

     call mpi_irecv(kg(1,1),ksx+2,MPI_INT,rank-1,11,MPI_COMM_WORLD,req3,ierr)
     call mpi_irecv(kg(1,ksy+2),ksx+2,MPI_INT,0,12,MPI_COMM_WORLD,req4,ierr)

  end if

     call mpi_wait(req1,stat1,ierr)
     call mpi_wait(req2,stat2,ierr)
     call mpi_wait(req3,stat3,ierr)
     call mpi_wait(req4,stat4,ierr)


end subroutine swap_kmc

