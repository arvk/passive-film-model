subroutine swap_electro
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y   ! Loop variables
  integer :: ierr
  integer :: stat(MPI_STATUS_SIZE)

  !! Triple number codes for MPI request variables
  !! a) First number denotes field: 1 - elpot
  !! b) Second number denotes sending/receiving: 1 - Sending, 2 - Receiving
  !! c) Third number denotes target: 0 - Lower rank, 1 - Higher rank (Rank 0 is assumed to be higher when data is transferred to it from rank N-1)
  integer :: q110, q111, q120, q121

  !! Double number codes for MPI tags
  !! a) First number denotes field: 1 - elpot
  !! b) Second number denotes sending/receiving and destination combo: 0 = Sending to higher or receiving from lower. 1 = Sending to lower rank or receiving from higher rank

  if ((rank.gt.0).and.(rank.lt.procs-1)) then

     call mpi_isend(elpot(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,q111,ierr)

     call mpi_isend(elpot(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,11,MPI_COMM_WORLD,q110,ierr)

     call mpi_irecv(elpot(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,q120,ierr)

     call mpi_irecv(elpot(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,11,MPI_COMM_WORLD,q121,ierr)

     call mpi_wait(q111,stat,ierr)
     call mpi_wait(q110,stat,ierr)
     call mpi_wait(q120,stat,ierr)
     call mpi_wait(q121,stat,ierr)

  elseif (rank.eq.0) then

     call mpi_isend(elpot(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,q111,ierr)

     call mpi_irecv(elpot(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,11,MPI_COMM_WORLD,q121,ierr)

     do x = 1,psx
        do y = 1,psy
           elpot(x,y,1) = elpot(x,y,2)
        end do
     end do

     call mpi_wait(q111,stat,ierr)
     call mpi_wait(q121,stat,ierr)

  else

     call mpi_isend(elpot(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,11,MPI_COMM_WORLD,q110,ierr)

     call mpi_irecv(elpot(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,q120,ierr)

     do x = 1,psx
        do y = 1,psy
           elpot(x,y,psz+2) = elpot(x,y,psz+1)
        end do
     end do

     call mpi_wait(q110,stat,ierr)
     call mpi_wait(q120,stat,ierr)

  end if


end subroutine swap_electro

