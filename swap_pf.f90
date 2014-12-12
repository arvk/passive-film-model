subroutine swap_pf
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y   ! Loop variables
  integer :: ierr
  integer :: stat(MPI_STATUS_SIZE)

  !! Triple number codes for MPI request variables
  !! a) First number denotes field: 1 - Met, 2 - Pht, 3 - Pyr, 4 - Env
  !! b) Second number denotes sending/receiving: 1 - Sending, 2 - Receiving
  !! c) Third number denotes target: 0 - Lower rank, 1 - Higher rank (Rank 0 is assumed to be higher when data is transferred to it from rank N-1)
  integer :: q110, q111, q120, q121
  integer :: q210, q211, q220, q221
  integer :: q310, q311, q320, q321
  integer :: q410, q411, q420, q421

  !! Double number codes for MPI tags
  !! a) First number denotes field: 1 - Met, 2 - Pht, 3 - Pyr, 4 - Env
  !! b) Second number denotes sending/receiving and destination combo: 0 = Sending to higher or receiving from lower. 1 = Sending to lower rank or receiving from higher rank

  if ((rank.gt.0).and.(rank.lt.procs-1)) then

     call mpi_isend(met(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,q111,ierr)
     call mpi_isend(pht(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,20,MPI_COMM_WORLD,q211,ierr)
     call mpi_isend(pyr(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,30,MPI_COMM_WORLD,q311,ierr)
     call mpi_isend(env(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,40,MPI_COMM_WORLD,q411,ierr)

     call mpi_isend(met(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,11,MPI_COMM_WORLD,q110,ierr)
     call mpi_isend(pht(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,21,MPI_COMM_WORLD,q210,ierr)
     call mpi_isend(pyr(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,31,MPI_COMM_WORLD,q310,ierr)
     call mpi_isend(env(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,41,MPI_COMM_WORLD,q410,ierr)

     call mpi_irecv(met(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,q120,ierr)
     call mpi_irecv(pht(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,20,MPI_COMM_WORLD,q220,ierr)
     call mpi_irecv(pyr(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,30,MPI_COMM_WORLD,q320,ierr)
     call mpi_irecv(env(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,40,MPI_COMM_WORLD,q420,ierr)

     call mpi_irecv(met(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,11,MPI_COMM_WORLD,q121,ierr)
     call mpi_irecv(pht(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,21,MPI_COMM_WORLD,q221,ierr)
     call mpi_irecv(pyr(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,31,MPI_COMM_WORLD,q321,ierr)
     call mpi_irecv(env(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,41,MPI_COMM_WORLD,q421,ierr)

     call mpi_wait(q111,stat,ierr); call mpi_wait(q211,stat,ierr); call mpi_wait(q311,stat,ierr); call mpi_wait(q411,stat,ierr) 
     call mpi_wait(q110,stat,ierr); call mpi_wait(q210,stat,ierr); call mpi_wait(q310,stat,ierr); call mpi_wait(q410,stat,ierr) 
     call mpi_wait(q120,stat,ierr); call mpi_wait(q220,stat,ierr); call mpi_wait(q320,stat,ierr); call mpi_wait(q420,stat,ierr) 
     call mpi_wait(q121,stat,ierr); call mpi_wait(q221,stat,ierr); call mpi_wait(q321,stat,ierr); call mpi_wait(q421,stat,ierr) 

  elseif (rank.eq.0) then

     call mpi_isend(met(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,q111,ierr)
     call mpi_isend(pht(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,20,MPI_COMM_WORLD,q211,ierr)
     call mpi_isend(pyr(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,30,MPI_COMM_WORLD,q311,ierr)
     call mpi_isend(env(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,40,MPI_COMM_WORLD,q411,ierr)

     call mpi_irecv(met(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,11,MPI_COMM_WORLD,q121,ierr)
     call mpi_irecv(pht(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,21,MPI_COMM_WORLD,q221,ierr)
     call mpi_irecv(pyr(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,31,MPI_COMM_WORLD,q321,ierr)
     call mpi_irecv(env(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,41,MPI_COMM_WORLD,q421,ierr)

     do x = 1,psx
        do y = 1,psy
           met(x,y,1) = met(x,y,2)
           pht(x,y,1) = pht(x,y,2)
           env(x,y,1) = env(x,y,2)
           pyr(x,y,1) = pyr(x,y,2)
        end do
     end do

     call mpi_wait(q111,stat,ierr); call mpi_wait(q211,stat,ierr); call mpi_wait(q311,stat,ierr); call mpi_wait(q411,stat,ierr) 
     call mpi_wait(q121,stat,ierr); call mpi_wait(q221,stat,ierr); call mpi_wait(q321,stat,ierr); call mpi_wait(q421,stat,ierr) 

  else

     call mpi_isend(met(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,11,MPI_COMM_WORLD,q110,ierr)
     call mpi_isend(pht(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,21,MPI_COMM_WORLD,q210,ierr)
     call mpi_isend(pyr(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,31,MPI_COMM_WORLD,q310,ierr)
     call mpi_isend(env(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,41,MPI_COMM_WORLD,q410,ierr)

     call mpi_irecv(met(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,q120,ierr)
     call mpi_irecv(pht(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,20,MPI_COMM_WORLD,q220,ierr)
     call mpi_irecv(pyr(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,30,MPI_COMM_WORLD,q320,ierr)
     call mpi_irecv(env(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,40,MPI_COMM_WORLD,q420,ierr)

     do x = 1,psx
        do y = 1,psy
           met(x,y,psz+2) = met(x,y,psz+1)
           pht(x,y,psz+2) = pht(x,y,psz+1)
           env(x,y,psz+2) = env(x,y,psz+1)
           pyr(x,y,psz+2) = pyr(x,y,psz+1)
        end do
     end do

     call mpi_wait(q110,stat,ierr); call mpi_wait(q210,stat,ierr); call mpi_wait(q310,stat,ierr); call mpi_wait(q410,stat,ierr) 
     call mpi_wait(q120,stat,ierr); call mpi_wait(q220,stat,ierr); call mpi_wait(q320,stat,ierr); call mpi_wait(q420,stat,ierr) 

  end if


end subroutine swap_pf

