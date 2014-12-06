subroutine swap_pf()
  use commondata
  implicit none
  include 'mpif.h'

  integer :: x, y   ! Loop variables
  integer :: ierr
  integer :: stat1(MPI_STATUS_SIZE),stat2(MPI_STATUS_SIZE),stat3(MPI_STATUS_SIZE),stat4(MPI_STATUS_SIZE),stat5(MPI_STATUS_SIZE)
  integer :: stat6(MPI_STATUS_SIZE),stat7(MPI_STATUS_SIZE),stat8(MPI_STATUS_SIZE),stat9(MPI_STATUS_SIZE),stat10(MPI_STATUS_SIZE)
  integer :: stat11(MPI_STATUS_SIZE),stat12(MPI_STATUS_SIZE),stat13(MPI_STATUS_SIZE),stat14(MPI_STATUS_SIZE),stat15(MPI_STATUS_SIZE)
  integer :: stat16(MPI_STATUS_SIZE),stat17(MPI_STATUS_SIZE),stat18(MPI_STATUS_SIZE),stat19(MPI_STATUS_SIZE),stat20(MPI_STATUS_SIZE)
  integer :: req1,req2,req3,req4,req5,req6,req7,req8,req9,req10
  integer :: req11,req12,req13,req14,req15,req16,req17,req18,req19,req20

  if ((rank.gt.0).and.(rank.lt.procs-1)) then
     call mpi_isend(met(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,1,MPI_COMM_WORLD,req1,ierr)
     call mpi_isend(pht(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,3,MPI_COMM_WORLD,req3,ierr)
     call mpi_isend(env(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,5,MPI_COMM_WORLD,req5,ierr)
     call mpi_isend(mu(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,7,MPI_COMM_WORLD,req7,ierr)
     call mpi_isend(ph(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,9,MPI_COMM_WORLD,req9,ierr)

     call mpi_isend(met(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,2,MPI_COMM_WORLD,req2,ierr)
     call mpi_isend(pht(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,4,MPI_COMM_WORLD,req4,ierr)
     call mpi_isend(env(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,6,MPI_COMM_WORLD,req6,ierr)
     call mpi_isend(mu(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,8,MPI_COMM_WORLD,req8,ierr)
     call mpi_isend(ph(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,req10,ierr)

     call mpi_irecv(met(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,1,MPI_COMM_WORLD,req11,ierr)
     call mpi_irecv(pht(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,3,MPI_COMM_WORLD,req13,ierr)
     call mpi_irecv(env(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,5,MPI_COMM_WORLD,req15,ierr)
     call mpi_irecv(mu(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,7,MPI_COMM_WORLD,req17,ierr)
     call mpi_irecv(ph(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,9,MPI_COMM_WORLD,req19,ierr)

     call mpi_irecv(met(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,2,MPI_COMM_WORLD,req12,ierr)
     call mpi_irecv(pht(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,4,MPI_COMM_WORLD,req14,ierr)
     call mpi_irecv(env(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,6,MPI_COMM_WORLD,req16,ierr)
     call mpi_irecv(mu(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,8,MPI_COMM_WORLD,req18,ierr)
     call mpi_irecv(ph(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,req20,ierr)

     call mpi_wait(req1,stat1,ierr)
     call mpi_wait(req2,stat1,ierr)
     call mpi_wait(req3,stat1,ierr)
     call mpi_wait(req4,stat1,ierr)
     call mpi_wait(req5,stat1,ierr)
     call mpi_wait(req6,stat1,ierr)
     call mpi_wait(req7,stat1,ierr)
     call mpi_wait(req8,stat1,ierr)
     call mpi_wait(req9,stat1,ierr)
     call mpi_wait(req10,stat1,ierr)
     call mpi_wait(req11,stat1,ierr)
     call mpi_wait(req12,stat1,ierr)
     call mpi_wait(req13,stat1,ierr)
     call mpi_wait(req14,stat1,ierr)
     call mpi_wait(req15,stat1,ierr)
     call mpi_wait(req16,stat1,ierr)
     call mpi_wait(req17,stat1,ierr)
     call mpi_wait(req18,stat1,ierr)
     call mpi_wait(req19,stat1,ierr)
     call mpi_wait(req20,stat1,ierr)

  elseif (rank.eq.0) then

     call mpi_irecv(met(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,2,MPI_COMM_WORLD,req12,ierr)
     call mpi_irecv(pht(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,4,MPI_COMM_WORLD,req14,ierr)
     call mpi_irecv(env(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,6,MPI_COMM_WORLD,req16,ierr)
     call mpi_irecv(mu(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,8,MPI_COMM_WORLD,req18,ierr)
     call mpi_irecv(ph(1,1,psz+2),psx*psy,MPI_DOUBLE_PRECISION,rank+1,10,MPI_COMM_WORLD,req20,ierr)

     call mpi_isend(met(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,1,MPI_COMM_WORLD,req1,ierr)
     call mpi_isend(pht(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,3,MPI_COMM_WORLD,req3,ierr)
     call mpi_isend(env(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,5,MPI_COMM_WORLD,req5,ierr)
     call mpi_isend(mu(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,7,MPI_COMM_WORLD,req7,ierr)
     call mpi_isend(ph(1,1,psz+1),psx*psy,MPI_DOUBLE_PRECISION,rank+1,9,MPI_COMM_WORLD,req9,ierr)

     do x = 1,psx
        do y = 1,psy
           met(x,y,1) = met(x,y,2)
           pht(x,y,1) = pht(x,y,2)
           env(x,y,1) = env(x,y,2)
           mu(x,y,1) = mu(x,y,2)
           ph(x,y,1) = ph(x,y,2)
        end do
     end do

     call mpi_wait(req1,stat1,ierr)
     call mpi_wait(req12,stat1,ierr)
     call mpi_wait(req3,stat1,ierr)
     call mpi_wait(req14,stat1,ierr)
     call mpi_wait(req5,stat1,ierr)
     call mpi_wait(req16,stat1,ierr)
     call mpi_wait(req7,stat1,ierr)
     call mpi_wait(req18,stat1,ierr)
     call mpi_wait(req9,stat1,ierr)
     call mpi_wait(req20,stat1,ierr)

  else

     call mpi_irecv(met(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,1,MPI_COMM_WORLD,req11,ierr)
     call mpi_irecv(pht(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,3,MPI_COMM_WORLD,req13,ierr)
     call mpi_irecv(env(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,5,MPI_COMM_WORLD,req15,ierr)
     call mpi_irecv(mu(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,7,MPI_COMM_WORLD,req17,ierr)
     call mpi_irecv(ph(1,1,1),psx*psy,MPI_DOUBLE_PRECISION,rank-1,9,MPI_COMM_WORLD,req19,ierr)

     call mpi_isend(met(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,2,MPI_COMM_WORLD,req2,ierr)
     call mpi_isend(pht(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,4,MPI_COMM_WORLD,req4,ierr)
     call mpi_isend(env(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,6,MPI_COMM_WORLD,req6,ierr)
     call mpi_isend(mu(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,8,MPI_COMM_WORLD,req8,ierr)
     call mpi_isend(ph(1,1,2),psx*psy,MPI_DOUBLE_PRECISION,rank-1,10,MPI_COMM_WORLD,req10,ierr)

     do x = 1,psx
        do y = 1,psy
           met(x,y,psz+2) = met(x,y,psz+1)
           pht(x,y,psz+2) = pht(x,y,psz+1)
           env(x,y,psz+2) = env(x,y,psz+1)
           mu(x,y,psz+2) = mu(x,y,psz+1)
           ph(x,y,psz+2) = ph(x,y,psz+1)
        end do
     end do

     call mpi_wait(req11,stat1,ierr)
     call mpi_wait(req2,stat1,ierr)
     call mpi_wait(req13,stat1,ierr)
     call mpi_wait(req4,stat1,ierr)
     call mpi_wait(req15,stat1,ierr)
     call mpi_wait(req6,stat1,ierr)
     call mpi_wait(req17,stat1,ierr)
     call mpi_wait(req8,stat1,ierr)
     call mpi_wait(req19,stat1,ierr)
     call mpi_wait(req10,stat1,ierr)

  end if


end subroutine swap_pf

