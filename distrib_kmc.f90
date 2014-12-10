subroutine distrib_kmc()
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: ierr,status(MPI_STATUS_SIZE)   

  call mpi_bcast(kg_g(1,1),ksx_g*ksy_g,MPI_INT,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(vfe_f_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(vfe_a_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(vs_f_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(vs_a_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(fes_diss_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(v_diff_g(1,1),ksx_g*ksy_g,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  do x = 1,ksx
     do y = 1,ksy
        kg(x+1,y+1) = kg_g(x,y+(rank*ksy_g/procs))
        vfe_f(x+1,y+1) = vfe_f_g(x,y+(rank*ksy_g/procs))
        vfe_a(x+1,y+1) = vfe_a_g(x,y+(rank*ksy_g/procs))
        vs_f(x+1,y+1) = vs_f_g(x,y+(rank*ksy_g/procs))
        vs_a(x+1,y+1) = vs_a_g(x,y+(rank*ksy_g/procs))
        fes_diss(x+1,y+1) = fes_diss_g(x,y+(rank*ksy_g/procs))
        v_diff(x+1,y+1) = v_diff_g(x,y+(rank*ksy_g/procs))
     end do
  end do

end subroutine distrib_kmc
