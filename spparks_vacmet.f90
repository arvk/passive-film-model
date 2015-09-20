subroutine spparks_vacmet(iter,simstate)
  use, intrinsic :: iso_c_binding
  use commondata
  use fields
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  integer :: ierr, my_rank
  integer (C_INT) :: myargc
  character*1 , target :: myargv
  type(c_ptr), target :: myspparks
  integer, dimension(psx_g,psy_g) :: interface_loc
  real*8, dimension(psx_g,psy_g) :: distance_interface_moved
  real*8, dimension(psx_g,psy_g) :: vac_form_bias
  integer :: x,y,z,xfine,yfine
  integer :: nint
  integer :: site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
  character*24 :: kmc_numel_string
  integer, dimension(psx_g*kg_scale,psy_g*kg_scale) :: fine_kmc_array
  real*8 :: average_from_fine
  logical :: already_run_once
  type(context) simstate
  integer, intent(in) :: iter  ! Iteration count
  integer, dimension(psx_g,psy_g) :: coarse_h2_evolved
  integer :: coarsex,coarsey
  PetscScalar, pointer :: statepointer(:,:,:,:)
  integer :: floor
  real*8 :: max, min

  interface

     subroutine spparks_close(instance) bind(C,NAME='spparks_close')
       use, intrinsic :: iso_c_binding
       import :: c_ptr
       implicit none
       type(c_ptr), value :: instance
     end subroutine spparks_close

     subroutine spparks_open(argc,argv,communicator,ptr) bind(C,NAME='spparks_open')
       use, intrinsic :: iso_c_binding
       import :: c_ptr
       implicit none
       integer (C_INT), value :: argc
       type(c_ptr), value :: argv
       type(c_ptr), value :: communicator
       type(c_ptr), value :: ptr
     end subroutine spparks_open

     subroutine spparks_open_no_mpi(argc,argv,ptr) bind(C,NAME='spparks_open_no_mpi')
       use, intrinsic :: iso_c_binding
       import :: c_ptr
       implicit none
       integer (C_INT), value :: argc
       type (c_ptr), value :: argv
       type (c_ptr), value :: ptr
     end subroutine spparks_open_no_mpi

     subroutine spparks_file(ptr,str) bind(C,NAME='spparks_file')
       use, intrinsic :: iso_c_binding
       import :: c_ptr
       implicit none
       type (c_ptr), value :: ptr
       character, dimension(*), intent(in) :: str
     end subroutine spparks_file

  end interface

  myargc = 0
  myargv = C_NULL_CHAR

  INQUIRE(FILE='input.metfilm', EXIST=already_run_once)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  if(isroot)then
     call system('cp vacmet.spparksscript.template vacmet.spparksscript')
     call system('echo "dump_modify mydump delay 1000" >> vacmet.spparksscript')
     call system('echo "run 1000" >> vacmet.spparksscript')
  end if
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'vacmet.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  ! if(isroot)then
  !    call system('rm -f input.metfilm log.spparks')
  !    ! READ RAW SPPARKS OUTPUT AND WRITE FORTRAN-READABLE OUTPUT
  !    write(kmc_numel_string,'(I24)') psx_g*psy_g*kg_scale*kg_scale*4
  !    call system('tail -n '//trim(kmc_numel_string)//' raw_metfilm_output | awk ''{printf " %7.7i %4.4i %4.4i %3.3i %3.3i %5.5i \n", $1,$2,$3,$4,$5,$6}'' > metfilm_spparks_output.template')
  !    call system('rm -f raw_metfilm_output')
  ! end if

  ! call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  ! if(isroot)then
  !    !  READ FORTRAN-READABLE INPUT AND STORE IN ARRAY FOR COUPLING

  !    call system('rm -f input.metfilm')
  !    open(unit = 666, file = 'input.metfilm', status = 'new')
  !    open (unit = 667, file = 'metfilm_spparks_output.template', status = 'old')

  !    write(666,*) 'Testing'
  !    write(666,*) '2 dimension'
  !    write(666,*) '0 ', psx_g*kg_scale ,' xlo xhi'
  !    write(666,*) '0 ', psy_g*kg_scale ,' ylo yhi'
  !    write(666,*) '-0.5 0.5 zlo zhi'
  !    write(666,*) psx_g*psy_g*kg_scale*kg_scale*4, ' sites'
  !    write(666,*) ' '
  !    write(666,*) 'Values'
  !    write(666,*) ' '

  !    do x = 1,psx_g*psy_g*kg_scale*kg_scale*4
  !       read(667,'(I8, I5, I5, I4, I4, I6)') site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
  !       write(666,'(I8,I4,I4,I6)') site_id, i1, i2, no_h2_already_evolved
  !    end do

  !    close(667)
  !    close(666)

  ! end if

  ! call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  ! coarse_h2_evolved = 0
  ! open (unit = 667, file = 'metfilm_spparks_output.template', status = 'old')
  ! do x = 1,psx_g*psy_g*kg_scale*kg_scale*4
  !    read(667,'(I8, I5, I5, I4, I4, I6)') site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
  !    coarsex = floor(1.0d0*partial_x/kg_scale)
  !    coarsey = floor(1.0d0*partial_y/kg_scale)
  !    coarse_h2_evolved(coarsex+1,coarsey+1) = coarse_h2_evolved(coarsex+1,coarsey+1) + no_h2_already_evolved
  ! end do
  ! close(667)

  ! call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  ! do x = simstate%startx , simstate%startx + simstate%widthx-1
  !    do y = simstate%starty , simstate%starty + simstate%widthy-1
  !       do z = simstate%startz , simstate%startz + simstate%widthz-2
  !          if (((statepointer(nmet,x,y,z).gt.0.5d0).and.(statepointer(nmkw,x,y,z).lt.0.5d0)) .and. ((statepointer(nmet,x,y,z+1).lt.0.5d0).and.(statepointer(nmkw,x,y,z+1).gt.0.5d0))) then
  !             statepointer(nvoi,x,y,z) = statepointer(nvoi,x,y,z) - (coarse_h2_evolved(x+1,y+1)/300.0d0)
  !             statepointer(nvoi,x,y,z) = min(max(statepointer(nvoi,x,y,z),0.0d0),1.0d0)
  !          end if
  !       end do
  !    end do
  ! end do
  ! call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)


  ! call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before file IO
  ! if(isroot)then
  !    call system('rm -f metfilm_spparks_output.template')
  ! end if
  ! call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before returning to PF routines

end subroutine spparks_vacmet
