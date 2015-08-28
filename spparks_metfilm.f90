subroutine spparks_metfilm()
  use, intrinsic :: iso_c_binding
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: ierr, my_rank

  integer (C_INT) :: myargc
  character*1 , target :: myargv
  type(c_ptr), target :: myspparks

  integer, dimension(psx,psy) :: interface_loc
  real*8, dimension(psx,psy) :: distance_interface_moved
  real*8, dimension(psx,psy) :: vac_form_bias
  integer :: x,y,z,xfine,yfine
  integer :: nint

  integer,dimension(psx*kg_scale*psy*kg_scale*4) :: site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved

  character*24 :: kmc_numel_string
  integer, dimension(psx*kg_scale,psy*kg_scale) :: fine_kmc_array
  real*8 :: average_from_fine

  logical :: already_run_once

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

  call gather_pf()

  if(isroot)then

     call system('cp metfilm.spparksscript.template metfilm.spparksscript')

     if (already_run_once) then
!        call system('echo "read_sites input.metfilm" >> metfilm.spparksscript')
     end if

     call system('echo "run 10000" >> metfilm.spparksscript')

  end if



  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'metfilm.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)







  if(isroot)then

     call system('rm -f input.metfilm log.spparks')

     ! READ RAW SPPARKS OUTPUT AND WRITE FORTRAN-READABLE OUTPUT
     write(kmc_numel_string,'(I24)') psx*psy*kg_scale*kg_scale*4
     call system('tail -n '//trim(kmc_numel_string)//' raw_metfilm_output | awk ''{printf " %7.7i %4.4i %4.4i %3.3i %3.3i %5.5i \n", $1,$2,$3,$4,$5,$6}'' > metfilm_spparks_output.template')
     call system('rm -f raw_metfilm_output')

     ! READ FORTRAN-READABLE INPUT AND STORE IN ARRAY FOR COUPLING
     open (unit = 667, file = 'metfilm_spparks_output.template', status = 'old')
     do x = 1,psx*psy*kg_scale*kg_scale*4
        read(667,'(I8, I5, I5, I4, I4, I6)') site_id(x), partial_x(x), partial_y(x), i1(x), i2(x), no_h2_already_evolved(x)
     end do
     close(667)
     call system('rm -f metfilm_spparks_output.template')

     ! WRITE SPPARKS-READABLE OUTPUT
     call system('rm -f input.metfilm')
     open(unit = 666, file = 'input.metfilm', status = 'new')
     write(666,*) 'Testing'
     write(666,*) '2 dimension'
     write(666,*) '0 ', psx*kg_scale ,' xlo xhi'
     write(666,*) '0 ', psy*kg_scale ,' ylo yhi'
     write(666,*) '-0.5 0.5 zlo zhi'
     write(666,*) psx*psy*kg_scale*kg_scale*4, ' sites'
     write(666,*) ' '
     write(666,*) 'Values'
     write(666,*) ' '

     do x = 1,psx*psy*kg_scale*kg_scale*4
        write(666,'(I8,I4,I4,I6)') site_id(x), i1(x), i2(x), no_h2_already_evolved(x)
     end do

     close(666)


     ! do x = 1,psx*psy*kg_scale*kg_scale*4
     !    fine_kmc_array(partial_x(x)+1,partial_y(x)+1) = fine_kmc_array(partial_x(x),partial_y(x)) + no_h2_already_evolved(x)
     ! end do    


  end if


end subroutine spparks_metfilm
