subroutine spparks_filmenv()
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
  real*8, dimension(psx,psy) :: vac_form_bias
  integer :: x,y,z,xfine,yfine

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

  call gather_pf()

  if(isroot)then
     interface_loc = 0

     do x = 1,psx
        do y = 1,psy
           do z = 1,psz_g
              if ((env_g(x,y,z) .lt. 5.0E-1).and.(env_g(x,y,z+1) .gt. 5.0E-1)) then
                 interface_loc(x,y) = z
                 vac_form_bias(x,y) = sin(x*y/5.0d0)+1.0d0
                 exit
              end if
           end do
        end do
     end do

     call system('rm input.filmenv')

     open(unit = 666, file = 'input.filmenv', status = 'new')
     write(666,*) 'Testing'
     write(666,*) '2 dimension'
     write(666,*) '0 ', psx*kg_scale ,' xlo xhi'
     write(666,*) '0 ', psy*kg_scale ,' ylo yhi'
     write(666,*) '-0.5 0.5 zlo zhi'
     write(666,*) psx*psy*kg_scale*kg_scale, ' sites'
     write(666,*) ' '
     write(666,*) 'Values'
     write(666,*) ' '

     do x = 0,psx-1
        do y = 0,psy-1
           do xfine = 0,kg_scale-1
              do yfine = 0,kg_scale-1
                 write(666,'(I8,I3,I7,F9.6)') 1+((x*kg_scale)+xfine)+(((y*kg_scale)+yfine)*(kg_scale*psx)), 0, interface_loc(x+1,y+1), vac_form_bias(x+1,y+1)
              end do
           end do
        end do
     end do

     close(666)

  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'filmenv.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)


end subroutine spparks_filmenv
