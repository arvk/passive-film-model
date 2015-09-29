subroutine kMC_filmdissolve(iter,simstate)
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

  PetscInt :: ierr, my_rank
  integer (C_INT) :: myargc
  character*1 , target :: myargv
  type(c_ptr), target :: myspparks
  PetscInt, dimension(psx_g,psy_g) :: interface_loc
  PetscScalar, dimension(psx_g,psy_g) :: distance_interface_moved
  PetscScalar, dimension(psx_g,psy_g) :: vac_form_bias
  PetscInt :: x,y,z,xfine,yfine
  PetscInt :: nint
  character*24 :: kmc_numel_string
  PetscInt, dimension(psx_g*kg_scale,psy_g*kg_scale) :: fine_kmc_array
  PetscScalar :: average_from_fine
  type(context) simstate
  PetscScalar no_of_env_cells
  PetscInt, intent(in) :: iter  ! Iteration count
  DM latticemu
  Vec onlymus_petscorder, onlymus_naturalorder, onlymus_maxstride, onlymus_alongz
  PetscScalar, pointer :: statepointer(:,:,:,:)


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


  call DMDACreate3D(MPI_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
       & DMDA_STENCIL_STAR,psx_g,psy_g,psz_g,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1, &
       & 1,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,latticemu,ierr)

  call DMCreateGlobalVector(latticemu,onlymus_petscorder,ierr)
  call VecStrideGather(simstate%slice,nenv,onlymus_petscorder,INSERT_VALUES,ierr)

  call DMDACreateNaturalVector(latticemu,onlymus_naturalorder,ierr)
  call DMDAGlobalToNaturalBegin(latticemu,onlymus_petscorder,INSERT_VALUES,onlymus_naturalorder,ierr)
  call DMDAGlobalToNaturalEnd(latticemu,onlymus_petscorder,INSERT_VALUES,onlymus_naturalorder,ierr)

  call VecCreate(MPI_COMM_WORLD,onlymus_maxstride,ierr)
  call VecSetBlockSize(onlymus_maxstride,psx_g*psy_g,ierr)
  call VecSetSizes(onlymus_maxstride,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(onlymus_maxstride,ierr)
  call VecCopy(onlymus_naturalorder,onlymus_maxstride,ierr)


  call VecCreate(MPI_COMM_WORLD,onlymus_alongz,ierr)
  call VecSetSizes(onlymus_alongz,PETSC_DECIDE,psz_g,ierr)
  call VecSetUp(onlymus_alongz,ierr)

     do x = 1,psx_g
        do y = 1,psy_g
           call VecStrideGather(onlymus_maxstride,(psx_g*(y-1))+(x-1),onlymus_alongz,INSERT_VALUES,ierr)
           call VecSum(onlymus_alongz,no_of_env_cells,ierr)
           interface_loc(x,y) = psz_g - no_of_env_cells
           vac_form_bias(x,y) = (1.0d0-sin(3.14d0*x*y/800.0d0))*0.20d0
        end do
     end do


  if(isroot)then

     call system('rm -f input.filmenv')

     open(unit = 666, file = 'input.filmenv', status = 'new')
     write(666,*) 'Testing'
     write(666,*) '2 dimension'
     write(666,*) '0 ', psx_g*kg_scale ,' xlo xhi'
     write(666,*) '0 ', psy_g*kg_scale ,' ylo yhi'
     write(666,*) '-0.5 0.5 zlo zhi'
     write(666,*) psx_g*psy_g*kg_scale*kg_scale, ' sites'
     write(666,*) ' '
     write(666,*) 'Values'
     write(666,*) ' '

     do x = 0,psx_g-1
        do y = 0,psy_g-1
           do xfine = 0,kg_scale-1
              do yfine = 0,kg_scale-1
                 write(666,'(I8,I3,I7,F19.12)') 1+((x*kg_scale)+xfine)+(((y*kg_scale)+yfine)*(kg_scale*psx_g)), 1, interface_loc(x+1,y+1)*kg_scale, vac_form_bias(x+1,y+1)
              end do
           end do
        end do
     end do

     close(666)

  end if




  if(isroot)then

     call system('rm -f filmenv.spparksscript')

     open(unit = 667, file = 'filmenv.spparksscript', status = 'new')
     write(667,*) 'seed 1273'
     write(667,*) 'app_style filmenv 0.25'
     write(667,*) 'dimension 2'
     write(667,*) 'boundary p p p'
     write(667,*) 'lattice sq/4n 1.0'
     write(667,'(A,F16.8)') 'temperature ', 8.61733034E-5*T
     write(667,*) 'region cell block 0 ',psx_g*kg_scale,' 0 ',psy_g*kg_scale,' -0.5 0.5'
     write(667,*) 'create_box cell'
     write(667,*) 'create_sites box'
     write(667,*) 'read_sites input.filmenv'
     write(667,*) 'solve_style tree'
     write(667,*) 'sector yes'
     write(667,*) 'diag_style energy stats yes'
     write(667,*) 'stats 5'
     write(667,*) 'dump couplingfe text 15 couplingfe x y i2'
     write(667,'(A,F16.8)') 'run ', dt*kmc_freq
     write(667,*) 'undump couplingfe'
     close(667)

  end if



  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'filmenv.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)



  if(isroot)then
     call system('rm -f input.filmenv log.spparks')
     write(kmc_numel_string,'(I24)') psx_g*psy_g*kg_scale*kg_scale
     call system('tail -n '//trim(kmc_numel_string)//' couplingfe | awk ''{printf " %5.5i %5.5i %5.5i \n", $1,$2,$3}'' > tocouple')
     call system('rm -f couplingfe')
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Wait for RANK 0 to write the 'toucouple' file before reading it

     open (unit = 667, file = 'tocouple', status = 'old')
     do x = 1,psx_g*psy_g*kg_scale*kg_scale
        read(667,'(I6, I6, I6)') xfine, yfine, fine_kmc_array(xfine+1,yfine+1)
     end do
     close(667)

     do x = 0,psx_g-1
        do y = 0,psy_g-1
           average_from_fine = 0.0d0
           do xfine = 0,kg_scale-1
              do yfine = 0,kg_scale-1
                 average_from_fine = average_from_fine + fine_kmc_array((x*kg_scale)+xfine+1,(y*kg_scale)+yfine+1)
              end do
           end do
           average_from_fine = average_from_fine/(kg_scale*kg_scale)
           distance_interface_moved(x+1,y+1) = interface_loc(x+1,y+1) - (average_from_fine/kg_scale)
           interface_loc(x+1,y+1) = floor(average_from_fine/kg_scale)
        end do
     end do


     call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

     do x = simstate%startx , simstate%startx + simstate%widthx-1
        do y = simstate%starty , simstate%starty + simstate%widthy-1
           do z = max(simstate%startz,interface_loc(x+1,y+1)) , simstate%startz + simstate%widthz-1
              statepointer(nmet,x,y,z) = 0.0d0
              statepointer(nmkw,x,y,z) = 0.0d0
              statepointer(npht,x,y,z) = 0.0d0
              statepointer(npyr,x,y,z) = 0.0d0
              statepointer(nenv,x,y,z) = 1.0d0
              statepointer(nmus,x,y,z) = avg_mu_env
              statepointer(nvoi,x,y,z) = 1.0d0
           end do
        end do
     end do

     call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)


  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Wait for ALL RANKS to complete coupling the kMC to PF before erasing the 'toucouple' file
  if(isroot) call system('rm -f tocouple')
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before going back to the PF routines

end subroutine kMC_filmdissolve
