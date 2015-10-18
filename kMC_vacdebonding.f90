subroutine kMC_vacdebonding(iter,simstate,metal_content_in_simcell,random_context)
  use, intrinsic :: iso_c_binding
  use commondata
  use fields
  use diffusion_constants
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
  PetscInt :: site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
  character*24 :: kmc_numel_string
  PetscInt, dimension(psx_g*kg_scale,psy_g*kg_scale) :: fine_kmc_array
  PetscScalar :: average_from_fine
  type(context) simstate
  PetscInt, intent(in) :: iter  ! Iteration count
  PetscInt, dimension(psx_g,psy_g) :: coarse_vac_config
  PetscInt :: coarsex,coarsey
  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscInt :: floor
  PetscScalar :: max, min
  PetscScalar, intent(inout) :: metal_content_in_simcell    !! Amount of metal phase in the simulation cell
  PetscScalar :: metal_content_in_simcell_last_timestep     !! Amount of metal phase in the simulation cell
  PetscScalar :: DeltaC
  PetscRandom, intent(inout) :: random_context !! Context to seed and generate random numbers
  PetscReal :: random_number  !! Pseudo random number generated from a PETSc context

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

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  metal_content_in_simcell_last_timestep = metal_content_in_simcell
  call VecStrideNorm(simstate%slice,nmet,NORM_1,metal_content_in_simcell,ierr)

  ! Calculate \DeltaC => Vacancy concentration at the metal-mkw interface
  DeltaC = ((metal_content_in_simcell_last_timestep - metal_content_in_simcell)*(10*dpf*dpf))/(D_Fe_met*psx_g*psy_g*dt*kmc_freq)

  call PetscRandomGetValueReal(random_context,random_number,ierr)

  if(isroot)then

     call system('rm -f vacmet.spparksscript')
     open(unit = 666, file = 'vacmet.spparksscript', status = 'new')
     write(666,*) 'seed ', floor(100000.0d0*random_number)+1
     write(666,*) 'app_style diffusion nonlinear hop'
     write(666,*) 'dimension 2'
     write(666,*) 'boundary p p p'
     write(666,*) 'lattice sq/4n 1.0'
     write(666,'(A,F16.8)') 'temperature ', 8.61733034E-5*T
     write(666,*) 'region cell block 0 ',psx_g*kg_scale,' 0 ',psy_g*kg_scale,' -0.5 0.5'
     write(666,*) 'create_box cell'
     write(666,*) 'create_sites box'
     write(666,*) 'set i1 value 2 region cell'
     write(666,*) 'ecoord 0  0.20'
     write(666,*) 'ecoord 1  0.15'
     write(666,*) 'ecoord 2  0.05'
     write(666,*) 'ecoord 3  0.00'
     write(666,*) 'ecoord 4 -0.05'
     write(666,*) 'barrier hop 0.05'
     write(666,*) 'solve_style tree'
     write(666,*) 'sector yes'
     write(666,'(A,F16.8,A)') 'dump mydump text ', dt*kmc_freq, ' raw_vacmet_output x y i1'

     do x = 0,psx_g-1
        do y = 0,psy_g-1
           call PetscRandomGetValueReal(random_context,random_number,ierr)
           write(666,'(A,F6.3,A,I5,A,I5,A,I5,A,I5)') ' set i1 value 1 fraction ', max(min(DeltaC*2.0d0*random_number,0.99d0),0.01d0), ' if x > ',x*kg_scale,  ' if x < ',(x+1)*kg_scale,  ' if y > ',y*kg_scale,  ' if y < ',(y+1)*kg_scale
        end do
     end do

     write(666,'(A,F16.8)') 'dump_modify mydump delay ', dt*kmc_freq
     write(666,'(A,F16.8)') 'run ', dt*kmc_freq
     close(666)

  end if
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'vacmet.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  if(isroot)then
     call system('rm -f input.vacmet log.spparks')
     ! READ RAW SPPARKS OUTPUT AND WRITE FORTRAN-READABLE OUTPUT
     write(kmc_numel_string,'(I24)') psx_g*psy_g*kg_scale*kg_scale
     call system('tail -n '//trim(kmc_numel_string)//' raw_vacmet_output | awk ''{printf " %4.4i %4.4i %3.3i \n", $1,$2,$3}'' > vacmet_spparks_output.template')
     call system('rm -f raw_vacmet_output')
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions


  call system('rm -f input.vacmet')
  open (unit = 667, file = 'vacmet_spparks_output.template', status = 'old')

  coarse_vac_config = 0.0d0

  do x = 1,psx_g*psy_g*kg_scale*kg_scale
     read(667,'(I5, I5, I4)') partial_x, partial_y, i1
     coarsex = floor(1.0d0*partial_x/kg_scale)
     coarsey = floor(1.0d0*partial_y/kg_scale)
     coarse_vac_config(coarsex+1,coarsey+1) = coarse_vac_config(coarsex+1,coarsey+1) + (i1/(kg_scale*kg_scale))
  end do

  close(667)

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  do x = simstate%startx , simstate%startx + simstate%widthx-1
     do y = simstate%starty , simstate%starty + simstate%widthy-1
        do z = simstate%startz + 1 , simstate%startz + simstate%widthz-1
           if (((statepointer(nmet,x,y,z-1).gt.0.5d0).and.(statepointer(nmkw,x,y,z-1).lt.0.5d0)) .and. ((statepointer(nmet,x,y,z).lt.0.5d0).and.(statepointer(nmkw,x,y,z).gt.0.5d0))) then
              statepointer(nvoi,x,y,z-1) = statepointer(nvoi,x,y,z-1) - coarse_vac_config(x+1,y+1)
              statepointer(nvoi,x,y,z-1) = min(max(statepointer(nvoi,x,y,z-1),0.0d0),1.0d0)
           end if
        end do
     end do
  end do
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)


  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before file IO
  if(isroot)then
     call system('rm -f vacmet_spparks_output.template')
  end if
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before returning to PF routines

end subroutine kMC_vacdebonding
