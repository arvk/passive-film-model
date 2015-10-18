subroutine kMC_filmdissolve(iter,simstate,metal_content_in_simcell,random_context)
  use, intrinsic :: iso_c_binding
  use commondata
  use fields
  use thermo_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
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
  PetscInt, intent(in) :: iter  ! Iteration count
  PetscScalar, pointer :: statepointer(:,:,:,:)
  VecScatter :: gatherslicetoroot
  Vec :: gatheredslice
  PetscScalar, pointer :: gatheredpointer(:)
  Vec :: slice_naturalorder
  PetscScalar :: env_at_lattice_site_above, env_at_this_lattice_site, pH_at_this_lattice_site, pot_at_this_lattice_site
  PetscScalar, dimension(0:(nphases-1)) :: raw_dissolution_rate, vac_form_prob
  PetscInt :: fesphase
  PetscScalar, intent(inout) :: metal_content_in_simcell    !! Amount of metal phase in the simulation cell
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









  call DMDACreateNaturalVector(simstate%lattval,slice_naturalorder,ierr)
  call DMDAGlobalToNaturalBegin(simstate%lattval,simstate%slice,INSERT_VALUES,slice_naturalorder,ierr)
  call DMDAGlobalToNaturalEnd(simstate%lattval,simstate%slice,INSERT_VALUES,slice_naturalorder,ierr)

  call VecScatterCreateToAll(simstate%slice,gatherslicetoroot,gatheredslice,ierr)
  call VecScatterBegin(gatherslicetoroot,slice_naturalorder,gatheredslice,INSERT_VALUES,SCATTER_FORWARD,ierr)
  call VecScatterEnd(gatherslicetoroot,slice_naturalorder,gatheredslice,INSERT_VALUES,SCATTER_FORWARD,ierr)

  call VecGetArrayF90(gatheredslice,gatheredpointer,ierr)

  vac_form_bias = 0.0d0
  interface_loc = 0.0d0

  do x = 1,psx_g
     do y = 1,psy_g
        do z = psz_g-1,1,-1
           env_at_this_lattice_site = gatheredpointer(nfields*(((z-1)*psx_g*psy_g)+((y-1)*psx_g)+(x-1))+nenv+1)
           env_at_lattice_site_above = gatheredpointer(nfields*(((z-1+1)*psx_g*psy_g)+((y-1)*psx_g)+(x-1))+nenv+1)
           if ((env_at_lattice_site_above-env_at_this_lattice_site).gt.0.1d0) then
              interface_loc(x,y) = z

              pH_at_this_lattice_site = gatheredpointer(nfields*(((z-1)*psx_g*psy_g)+((y-1)*psx_g)+(x-1))+npH+1)
              pot_at_this_lattice_site = gatheredpointer(nfields*(((z-1)*psx_g*psy_g)+((y-1)*psx_g)+(x-1))+npot+1)

              ! Define raw dissolution rates (in nm/s)
              !---------------------------------------
              if (pH_at_this_lattice_site .gt. 1E-2) then
                 raw_dissolution_rate(nmet) = 1.781E7*((1E-14/pH_at_this_lattice_site)**0.6)*exp((0.85*Faraday_constant/(R*T))*(pot_at_this_lattice_site+0.45))
              else
                 raw_dissolution_rate(nmet) = 1.781E7*((1E-14/1E-2)**0.6)*exp((0.85*Faraday_constant/(R*T))*(pot_at_this_lattice_site+0.45))
              end if
              ! REF = Electrodissolution Kinetics of Iron in Chloride Solutions by Robert J. Chin* and Ken Nobe, JECS Vol 119, No. 11, Nov. 1972
              raw_dissolution_rate(nmkw) = 0.015d0
              ! REF = CORROSION MECHANISMS AND MATERIAL PERFORMANCE IN ENVIRONMENTS CONTAINING HYDROGEN SULFIDE AND ELEMENTAL SULFUR, Liane Smith and Bruce Craig, SACNUC Workshop 22nd and 23rd October, 2008, Brussels and References therein
              raw_dissolution_rate(npht) = 289.15*exp(0.0d0-(65900/(R*T)))*(10**(-1.46*pH_at_this_lattice_site))
              ! REF = "Pyrrhotite dissolution in acidic media" Chirita. P. et.al. Applied Geochemistry 41 (2014) 1-10. DOI: 10.1016/j.apgeochem.2013.11.013
              raw_dissolution_rate(npyr) = 0.00017244d0
              ! REF = "Interferometric study of pyrite surface reactivity in acidic conditions", Asta, M. P. et.al. American Mineralogist, Volume 93, pages 508–519, 2008
              raw_dissolution_rate(nenv) = 0.0d0

              ! Convert dissolution rates to probabilities for kMC
              do fesphase = 0,(nphases-1)
                 call PetscRandomGetValueReal(random_context,random_number,ierr)
                 vac_form_prob(fesphase) = max(((2.0d0*random_number*raw_dissolution_rate(fesphase))-0.0055d0)/0.00077d0,0.0d0)
              end do

              ! Final probabilities are weighted by phase fractions at the interface
              do fesphase = 0,(nphases-1)
                 vac_form_bias(x,y) = vac_form_bias(x,y) + (vac_form_prob(fesphase)*(gatheredpointer(nfields*(((z-1)*psx_g*psy_g)+((y-1)*psx_g)+(x-1))+fesphase+1)/(1.0d0-env_at_this_lattice_site)))
              end do

              exit
           end if
        end do
     end do
  end do

  call VecRestoreArrayF90(gatheredslice,gatheredpointer,ierr)

  call VecScatterDestroy(gatherslicetoroot,ierr)
  call VecDestroy(gatheredslice,ierr)
  call VecDestroy(slice_naturalorder,ierr)

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


  call PetscRandomGetValueReal(random_context,random_number,ierr)

  if(isroot)then

     call system('rm -f filmenv.spparksscript')

     open(unit = 667, file = 'filmenv.spparksscript', status = 'new')
     write(667,*) 'seed ', floor(100000.0d0*random_number)+1
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

        do z=simstate%startz,simstate%startz+simstate%widthz-1
           if ((statepointer(nenv,x,y,min(z+3,simstate%startz+simstate%widthz-1))-statepointer(nenv,x,y,max(z,simstate%startz))).gt.0.1d0) then
              statepointer(npH,x,y,min(z+3,simstate%startz+simstate%widthz-1)) = max(statepointer(npH,x,y,min(z+3,simstate%startz+simstate%widthz-1)) - (distance_interface_moved(x+1,y+1)*rhoFe),10**(0-pH_in))
           end if
        end do

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


  call VecStrideNorm(simstate%slice,nmet,NORM_1,metal_content_in_simcell,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Wait for ALL RANKS to complete coupling the kMC to PF before erasing the 'toucouple' file
  if(isroot) call system('rm -f tocouple')
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before going back to the PF routines

end subroutine kMC_filmdissolve
