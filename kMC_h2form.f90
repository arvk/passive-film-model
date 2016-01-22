subroutine kMC_h2form(iter,simstate,random_context)
  use, intrinsic :: iso_c_binding
  use commondata
  use fields
  use thermo_constants
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
  logical :: already_run_once
  type(context) simstate
  PetscInt, intent(in) :: iter  ! Iteration count
  PetscInt, dimension(psx_g,psy_g) :: coarse_h2_evolved
  PetscInt :: coarsex,coarsey
  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscInt :: floor
  PetscScalar :: max, min
  PetscRandom, intent(inout) :: random_context !! Context to seed and generate random numbers
  PetscReal :: random_number  !! Pseudo random number generated from a PETSc context

  PetscScalar :: mkw_modulus  !! Elastic modulus of mackinawite. Source: Own work
  PetscScalar :: mkw_surf_energy !! Surface energy of mackinawite Source: Own work
  PetscScalar :: blist_size !! Lateral size of the H2 blister
  PetscScalar :: blist_rad_curv !! Radius of curvature of the blister
  PetscScalar :: mkw_thickness !! Thickness of the mackinawite film that undergoes blistering
  PetscScalar :: pht_thickness !! Thickness of the pyrrhotite film that undergoes blistering
  PetscScalar :: pyr_thickness !! Thickness of the pyrite film that undergoes blistering
  PetscScalar :: film_thickness !! Total thickness of the composite FeS_x film that undergoes blistering
  PetscScalar :: threshold_P !! Critical pressure for stable blister formation
  PetscScalar :: delamination_height !! Distance between two (001) planes of mackinawite before delamination occurs
  PetscScalar :: vol_spher_cap !! Volume of the blister. Bilster is assumed to be a spherical cap
  PetscScalar :: threshold_H2_molecules !! Number of H2 molecules inside a stable blister
  PetscScalar :: Pi = 3.14159265d0 !! \(\pi\)

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

  call PetscRandomGetValueReal(random_context,random_number,ierr)

  if(isroot)then
     call system('rm -f metfilm.spparksscript')
     open(unit = 666, file = 'metfilm.spparksscript', status = 'new')
     write(666,*) 'seed ', floor(100000.0d0*random_number)+1
     write(666,*) 'app_style metfilm'
     write(666,*) 'dimension 2'
     write(666,*) 'boundary p p p'
     write(666,*) 'lattice mkw 1.0'
     write(666,'(A,F16.8)') 'temperature ', 8.61733034E-5*T
     write(666,*) 'region cell block 0 ',psx_g*kg_scale,' 0 ',psy_g*kg_scale,' -0.5 0.5'
     write(666,*) 'create_box cell'
     write(666,*) '# BASIS : 1=Fe, 2=Fe, 3=S_up, 4=S_down'
     write(666,*) '# SITE IDS (I2) : 1=Fe, 2=Sup-H, 3=Fe-H2, 4=V_Fe, 5=S_up, 6=S_down'
     write(666,*) 'create_sites box value i2 2 basis 1 1 basis 2 1 basis 3 5 basis 4 5'
     write(666,*) 'set i1 value 2'
     write(666,*) 'set i2 value 2 if i2 = 5 fraction 0.1'
     write(666,*) 'set i2 value 4 if i2 = 1 fraction 0.01'
     write(666,*) 'solve_style tree'
     write(666,*) 'sector yes'
     write(666,*) 'event 2 oct oct fe sup 0.32355940 fe h  # Adsorption'
     write(666,*) 'event 2 oct oct fe h 0.10 fe sup  # Desorption'
     write(666,*) 'event 2 oct oct h sup 0.47584796 sup h  # H glide along plane'
     write(666,*) 'event 3 oct oct oct fe h h 0.68543162 fe sup sup'
     write(666,*) 'event 3 oct oct oct fe sup sup 0.58708402 fe h h'
     write(666,*) 'event 3 oct oct oct vac h h 0.73076749 vac sup sup'
     write(666,*) 'event 3 oct oct oct vac sup sup 0.79569354 vac h h'
     write(666,*) 'diag_style metfilm stats yes list events d1 d2 d3 t1 t2 t3 t4'
     write(666,'(A,F16.8,A)') 'dump raw_metfilm_output text ', dt*kmc_freq, ' raw_metfilm_output id x y i1 i2 i3'
     write(666,'(A,F16.8)') 'dump_modify raw_metfilm_output delay ', dt*kmc_freq
     if (already_run_once) write(666,*) 'read_sites input.metfilm'
     write(666,'(A,F16.8,A)') 'run ', dt*kmc_freq, ' pre yes post no'
     close(666)
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  call spparks_open(myargc,C_LOC(myargv),C_LOC(MPI_COMM_WORLD),C_LOC(myspparks))
  call spparks_file(myspparks,'metfilm.spparksscript'//C_NULL_CHAR)
  call spparks_close(myspparks)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  if(isroot)then
     call system('rm -f input.metfilm log.spparks')
     ! READ RAW SPPARKS OUTPUT AND WRITE FORTRAN-READABLE OUTPUT
     write(kmc_numel_string,'(I24)') psx_g*psy_g*kg_scale*kg_scale*4
     call system('tail -n '//trim(kmc_numel_string)//' raw_metfilm_output | awk ''{printf " %7.7i %4.4i %4.4i %3.3i %3.3i %5.5i \n", $1,$2,$3,$4,$5,$6}'' > metfilm_spparks_output.template')
     call system('rm -f raw_metfilm_output')
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  if(isroot)then
     !  READ FORTRAN-READABLE INPUT AND STORE IN ARRAY FOR COUPLING

     call system('rm -f input.metfilm')
     open(unit = 666, file = 'input.metfilm', status = 'new')
     open (unit = 667, file = 'metfilm_spparks_output.template', status = 'old')

     write(666,*) 'Testing'
     write(666,*) '2 dimension'
     write(666,*) '0 ', psx_g*kg_scale ,' xlo xhi'
     write(666,*) '0 ', psy_g*kg_scale ,' ylo yhi'
     write(666,*) '-0.5 0.5 zlo zhi'
     write(666,*) psx_g*psy_g*kg_scale*kg_scale*4, ' sites'
     write(666,*) ' '
     write(666,*) 'Values'
     write(666,*) ' '

     do x = 1,psx_g*psy_g*kg_scale*kg_scale*4
        read(667,'(I8, I5, I5, I4, I4, I6)') site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
        write(666,'(I8,I4,I4,I6)') site_id, i1, i2, no_h2_already_evolved
     end do

     close(667)
     close(666)

  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning SPPARKS functions

  coarse_h2_evolved = 0
  open (unit = 667, file = 'metfilm_spparks_output.template', status = 'old')
  do x = 1,psx_g*psy_g*kg_scale*kg_scale*4
     read(667,'(I8, I5, I5, I4, I4, I6)') site_id, partial_x, partial_y, i1, i2, no_h2_already_evolved
     coarsex = floor(1.0d0*partial_x/kg_scale)
     coarsey = floor(1.0d0*partial_y/kg_scale)
     coarse_h2_evolved(coarsex+1,coarsey+1) = coarse_h2_evolved(coarsex+1,coarsey+1) + no_h2_already_evolved
  end do
  close(667)


  ! Calculate scaling factors

  ! Equations to use
  ! ================
  ! Threshold Pressure = 4*sigma*(t^2)/(3*R_c^2) + 2*gamma/R_c
  !    where, sigma = Elastic modulus, t = thickness of film, R_c = blister radius, gamma = surface energy

  ! Number of H2 molecules = N_a * PV/RT
  !    where, N_a = Avogadro number, P = Pressure, V = Volume of blister, R = gas constant, T = temperature

  ! REFERENCE: In situ study of the initiation of hydrogen bubbles at the aluminium metal/oxide interface, Nature Materials 14, 899-903 (2015) DOI: 10.1038/nmat4336
  ! REFERENCE: Hydrogen bubbles in metals, J.B. Condon, T. Schober, Journal of Nuclear Materials, Volume 207, December 1993, Pages 1-24, DOI: 10.1016/0022-3115(93)90244-S

  ! Calculate parameters
  mkw_modulus = 15E9    ! in Pa
  mkw_surf_energy = 0.1 ! in J/m^2

  blist_size = 25E-9    ! in m
  ! Size of incipient corrosion blister. REFERENCE: In situ study of the initiation of hydrogen bubbles at the aluminium metal/oxide interface, Nature Materials 14, 899-903 (2015) DOI: 10.1038/nmat4336

  call VecStrideNorm(simstate%slice,nmkw,NORM_1,mkw_thickness,ierr) ; mkw_thickness = (mkw_thickness/(psx_g*psy_g))*dpf ! in m
  call VecStrideNorm(simstate%slice,npht,NORM_1,pht_thickness,ierr) ; pht_thickness = (pht_thickness/(psx_g*psy_g))*dpf ! in m
  call VecStrideNorm(simstate%slice,npyr,NORM_1,pyr_thickness,ierr) ; pyr_thickness = (pyr_thickness/(psx_g*psy_g))*dpf ! in m
  film_thickness = mkw_thickness + pht_thickness + pyr_thickness
  delamination_height = 5.5E-10 ! in m
  blist_rad_curv = ((blist_size*blist_size)-(delamination_height*delamination_height))/(2*delamination_height) ! in m
  vol_spher_cap = Pi * delamination_height*((3*blist_size*blist_size)+(delamination_height*delamination_height))/6 ! in m

  threshold_P = ((4*mkw_modulus*film_thickness*film_thickness)/(3*blist_rad_curv*blist_rad_curv)) + ((2*mkw_surf_energy)/blist_rad_curv)
  threshold_H2_molecules = 6.022E23 * threshold_P * vol_spher_cap / (R*T)

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  do x = simstate%startx , simstate%startx + simstate%widthx-1
     do y = simstate%starty , simstate%starty + simstate%widthy-1
        do z = simstate%startz , simstate%startz + simstate%widthz-2
           if ((coarse_h2_evolved(x+1,y+1)/threshold_H2_molecules).gt.0.05d0) then
              if (((statepointer(nmkw,x,y,z).lt.0.1d0)) .and. ((statepointer(nmkw,x,y,z+1).gt.0.1d0))) then
                 statepointer(nvoi,x,y,z) = statepointer(nvoi,x,y,z) - (coarse_h2_evolved(x+1,y+1)/threshold_H2_molecules)
                 statepointer(nvoi,x,y,z) = min(max(statepointer(nvoi,x,y,z),0.0d0),1.0d0)
              end if
           end if
        end do
     end do
  end do
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before file IO
  if(isroot)then
     call system('rm -f metfilm_spparks_output.template')
  end if
  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before returning to PF routines

end subroutine kMC_h2form
