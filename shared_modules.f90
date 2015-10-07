!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module fields
  implicit none
  save
#include <finclude/petscsysdef.h>


  PetscInt :: psx_g,psy_g,psz_g !! Number of gridpoints in the simulation cell along 3 directions

  ! Global phase fields
  PetscScalar, dimension(:,:,:), allocatable :: met_g, mkw_g, pht_g, pyr_g, env_g !! Global PF grid for metal, pyrrhotite, pyrite and environment
  PetscScalar, dimension(:,:,:), allocatable :: mu_g                              !! Global mu_S grid
  PetscScalar, dimension(:,:,:), allocatable :: opyr_g                            !! Global orientation field for pyrite
  PetscScalar, dimension(:,:,:), allocatable :: elpot_g                           !! Global electrical potential

      type context
#include <finclude/petscdmdef.h>
         DM lattval
         Vec slice, exslice   !! All field variables that define the system state
         PetscInt :: startx,starty,startz !! Coordinates that define the bottom-left corner of the local part of the simulation cell
         PetscInt :: widthx,widthy,widthz !! Size of the local part of the simulation cell
      end type context

end module fields

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module commondata
  implicit none
#include <finclude/petscsysdef.h>
  save

  ! MPI parameters
  PetscBool :: isroot     !! Is this processor 0?
  PetscInt :: rank       !! What is the rank of the current processor?
  PetscInt :: procs      !! How many processors is the program running on?

  ! Input parameters
  character*1 :: isrestart     !! Is the calculation a restarted one?
  PetscInt :: nomc              !! Number of PF iterations
  PetscInt :: T                 !! Temperature of the simulation box
  PetscScalar :: pH_in              !! Scalar pH input
  PetscInt :: num_images             !! Number of output files
  PetscScalar :: metal_potential    !! Electric potential of the metal
  PetscBool :: include_dissolve  !! Include film dissolution
  PetscBool :: include_electro   !! Include potential distribution


  ! Simulation parameters
  PetscScalar, parameter :: dpf = 5E-9       !! Phase-field grid size
  PetscScalar :: dt                          !! Timestep for PF evolution
  PetscScalar :: avg_mu_env                  !! Sulfur chemical potential in the environment
  PetscScalar :: sulfidation_rate            !! Sulfidation rate / Film growth rate in m/s
  PetscInt :: kg_scale = 5               !! Number of KMC grid points per PF grid
  PetscInt :: kmc_freq = 100             !! Frequency with which kMC calculations are done
  PetscInt :: stat_freq = 200            !! Simulation stats are reported once every stat_freq iterations

  PetscInt, dimension(:), allocatable :: seed  ! Seed for PRNG

  PetscInt, parameter :: nfields = 10 !! Number of field variables at each grid point
  PetscInt, parameter :: nphases = 5 !! Number of FeS phases
  PetscInt, parameter :: nmet = 0 !! Index for the metal phase fraction
  PetscInt, parameter :: nmkw = 1 !! Index for the mackinawite phase fraction
  PetscInt, parameter :: npht = 2 !! Index for the pyrrhotite phase fraction
  PetscInt, parameter :: npyr = 3 !! Index for the pyrite phase fraction
  PetscInt, parameter :: nenv = 4 !! Index for the environment phase fraction
  PetscInt, parameter :: nmus = 5 !! Index for the sulfur chemical potential field
  PetscInt, parameter :: npH = 6 !! Index for the pH field
  PetscInt, parameter :: nang = 7 !! Index for the pyrite grain shape and orientation field
  PetscInt, parameter :: npot = 8 !! Index for the electrical potential field
  PetscInt, parameter :: nvoi = 9 !! Index for the void field

end module commondata

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module thermo_constants
  implicit none
#include <finclude/petscsysdef.h>
  save

  PetscScalar :: mus_met_mkw_eqb, mus_mkw_pht_eqb, mus_pht_pyr_eqb !! mu_S boundary between phases
  PetscScalar :: rho_met, rho_mkw, rho_pht, rho_pyr, rho_env       !! Sulfur density
  PetscScalar, parameter :: R = 8.3144621

  ! Double well potential heights
  PetscScalar, parameter :: double_well_barrier = 20.0d0 ! in J/mol
  PetscScalar, parameter :: hill = (16.0d0/3.0d0)*double_well_barrier
  !! Magnitude of the barrier in the double-well potential

  PetscScalar, allocatable :: Mob_pf(:,:) !! Field mobilities

  PetscScalar, allocatable :: sigma(:,:)  !! Energy of interfaces between phases
  PetscScalar, parameter :: sigma_pyr_0 = 1E-16 !! Energy of interfaces between pyrite and other phases

  PetscScalar, allocatable :: w_pf(:)  !! Free energy of different FeS phases

  PetscScalar :: gb_S = 0.0d0 !! Stability of boundary between different pyrite grains

  PetscScalar :: epsilon0  !! Vacuum permittivity
  PetscScalar, allocatable :: permittivity(:)  !! Relative permittivity

  PetscScalar, allocatable :: drho_dmu(:)  !! Derivative of sulfur concentration with chemical potential

  PetscScalar, allocatable :: rhoS(:)  !! Sulfur density in different FeS phases

  PetscScalar, allocatable :: sulf_rate_gas(:) !! Collected sulfidation rates for different FeS phases in gaseous environments
  PetscScalar, allocatable :: sulf_rate_liq(:) !! Collected sulfidation rates for different FeS phases in gaseous environments
  PetscScalar, allocatable :: sulf_rate(:)

end module thermo_constants

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module diffusion_constants
  implicit none
#include <finclude/petscsysdef.h>
  save

  PetscScalar :: D_S_met, D_S_mkw, D_S_pht, D_S_pyr, D_S_env  !! Diffusivity of sulfur in different FeS phases
  PetscScalar :: D_Fe_met, D_Fe_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_env  !! Diffusivity of iron in different FeS phases
  PetscScalar :: D_H_env, D_H_met, D_H_pht, D_H_pyr  !! Diffusivity of hydrogen in different FeS phases

end module diffusion_constants

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!
