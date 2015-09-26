!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module fields
  implicit none
  save

  integer :: psx_g,psy_g,psz_g !! Number of gridpoints in the simulation cell along 3 directions

  ! Global phase fields
  real*8, dimension(:,:,:), allocatable :: met_g, mkw_g, pht_g, pyr_g, env_g !! Global PF grid for metal, pyrrhotite, pyrite and environment
  real*8, dimension(:,:,:), allocatable :: mu_g                              !! Global mu_S grid
  real*8, dimension(:,:,:), allocatable :: opyr_g                            !! Global orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: elpot_g                           !! Global electrical potential

      type context
#include <finclude/petscdmdef.h>
         DM lattval
         Vec slice, exslice   !! All field variables that define the system state
         integer :: startx,starty,startz !! Coordinates that define the bottom-left corner of the local part of the simulation cell
         integer :: widthx,widthy,widthz !! Size of the local part of the simulation cell
      end type context

end module fields

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module commondata
  implicit none
  save

  ! MPI parameters
  logical :: isroot     !! Is this processor 0?
  integer :: rank       !! What is the rank of the current processor?
  integer :: procs      !! How many processors is the program running on?

  ! Input parameters
  character*1 :: isrestart     !! Is the calculation a restarted one?
  integer :: nomc              !! Number of PF iterations
  integer :: T                 !! Temperature of the simulation box
  real*8 :: pH_in              !! Scalar pH input
  integer :: noimg             !! Number of output files
  real*8 :: metal_potential    !! Electric potential of the metal
  logical :: include_dissolve  !! Include film dissolution
  logical :: include_electro   !! Include potential distribution


  ! Simulation parameters
  real*8, parameter :: dpf = 5E-9       !! Phase-field grid size
  real*8 :: dt                          !! Timestep for PF evolution
  real*8 :: avg_mu_env                  !! Sulfur chemical potential in the environment
  real*8 :: sulfidation_rate            !! Sulfidation rate / Film growth rate in m/s
  integer :: kg_scale = 5               !! Number of KMC grid points per PF grid
  integer :: kmc_freq = 100             !! Frequency with which kMC calculations are done

  integer, dimension(:), allocatable :: seed  ! Seed for PRNG

  integer, parameter :: nfields = 10 !! Number of field variables at each grid point
  integer, parameter :: nphases = 5 !! Number of FeS phases
  integer, parameter :: nmet = 0 !! Index for the metal phase fraction
  integer, parameter :: nmkw = 1 !! Index for the mackinawite phase fraction
  integer, parameter :: npht = 2 !! Index for the pyrrhotite phase fraction
  integer, parameter :: npyr = 3 !! Index for the pyrite phase fraction
  integer, parameter :: nenv = 4 !! Index for the environment phase fraction
  integer, parameter :: nmus = 5 !! Index for the sulfur chemical potential field
  integer, parameter :: npH = 6 !! Index for the pH field
  integer, parameter :: nang = 7 !! Index for the pyrite grain shape and orientation field
  integer, parameter :: npot = 8 !! Index for the electrical potential field
  integer, parameter :: nvoi = 9 !! Index for the void field

end module commondata

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module thermo_constants
  implicit none
  save

  real*8 :: mus_met_mkw_eqb, mus_mkw_pht_eqb, mus_pht_pyr_eqb !! mu_S boundary between phases
  real*8 :: rho_met, rho_mkw, rho_pht, rho_pyr, rho_env       !! Sulfur density
  real*8, parameter :: R = 8.3144621

  ! Double well potential heights
  real*8, parameter :: double_well_barrier = 100.0d0 ! in J/mol
  real*8, parameter :: hill = (16.0d0/3.0d0)*double_well_barrier
  !! Magnitude of the barrier in the double-well potential

  real*8, allocatable :: Mob_pf(:,:) !! Field mobilities

  real*8, allocatable :: sigma(:,:)  !! Energy of interfaces between phases
  real*8, parameter :: sigma_pyr_0 = 1E-12 !! Energy of interfaces between pyrite and other phases

  real*8, allocatable :: w_pf(:)  !! Free energy of different FeS phases

  real*8 :: gb_S = 0.0d0 !! Stability of boundary between different pyrite grains

  real*8 :: epsilon0  !! Vacuum permittivity
  real*8, allocatable :: permittivity(:)  !! Relative permittivity

  real*8, allocatable :: drho_dmu(:)  !! Derivative of sulfur concentration with chemical potential

  real*8, allocatable :: rhoS(:)  !! Sulfur density in different FeS phases

  real*8, allocatable :: sulf_rate_gas(:) !! Collected sulfidation rates for different FeS phases in gaseous environments
  real*8, allocatable :: sulf_rate_liq(:) !! Collected sulfidation rates for different FeS phases in gaseous environments
  real*8, allocatable :: sulf_rate(:)

end module thermo_constants

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

module diffusion_constants
  implicit none
  save

  real*8 :: D_S_met, D_S_mkw, D_S_pht, D_S_pyr, D_S_env  !! Diffusivity of sulfur in different FeS phases
  real*8 :: D_Fe_met, D_Fe_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_env  !! Diffusivity of iron in different FeS phases
  real*8 :: D_H_env, D_H_met, D_H_pht, D_H_pyr  !! Diffusivity of hydrogen in different FeS phases

end module diffusion_constants

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!
