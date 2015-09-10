!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module fields
  implicit none
  save

  integer :: psx_g,psy_g,psz_g ! Number of PF gridpoints in 3 directions
  integer :: psx,psy,psz       ! Number of PF gridpoints in 3 directions

  ! Local phase fields
  real*8, dimension(:,:,:), allocatable :: met,mkw,pht,pyr,env ! Local PF grid for metal, pyrrhotite, pyrite and environment
  real*8, dimension(:,:,:), allocatable :: mu                  ! Local mu_S grid
  real*8, dimension(:,:,:), allocatable :: opyr                ! Local orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: voids               ! Local void fraction
  real*8, dimension(:,:,:), allocatable :: elpot               ! Local electrical potential

  ! Global phase fields
  real*8, dimension(:,:,:), allocatable :: met_g, mkw_g, pht_g, pyr_g, env_g ! Global PF grid for metal, pyrrhotite, pyrite and environment
  real*8, dimension(:,:,:), allocatable :: mu_g                              ! Global mu_S grid
  real*8, dimension(:,:,:), allocatable :: opyr_g                            ! Global orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: elpot_g                           ! Global electrical potential

  real*8, dimension(:,:,:), allocatable :: dmet_dt, dmkw_dt, dpht_dt, dpyr_dt, denv_dt, dmu_dt
  real*8, dimension(:,:,:), allocatable :: newmet, newmkw, newpht, newpyr, newenv


      type context
#include <finclude/petscdmdef.h>
         DM lattval, exlattval
         integer :: startx,starty,startz
         integer :: widthx,widthy,widthz
      end type context

end module fields

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module commondata
  implicit none
  save

  ! MPI parameters
  logical :: isroot
  integer :: rank, procs

  ! Input parameters
  character*1 :: isrestart     ! Is the calculation a restarted one?
  integer :: nomc              ! Number of PF iterations
  integer :: T                 ! Temperature of the simulation box
  real*8 :: pH_in              ! Scalar pH input
  integer :: noimg             ! Number of output files
  real*8 :: metal_potential    ! Electric potential of the metal
  logical :: include_dissolve  ! Include film dissolution
  logical :: include_electro   ! Include potential distribution


  ! Simulation parameters
  real*8, parameter :: dpf = 5E-9       ! Phase-field grid size
  real*8 :: dt                          ! Timestep for PF evolution
  real*8 :: avg_mu_env                  ! Sulfur chemical potential in the environment
  real*8 :: sulfidation_rate            ! Sulfidation rate / Film growth rate in m/s
  integer, parameter :: ghost_width = 2 ! Number of ghost nodes (in the z direction)
  integer :: swap_freq_pf = 5           ! Frequency with which MPI swaps are conducted for PF solving
  integer :: swap_freq_kmc = 10         ! Frequency with which MPI swaps are conducted for KMC solving
  integer :: freq_scale = 1750000000    ! KMC information is transferred every freq_scale steps
  integer :: kg_scale = 5               ! Number of KMC grid points per PF grid

  integer, dimension(:), allocatable :: seed  ! Seed for PRNG

  integer, parameter :: nfields = 9
  integer, parameter :: nphases = 5
  integer, parameter :: nmet = 0
  integer, parameter :: nmkw = 1
  integer, parameter :: npht = 2
  integer, parameter :: npyr = 3
  integer, parameter :: nenv = 4
  integer, parameter :: nmus = 5
  integer, parameter :: npH = 6
  integer, parameter :: nang = 7
  integer, parameter :: npot = 8

end module commondata

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module thermo_constants
  implicit none
  save

  real*8 :: mus_met_mkw_eqb, mus_mkw_pht_eqb, mus_pht_pyr_eqb ! mu_S boundary between phases
  real*8 :: rho_met, rho_mkw, rho_pht, rho_pyr, rho_env       ! Sulfur density
  real*8, parameter :: R = 8.3144621

  ! Surface energy sigmas
  real*8, parameter :: sigma_mkw_env= 1E-12, sigma_env_mkw= 1E-12
  real*8, parameter :: sigma_mkw_met= 1E-12, sigma_met_mkw= 1E-12
  real*8, parameter :: sigma_mkw_pht= 1E-12, sigma_pht_mkw= 1E-12
  real*8, parameter :: sigma_mkw_pyr_0 = 1E-12, sigma_pyr_mkw_0 = 1E-12
  real*8, parameter :: sigma_pht_env= 1E-12, sigma_env_pht= 1E-12
  real*8, parameter :: sigma_pht_met= 1E-12, sigma_met_pht= 1E-12
  real*8, parameter :: sigma_pht_pyr_0 = 1E-12, sigma_pyr_pht_0 = 1E-12
  real*8, parameter :: sigma_met_env= 1E-12, sigma_env_met= 1E-12
  real*8, parameter :: sigma_met_pyr_0 = 1E-12, sigma_pyr_met_0 = 1E-12
  real*8, parameter :: sigma_env_pyr_0 = 1E-12, sigma_pyr_env_0 = 1E-12

  ! Field mobilities
  real*8, parameter :: M_pht_met = 3.00E-08, M_met_pht = 3.00E-08
  real*8, parameter :: M_mkw_met = 1.40E-06, M_met_mkw = 1.40E-06
  real*8, parameter :: M_met_pyr = 3.00E-08, M_pyr_met = 3.00E-08
  real*8, parameter :: M_pht_pyr = 5.00E-07, M_pyr_pht = 5.00E-07
  real*8, parameter :: M_mkw_pyr = 3.00E-08, M_pyr_mkw = 3.00E-08
  real*8, parameter :: M_pht_mkw = 6.00E-07, M_mkw_pht = 6.00E-07
  real*8, parameter :: M_pht_env = 4.00E-15, M_env_pht = 4.00E-15
  real*8, parameter :: M_mkw_env = 4.00E-15, M_env_mkw = 4.00E-15
  real*8, parameter :: M_met_env = 4.00E-15, M_env_met = 4.00E-15
  real*8, parameter :: M_env_pyr = 4.00E-15, M_pyr_env = 4.00E-15

  ! Double well potential heights
  real*8, parameter :: double_well_barrier = 100.0d0 ! in J/mol
  real*8, parameter :: hill = (16.0d0/3.0d0)*double_well_barrier
  real*8 :: hill_met_mkw, hill_met_pht, hill_met_pyr, hill_met_env
  real*8 :: hill_mkw_met, hill_mkw_pht, hill_mkw_pyr, hill_mkw_env
  real*8 :: hill_pht_met, hill_pht_mkw, hill_pht_pyr, hill_pht_env
  real*8 :: hill_pyr_met, hill_pyr_mkw, hill_pyr_pht, hill_pyr_env
  real*8 :: hill_env_met, hill_env_mkw, hill_env_pht, hill_env_pyr

  ! Field mobilities
  real*8, allocatable :: Mob_pf(:,:)

  ! Phase surface energies
  real*8, allocatable :: sigma(:,:)
  real*8, parameter :: sigma_pyr_0 = 1E-12

  ! Relative phase stabilities
  real*8, allocatable :: w_pf(:)

  ! Grain boundary stability
  real*8 :: gb_S = 0.0d0

  ! Relative permittivities of FeS phases
  real*8 :: epsilon0
  real*8, allocatable :: permittivity(:)

  ! Derivative of sulfur concentration with chemical potential
  real*8, allocatable :: drho_dmu(:)

  ! Sulfur density in different phases
  real*8, allocatable :: rhoS(:)

  ! Collected sulfidation rates for different FeS phases and environments
  real*8, allocatable :: sulf_rate_gas(:)
  real*8, allocatable :: sulf_rate_liq(:)
  real*8, allocatable :: sulf_rate(:)

end module thermo_constants

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module diffusion_constants
  implicit none
  save

  real*8 :: D_S_met, D_S_mkw, D_S_pht, D_S_pyr, D_S_env
  real*8 :: D_Fe_met, D_Fe_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_env
  real*8 :: D_H_env, D_H_met, D_H_pht, D_H_pyr

end module diffusion_constants

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module gradients
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: delypyr,delzpyr

end module gradients

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!
