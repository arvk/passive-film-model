!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module fields
  implicit none
  save

  integer :: psx_g,psy_g,psz_g ! Number of PF gridpoints in 3 directions
  integer :: psx,psy,psz       ! Number of PF gridpoints in 3 directions

  ! Local phase fields
  real*8, dimension(:,:,:), allocatable :: met,mkw,pht,pyr,env ! Local PF grid for metal, pyrrhotite, pyrite and environment
  real*8, dimension(:,:,:), allocatable :: mu                  ! Local mu_S grid
  real*8, dimension(:,:,:), allocatable :: pH                  ! Local pH grid
  real*8, dimension(:,:,:), allocatable :: opyr                ! Local orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: voids               ! Local void fraction
  real*8, dimension(:,:,:), allocatable :: elpot               ! Local electrical potential

  ! Global phase fields
  real*8, dimension(:,:,:), allocatable :: met_g, mkw_g, pht_g, pyr_g, env_g ! Global PF grid for metal, pyrrhotite, pyrite and environment
  real*8, dimension(:,:,:), allocatable :: mu_g                              ! Global mu_S grid
  real*8, dimension(:,:,:), allocatable :: pH_g                              ! Global pH grid
  real*8, dimension(:,:,:), allocatable :: opyr_g                            ! Global orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: elpot_g                           ! Global electrical potential

  real*8, dimension(:,:,:), allocatable :: dmet_dt, dmkw_dt, dpht_dt, dpyr_dt, denv_dt, dmu_dt
  real*8, dimension(:,:,:), allocatable :: newmet, newmkw, newpht, newpyr, newenv

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
  real*8 :: metal_amount                ! Amount of metal in the simulation cell
  real*8 :: sulfidation_rate            ! Sulfidation rate / Film growth rate in m/s
  integer, parameter :: ghost_width = 2 ! Number of ghost nodes (in the z direction)
  integer :: swap_freq_pf = 5           ! Frequency with which MPI swaps are conducted for PF solving
  integer :: swap_freq_kmc = 10         ! Frequency with which MPI swaps are conducted for KMC solving
  integer :: freq_scale = 1750000000    ! KMC information is transferred every freq_scale steps
  integer :: kg_scale = 5               ! Number of KMC grid points per PF grid

  integer, dimension(:), allocatable :: seed  ! Seed for PRNG

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
  real*8 :: hill_met_mkw, hill_met_pht, hill_met_pyr, hill_met_env
  real*8 :: hill_mkw_met, hill_mkw_pht, hill_mkw_pyr, hill_mkw_env
  real*8 :: hill_pht_met, hill_pht_mkw, hill_pht_pyr, hill_pht_env
  real*8 :: hill_pyr_met, hill_pyr_mkw, hill_pyr_pht, hill_pyr_env
  real*8 :: hill_env_met, hill_env_mkw, hill_env_pht, hill_env_pyr

end module thermo_constants

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module kmc_data
  implicit none
  save

  integer :: ksx,ksy           ! Number of local kMC gridpoints
  integer :: ksx_g,ksy_g       ! Number of global kMC gridpoints

  integer, dimension(:,:), allocatable :: kg      ! Global kMC grid
  integer, dimension(:,:), allocatable :: kg_g    ! Local kMC grid
  integer, dimension(:,:), allocatable :: kg_recv ! Empty array to receive MPI requests

  type :: prol
     integer :: fx,fy,tx,ty,from,to
     real*8 :: prob
  end type prol

  type(prol), dimension(:,:), allocatable :: plist

!!!!------------------------------------------------------!!!!
!!!!              KMC Process probabilities               !!!!
  ! Global                                                !!!!
  real*8, dimension(:,:), allocatable :: vfe_f_g,vfe_a_g  !!!!
  real*8, dimension(:,:), allocatable :: vs_f_g,vs_a_g    !!!!
  real*8, dimension(:,:), allocatable :: fes_diss_g       !!!!
  real*8, dimension(:,:), allocatable :: v_diff_g         !!!!
  ! Local                                                 !!!!
  real*8, dimension(:,:), allocatable :: vfe_f,vfe_a      !!!!
  real*8, dimension(:,:), allocatable :: vs_f,vs_a        !!!!
  real*8, dimension(:,:), allocatable :: fes_diss         !!!!
  real*8, dimension(:,:), allocatable :: v_diff           !!!!
!!!!------------------------------------------------------!!!!

end module kmc_data

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module diffusion_constants
  implicit none
  save

  real*8 :: D_S_met, D_S_mkw, D_S_pht, D_S_pyr, D_S_env
  real*8 :: D_Fe_met, D_Fe_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_env
  real*8 :: D_H_met, D_H_mkw, D_H_pht, D_H_pyr, D_H_env

end module diffusion_constants

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

module gradients
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: delypyr,delzpyr

end module gradients

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!
