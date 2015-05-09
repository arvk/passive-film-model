!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module fields
  implicit none
  save

  !! Define the extent of phase fields
  integer :: psx_g,psy_g,psz_g ! Number of PF gridpoints in 3 directions
  integer :: psx,psy,psz       ! Number of PF gridpoints in 3 directions

  !! Define local phase fields
  real*8, dimension(:,:,:), allocatable :: met,mkw,pht,pyr,env ! Local PF grid for metal, pyrrhotite, pyrite and environment respectively
  real*8, dimension(:,:,:), allocatable :: mu  ! Local mu_S grid
  real*8, dimension(:,:,:), allocatable :: opyr  ! Local orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: voids ! Local void fraction
  real*8, dimension(:,:,:), allocatable :: elpot ! Local electrical potential

  !! Define global phase fields
  real*8, dimension(:,:,:), allocatable :: met_g, mkw_g, pht_g, pyr_g, env_g ! Local PF grid for metal, pyrrhotite, pyrite and environment respectively
  real*8, dimension(:,:,:), allocatable :: mu_g ! Local mu_S grid
  real*8, dimension(:,:,:), allocatable :: opyr_g  ! Global orientation field for pyrite
  real*8, dimension(:,:,:), allocatable :: elpot_g ! Global electrical potential

  real*8, dimension(:,:,:), allocatable :: dmet_dt, dmkw_dt, dpht_dt, dpyr_dt, denv_dt, dmu_dt
  real*8, dimension(:,:,:), allocatable :: newmet, newmkw, newpht, newpyr, newenv

end module fields


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module commondata
  implicit none
  save

  ! MPI PARAMETERS
  logical :: isroot
  integer :: rank, procs

  ! INPUT PARAMETERS
  character*1 :: isrestart     ! Is the calculation a restarted one?
  integer :: nomc              ! Number of PF iterations
  integer :: T                 ! Temperature of the simulation box
  real*8 :: pH_in              ! Scalar pH input
  integer :: noimg             ! Number of output files
  real*8 :: metal_potential    ! Electric potential of the metal
  logical :: include_dissolve  ! Include film dissolution
  logical :: include_electro   ! Include potential distribution


  ! INTERNAL PARAMETERS
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

end module commondata


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module thermo_constants
  implicit none
  save

  real*8 :: mus_met_mkw_eqb, mus_mkw_pht_eqb, mus_pht_pyr_eqb
  real*8 :: rho_met, rho_mkw, rho_pht, rho_pyr, rho_env    !! Sulfur density in different phases
  real*8, parameter :: R = 8.3144621   

  !! Surface energy sigmas
  real*8 :: sigma_mkw_env= 1E-7, sigma_env_mkw= 1E-7
  real*8 :: sigma_mkw_met= 1E-7, sigma_met_mkw= 1E-7
  real*8 :: sigma_mkw_pht= 1E-7, sigma_pht_mkw= 1E-7
  real*8 :: sigma_mkw_pyr_0 = 1E-7, sigma_pyr_mkw_0 = 1E-7
  real*8 :: sigma_pht_env= 1E-7, sigma_env_pht= 1E-7
  real*8 :: sigma_pht_met= 1E-7, sigma_met_pht= 1E-7
  real*8 :: sigma_pht_pyr_0 = 1E-7, sigma_pyr_pht_0 = 1E-7
  real*8 :: sigma_met_env= 1E-7, sigma_env_met= 1E-7
  real*8 :: sigma_met_pyr_0 = 1E-7, sigma_pyr_met_0 = 1E-7
  real*8 :: sigma_env_pyr_0 = 1E-7, sigma_pyr_env_0 = 1E-7

  !! Field mobilities
  real*8 :: M_pht_met = 1.0E-12, M_met_pht = 1.0E-12
  real*8 :: M_mkw_met = 1.0E-12, M_met_mkw = 1.0E-12
  real*8 :: M_met_pyr = 1.0E-12, M_pyr_met = 1.0E-12
  real*8 :: M_pht_pyr = 1.0E-12, M_pyr_pht = 1.0E-12
  real*8 :: M_mkw_pyr = 1.0E-12, M_pyr_mkw = 1.0E-12
  real*8 :: M_pht_mkw = 1.0E-12, M_mkw_pht = 1.0E-12
  real*8 :: M_pht_env = 1.0E-15, M_env_pht = 1.0E-15
  real*8 :: M_mkw_env = 1.0E-15, M_env_mkw = 1.0E-15
  real*8 :: M_met_env = 1.0E-15, M_env_met = 1.0E-15
  real*8 :: M_env_pyr = 1.0E-15, M_pyr_env = 1.0E-15

  !! Double well potential heights
  real*8 :: double_well_barrier = 1000.0d0 ! in J/mol
  real*8 :: hill_met_mkw, hill_met_pht, hill_met_pyr, hill_met_env
  real*8 :: hill_mkw_met, hill_mkw_pht, hill_mkw_pyr, hill_mkw_env
  real*8 :: hill_pht_met, hill_pht_mkw, hill_pht_pyr, hill_pht_env
  real*8 :: hill_pyr_met, hill_pyr_mkw, hill_pyr_pht, hill_pyr_env
  real*8 :: hill_env_met, hill_env_mkw, hill_env_pht, hill_env_pyr

end module thermo_constants


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module kmc_data
  implicit none
  save

  integer :: ksx,ksy           ! Number of local kMC gridpoints
  integer :: ksx_g,ksy_g       ! Number of global kMC gridpoints


  integer, dimension(:,:), allocatable :: kg      ! Global kMC grid
  integer, dimension(:,:), allocatable :: kg_g    ! Local kMC grid
  integer, dimension(:,:), allocatable :: kg_recv

  type :: prol                       
     integer :: fx,fy,tx,ty,from,to
     real*8 :: prob                  
  end type prol

  type(prol), dimension(:,:), allocatable :: plist

!!!!------------------------------------------------------!!!!
!!!!              KMC Process probabilities               !!!!
  !! GLOBAL
  real*8, dimension(:,:), allocatable :: vfe_f_g,vfe_a_g  !!!!
  real*8, dimension(:,:), allocatable :: vs_f_g,vs_a_g    !!!!
  real*8, dimension(:,:), allocatable :: fes_diss_g       !!!!
  real*8, dimension(:,:), allocatable :: v_diff_g         !!!!
  !! LOCAL
  real*8, dimension(:,:), allocatable :: vfe_f,vfe_a      !!!!
  real*8, dimension(:,:), allocatable :: vs_f,vs_a        !!!!
  real*8, dimension(:,:), allocatable :: fes_diss         !!!!
  real*8, dimension(:,:), allocatable :: v_diff           !!!!
!!!!------------------------------------------------------!!!!

end module kmc_data


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module diffusion_constants
  implicit none
  save

  real*8 :: D_S_met, D_S_mkw, D_S_pht, D_S_pyr, D_S_env
  real*8 :: D_Fe_met, D_Fe_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_env

end module diffusion_constants


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!


module gradients
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: delypyr,delzpyr

end module gradients


!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!=!
