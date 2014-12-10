module fields
  implicit none
  save

  !! Define the extent of phase fields
  integer :: psx_g,psy_g,psz_g ! Number of PF gridpoints in 3 directions
  integer :: psx,psy,psz       ! Number of PF gridpoints in 3 directions
  
  !! Define local phase fields
  real*8, dimension(:,:,:), allocatable :: pht ! PF grid for pyrrhotite
  real*8, dimension(:,:,:), allocatable :: env ! PF grid for environment (liquid or gas)
  real*8, dimension(:,:,:), allocatable :: met ! PF grid for metal
  real*8, dimension(:,:,:), allocatable :: pyr ! PF grid for pyrite
  real*8, dimension(:,:,:), allocatable :: mu  ! mu_S grid
  real*8, dimension(:,:,:), allocatable :: ph  ! pH matrix in the simulation cell

  !! Define global phase fields
  real*8, dimension(:,:,:), allocatable :: met_g 
  real*8, dimension(:,:,:), allocatable :: pht_g 
  real*8, dimension(:,:,:), allocatable :: env_g 
  real*8, dimension(:,:,:), allocatable :: pyr_g 
  real*8, dimension(:,:,:), allocatable :: mu_g 
  real*8, dimension(:,:,:), allocatable :: ph_g 

end module fields





module commondata
  implicit none
  save

  character*1 :: isrestart     ! Is the calculation a restarted one?
  integer :: T                 ! Temperature of the simulation box
  real*8 :: pH_in              ! Scalar pH input
  integer :: rank, procs

  integer :: nomc                         ! Number of kMC steps
  integer :: noimg = 20                   ! Number of output files
  integer :: freq_scale = 1750000000      ! KMC information is transferred every freq_scale steps
  integer :: kg_scale = 5                 ! Number of KMC grid points per PF grid
  real*8 :: dpf = 5E-9
  real*8 :: R = 8.3144621   

  integer :: swap_freq_pf = 5             ! Frequency with which MPI swaps are conducted for PF solving
  integer :: swap_freq_kmc = 10           ! Frequency with which MPI swaps are conducted for KMC solving

  real*8 :: avg_mu_pht
  real*8 :: avg_mu_env
  real*8 :: avg_mu_met

  real*8 :: max_mu, min_mu

  logical :: isroot
  real*8 :: dt

end module commondata




module laplacians
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: del2pht
  real*8, dimension(:,:,:), allocatable :: del2env
  real*8, dimension(:,:,:), allocatable :: del2met
  real*8, dimension(:,:,:), allocatable :: del2pyr
  real*8, dimension(:,:,:), allocatable :: del2mu
  real*8, dimension(:,:,:), allocatable :: del2ph

end module laplacians



module thermo_constants
  implicit none
  save
  
  real*8 :: mus_met_pht_eqb, mus_pht_pyr_eqb

end module thermo_constants




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



module field_evolution
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: dpht_dt, denv_dt, dmet_dt, dpyr_dt, dmu_dt, dph_dt
  real*8, dimension(:,:,:), allocatable :: newpht, newenv, newmet, newpyr, oldmu, newph

end module field_evolution


module diffusion_constants
  implicit none
  save

  real*8 :: D_S_met, D_S_env, D_S_pht, D_S_pyr
  real*8 :: D_Fe_met, D_Fe_env, D_Fe_pht, D_Fe_pyr

end module diffusion_constants

