module commondata
  implicit none
  save

  character*1 :: isrestart     ! Is the calculation a restarted one?
  integer :: ksx,ksy           ! Number of kMC gridpoints in 3 directions
  integer :: ksx_g,ksy_g       ! Number of kMC gridpoints in 3 directions
  integer :: psx_g,psy_g,psz_g ! Number of PF gridpoints in 3 directions
  integer :: psx,psy,psz       ! Number of PF gridpoints in 3 directions
  integer :: T                 ! Temperature of the simulation box
  real*8 :: pH_in              ! Scalar pH input
  integer :: rank, procs

  real*8, dimension(:,:,:), allocatable :: pht ! PF grid for pyrrhotite
  real*8, dimension(:,:,:), allocatable :: env ! PF grid for environment (liquid or gas)
  real*8, dimension(:,:,:), allocatable :: met ! PF grid for metal
  real*8, dimension(:,:,:), allocatable :: mu  ! mu_S grid
  real*8, dimension(:,:,:), allocatable :: ph  ! pH matrix in the simulation cell

  real*8, dimension(:,:,:), allocatable :: met_g 
  real*8, dimension(:,:,:), allocatable :: pht_g 
  real*8, dimension(:,:,:), allocatable :: env_g 
  real*8, dimension(:,:,:), allocatable :: mu_g 
  real*8, dimension(:,:,:), allocatable :: ph_g 

  integer, dimension(:,:), allocatable :: kg
  integer, dimension(:,:), allocatable :: kg_g
  integer, dimension(:,:), allocatable :: kg_recv


  integer :: nomc                         ! Number of kMC steps
  integer :: noimg = 50                   ! Number of output files
  integer :: freq_scale = 1750            ! KMC information is transferred every freq_scale steps
  integer :: kg_scale = 5                 ! Number of KMC grid points per PF grid
  real*8 :: dpf = 5E-9
  real*8 :: R = 8.3144621   

  integer :: swap_freq_pf = 5             ! Frequency with which MPI swaps are conducted for PF solving
  integer :: swap_freq_kmc = 10           ! Frequency with which MPI swaps are conducted for KMC solving

  real*8 :: avg_mu_pht
  real*8 :: avg_mu_env
  real*8 :: avg_mu_met

  real*8 :: max_mu, min_mu

  logical :: isparent

end module commondata




module laplacians
  implicit none
  save

  real*8, dimension(:,:,:), allocatable :: del2pht
  real*8, dimension(:,:,:), allocatable :: del2env
  real*8, dimension(:,:,:), allocatable :: del2met
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

  real*8, dimension(:,:,:), allocatable :: dpht_dt, denv_dt, dmet_dt, dmu_dt, dph_dt
  real*8, dimension(:,:,:), allocatable :: newpht, newenv, newmet, oldmu, newph

end module field_evolution
