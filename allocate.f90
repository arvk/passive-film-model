subroutine allocate_matrices()
  use commondata
  use fields
  use gradients
  use thermo_constants
  implicit none
  include 'mpif.h'
  !!####Allocate all matrices responsible for storing the field variables

  integer :: seed_size !! Size of the seed for the PRNG

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  call random_seed(size=seed_size)
  allocate(seed(seed_size)) ! Allocate PRNG Seed

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate global matrices
  allocate(met_g(psx_g,psy_g,psz_g))
  allocate(mkw_g(psx_g,psy_g,psz_g))
  allocate(pht_g(psx_g,psy_g,psz_g))
  allocate(pyr_g(psx_g,psy_g,psz_g))
  allocate(env_g(psx_g,psy_g,psz_g))
  allocate(mu_g(psx_g,psy_g,psz_g))
  allocate(opyr_g(psx_g,psy_g,psz_g))
  allocate(elpot_g(psx_g,psy_g,psz_g))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate gradient matrices
  allocate(delypyr(psx,psy,psz+(2*ghost_width)))
  allocate(delzpyr(psx,psy,psz+(2*ghost_width)))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate time-evolution matrices
  allocate(dmet_dt(psx,psy,psz+(2*ghost_width))) ; dmet_dt = 0.0d0
  allocate(dmkw_dt(psx,psy,psz+(2*ghost_width))) ; dmkw_dt = 0.0d0
  allocate(dpht_dt(psx,psy,psz+(2*ghost_width))) ; dpht_dt = 0.0d0
  allocate(dpyr_dt(psx,psy,psz+(2*ghost_width))) ; dpyr_dt = 0.0d0
  allocate(denv_dt(psx,psy,psz+(2*ghost_width))) ; denv_dt = 0.0d0
  allocate(dmu_dt(psx,psy,psz+(2*ghost_width)))  ; dmu_dt = 0.0d0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate field mobilities
  allocate(Mob_pf(0:(nphases-1),0:(nphases-1)))

  ! Allocate surface energies
  allocate(sigma(0:(nphases-1),0:(nphases-1)))

  ! Allocate realtive phase stabilities
  allocate(w_pf(0:(nphases-1)))

  ! Allocate realtive phase permittivities
  allocate(permittivity(0:(nphases-1)))

  ! Derivative of sulfur concentration with chemical potential
  allocate(drho_dmu(0:(nphases-1)))

  ! Sulfur density in different phases
  allocate(rhoS(0:(nphases-1)))

  ! Collected sulfidation rates for different FeS phases and environments
  allocate(sulf_rate_gas(0:(nphases-1)))
  allocate(sulf_rate_liq(0:(nphases-1)))
  allocate(sulf_rate(0:(nphases-1)))

end subroutine allocate_matrices
