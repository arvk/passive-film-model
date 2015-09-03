subroutine allocate_matrices()
  use commondata
  use fields
  use gradients
  use thermo_constants
  implicit none
  include 'mpif.h'

  integer :: seed_size ! For seeding the PRNG

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call random_seed(size=seed_size)
  allocate(seed(seed_size)) ! Allocate PRNG Seed

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate global matrices
  allocate(met_g(psx_g,psy_g,psz_g))
  allocate(mkw_g(psx_g,psy_g,psz_g))
  allocate(pht_g(psx_g,psy_g,psz_g))
  allocate(pyr_g(psx_g,psy_g,psz_g))
  allocate(env_g(psx_g,psy_g,psz_g))
  allocate(mu_g(psx_g,psy_g,psz_g))
  allocate(opyr_g(psx_g,psy_g,psz_g))
  allocate(elpot_g(psx_g,psy_g,psz_g))

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate local matrices
  allocate(met(psx,psy,psz+(2*ghost_width)))
  allocate(mkw(psx,psy,psz+(2*ghost_width)))
  allocate(pht(psx,psy,psz+(2*ghost_width)))
  allocate(pyr(psx,psy,psz+(2*ghost_width)))
  allocate(env(psx,psy,psz+(2*ghost_width)))
  allocate(mu(psx,psy,psz+(2*ghost_width)))
  allocate(opyr(psx,psy,psz+(2*ghost_width)))
  allocate(voids(psx,psy,psz+(2*ghost_width))) ; voids = 0.0d0
  allocate(elpot(psx,psy,psz+(2*ghost_width)))

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate gradient matrices
  allocate(delypyr(psx,psy,psz+(2*ghost_width)))
  allocate(delzpyr(psx,psy,psz+(2*ghost_width)))

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate time-evolution matrices
  allocate(dmet_dt(psx,psy,psz+(2*ghost_width))) ; dmet_dt = 0.0d0
  allocate(dmkw_dt(psx,psy,psz+(2*ghost_width))) ; dmkw_dt = 0.0d0
  allocate(dpht_dt(psx,psy,psz+(2*ghost_width))) ; dpht_dt = 0.0d0
  allocate(dpyr_dt(psx,psy,psz+(2*ghost_width))) ; dpyr_dt = 0.0d0
  allocate(denv_dt(psx,psy,psz+(2*ghost_width))) ; denv_dt = 0.0d0
  allocate(dmu_dt(psx,psy,psz+(2*ghost_width)))  ; dmu_dt = 0.0d0

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate local matrices for updated phase-fields
  allocate(newmet(psx,psy,psz+(2*ghost_width)))
  allocate(newmkw(psx,psy,psz+(2*ghost_width)))
  allocate(newpht(psx,psy,psz+(2*ghost_width)))
  allocate(newpyr(psx,psy,psz+(2*ghost_width)))
  allocate(newenv(psx,psy,psz+(2*ghost_width)))

end subroutine allocate_matrices
