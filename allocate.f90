subroutine allocate_matrices()
  use commondata
  use fields
  use gradients
  use thermo_constants
  use kmc_data
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

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Allocate kMC grids
  allocate(kg_g(psx_g*kg_scale, psy_g*kg_scale))                  ! Global kMC grid on Root
  allocate(kg((psx*kg_scale)+2, (psy*kg_scale/procs)+2)) ; kg = 1 ! Local kMC grid slices
  allocate(kg_recv((psx_g*kg_scale)+2, psy_g*kg_scale))           ! Global matric to gather local kMC grids

  allocate(plist((ksx+2)*(ksy+2),12)) ! Allocate process-list

!! Allocate process probabilities
  allocate(vfe_f_g(ksx_g, ksy_g))    ; vfe_f_g = 0    ! V_Fe formation
  allocate(vfe_a_g(ksx_g, ksy_g))    ; vfe_a_g = 0    ! V_Fe annihilation
  allocate(vs_f_g(ksx_g, ksy_g))     ; vs_f_g = 0     ! V_S formation
  allocate(vs_a_g(ksx_g, ksy_g))     ; vs_a_g = 0     ! V_S annihilation
  allocate(fes_diss_g(ksx_g, ksy_g)) ; fes_diss_g = 0 ! Fe+S dissolution
  allocate(v_diff_g(ksx_g, ksy_g))   ; v_diff_g = 0   ! V_Fe/V_S diffusion

  allocate(vfe_f(ksx+2, ksy+2))      ; vfe_f = 0      ! V_Fe formation
  allocate(vfe_a(ksx+2, ksy+2))      ; vfe_a = 0      ! V_Fe annihilation
  allocate(vs_f(ksx+2, ksy+2))       ; vs_f = 0       ! V_S formation
  allocate(vs_a(ksx+2, ksy+2))       ; vs_a = 0       ! V_S annihilation
  allocate(fes_diss(ksx+2, ksy+2))   ; fes_diss = 0   ! Fe+S dissolution
  allocate(v_diff(ksx+2, ksy+2))     ; v_diff = 0     ! V_Fe/V_S diffusion

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

end subroutine allocate_matrices
