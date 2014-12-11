subroutine allocate_matrices()
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use field_evolution
  use kmc_data
  implicit none
  include 'mpif.h'


!! ALLOCATE PF VARIABLES
!! =====================

!! Allocate global matrices  
  allocate(met_g(psx_g,psy_g,psz_g)); allocate(pht_g(psx_g,psy_g,psz_g)); allocate(pyr_g(psx_g,psy_g,psz_g)); allocate(env_g(psx_g,psy_g,psz_g))
  allocate(mu_g(psx_g,psy_g,psz_g))
  allocate(opyr_g(psx_g,psy_g,psz_g))


!! Allocate local matrices
  allocate(met(psx,psy,psz+2)); allocate(pht(psx,psy,psz+2)); allocate(pyr(psx,psy,psz+2)); allocate(env(psx,psy,psz+2)) 
  allocate(mu(psx,psy,psz+2)) 
  allocate(opyr(psx,psy,psz+2)) 

  
!! Allocate Laplacian matrices
  allocate(del2met(psx,psy,psz+2)); allocate(del2pht(psx,psy,psz+2)); allocate(del2pyr(psx,psy,psz+2)); allocate(del2env(psx,psy,psz+2))
  allocate(del2mu(psx,psy,psz+2))


!! Allocate time-evolution matrices
  allocate(dmet_dt(psx,psy,psz+2)); allocate(dpht_dt(psx,psy,psz+2)); allocate(dpyr_dt(psx,psy,psz+2)); allocate(denv_dt(psx,psy,psz+2))
  allocate(dmu_dt(psx,psy,psz+2))
  
  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dpyr_dt = 0.0d0 ; dmu_dt = 0.0d0


!! Allocate local matrices for updated phase-fields
  allocate(newmet(psx,psy,psz+2)); allocate(newpht(psx,psy,psz+2)); allocate(newpyr(psx,psy,psz+2)); allocate(newenv(psx,psy,psz+2))





!! ALLOCATE KMC VARIABLES
!! ======================

!! Allocate kMC grids
  allocate(kg_g(psx_g*kg_scale, psy_g*kg_scale))         ! Global kMC grid on Root
  allocate(kg((psx*kg_scale)+2, (psy*kg_scale/procs)+2)) ! Local kMC grid slices
  kg = 1
  allocate(kg_recv((psx_g*kg_scale)+2, psy_g*kg_scale))  ! Global matric to gather local kMC grids


!! Allocate process-list
  allocate(plist((ksx+2)*(ksy+2),12))


!! Allocate process probabilities     
  allocate(vfe_f_g(ksx_g, ksy_g))    ! V_Fe formation 
  allocate(vfe_a_g(ksx_g, ksy_g))    ! V_Fe annihilation  
  allocate(vs_f_g(ksx_g, ksy_g))     ! V_S formation   
  allocate(vs_a_g(ksx_g, ksy_g))     ! V_S annihilation
  allocate(fes_diss_g(ksx_g, ksy_g)) ! Fe+S dissolution
  allocate(v_diff_g(ksx_g, ksy_g))   ! V_Fe/V_S diffusion
                                       
  allocate(vfe_f(ksx+2, ksy+2))      ! V_Fe formation 
  allocate(vfe_a(ksx+2, ksy+2))      ! V_Fe annihilation  
  allocate(vs_f(ksx+2, ksy+2))       ! V_S formation   
  allocate(vs_a(ksx+2, ksy+2))       ! V_S annihilation
  allocate(fes_diss(ksx+2, ksy+2))   ! Fe+S dissolution
  allocate(v_diff(ksx+2, ksy+2))     ! V_Fe/V_S diffusion

  vfe_f_g = 0 ; vfe_a_g = 0 ; vs_f_g = 0 ; vs_a_g = 0 ; fes_diss_g = 0 ; v_diff_g = 0
  vfe_f = 0 ; vfe_a = 0 ; vs_f = 0 ; vs_a = 0 ; fes_diss = 0 ; v_diff = 0


end subroutine allocate_matrices
