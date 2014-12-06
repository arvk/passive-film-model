subroutine allocate_matrices()
  use commondata
  use laplacians
  use thermo_constants
  use field_evolution
  use kmc_data
  implicit none
  include 'mpif.h'

  allocate(mu(psx,psy,psz+2)) ; allocate(met(psx,psy,psz+2))
  allocate(pht(psx,psy,psz+2)) ; allocate(env(psx,psy,psz+2)) 
  allocate(ph(psx,psy,psz+2))

  allocate(met_g(psx_g,psy_g,psz_g))
  allocate(pht_g(psx_g,psy_g,psz_g))
  allocate(env_g(psx_g,psy_g,psz_g))
  allocate(mu_g(psx_g,psy_g,psz_g))
  allocate(ph_g(psx_g,psy_g,psz_g))

  allocate(del2pht(psx,psy,psz+2)) 
  allocate(del2env(psx,psy,psz+2))
  allocate(del2met(psx,psy,psz+2))
  allocate(del2mu(psx,psy,psz+2))
  allocate(del2ph(psx,psy,psz+2))

  allocate(dpht_dt(psx,psy,psz+2)); allocate(denv_dt(psx,psy,psz+2)) 
  allocate(dmet_dt(psx,psy,psz+2)); allocate(dmu_dt(psx,psy,psz+2))
  allocate(dph_dt(psx,psy,psz+2))

  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dmu_dt = 0.0d0 ; dph_dt = 0.0d0

!  allocate(plist(psx*kg_scale*psy*kg_scale,12))
  allocate(plist((ksx+2)*(ksy+2),12))
  allocate(kg_g(psx_g*kg_scale, psy_g*kg_scale))
  allocate(kg_recv((psx_g*kg_scale)+2, psy_g*kg_scale))
  allocate(kg((psx*kg_scale)+2, (psy*kg_scale/procs)+2))
  
  kg = 1

  allocate(newpht(psx,psy,psz+2)) ; allocate(newmet(psx,psy,psz+2))
  allocate(newenv(psx,psy,psz+2)) ; allocate(oldmu(psx,psy,psz+2))
  allocate(newph(psx,psy,psz+2)) 
 
  oldmu = 0.0d0

!!!!-----------------------------------!!!!
!!!!     KMC Process probabilities     !!!!
  allocate(vfe_f_g(ksx_g, ksy_g))      !!!!
  allocate(vfe_a_g(ksx_g, ksy_g))      !!!!
  allocate(vs_f_g(ksx_g, ksy_g))       !!!!
  allocate(vs_a_g(ksx_g, ksy_g))       !!!!
  allocate(fes_diss_g(ksx_g, ksy_g))   !!!!
  allocate(v_diff_g(ksx_g, ksy_g))     !!!!
                                       !!!!
  allocate(vfe_f(ksx+2, ksy+2))        !!!!
  allocate(vfe_a(ksx+2, ksy+2))        !!!!
  allocate(vs_f(ksx+2, ksy+2))         !!!!
  allocate(vs_a(ksx+2, ksy+2))         !!!!
  allocate(fes_diss(ksx+2, ksy+2))     !!!!
  allocate(v_diff(ksx+2, ksy+2))       !!!!
!!!!-----------------------------------!!!!

vfe_f_g = 0 ; vfe_a_g = 0 ; vs_f_g = 0 ; vs_a_g = 0 ; fes_diss_g = 0 ; v_diff_g = 0
vfe_f = 0 ; vfe_a = 0 ; vs_f = 0 ; vs_a = 0 ; fes_diss = 0 ; v_diff = 0


end subroutine allocate_matrices
