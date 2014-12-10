subroutine estimate_timestep
  use fields
  use diffusion_constants

  real*8 :: dt_stable_diffusion    ! Maximum stable forward-euler timestep for integrating the mu field
  real*8 :: dt_stable_phase_field  ! Maximum stable forward-euler timestep for integrating the phase field

  dt_stable_diffusion = (dpf*dpf)/(2*dmax1(D_S_pht, D_S_pyr, D_Fe_pht, D_Fe_pyr))
  dt_stable_pf = -(dpf*dpf/8)/(dmin1(M_pht_env*sigma_pht_env,M_pht_met*sigma_pht_met,M_pht_pyr*sigma_pht_pyr,M_met_env*sigma_met_env,M_met_pyr*sigma_met_pyr,M_env_pyr*sigma_env_pyr))

  dt = min(dt_stable_pf,dt_stable_diffusion)/2.0d0
  write(6,*) dt

end subroutine estimate_timestep
