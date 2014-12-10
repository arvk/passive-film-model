subroutine estimate_timestep
  use fields
  use commondata
  use diffusion_constants
  implicit none

  real*8 :: dt_stable_diffusion    ! Maximum stable forward-euler timestep for integrating the mu field
  real*8 :: dt_stable_phase_field  ! Maximum stable forward-euler timestep for integrating the phase field

  dt_stable_diffusion = (dpf*dpf)/(2*max(D_S_pht, D_S_pyr, D_Fe_pht, D_Fe_pyr)) ! Standard Forward-Euler criterion
  dt_stable_phase_field = -(dpf*dpf)/(8*min(M_pht_env*sigma_pht_env,M_pht_met*sigma_pht_met,M_pht_pyr*sigma_pht_pyr,M_met_env*sigma_met_env,M_met_pyr*sigma_met_pyr,M_env_pyr*sigma_env_pyr))
  !! The prefactor before the laplacian for PF evolution is the sum of M*sigma for each phase. So we conservatively use 4*max(M*sigma) to estimate the stable time. Also, the sigmas are negative, so we use "-min" instead of "max"

  dt = min(dt_stable_phase_field,dt_stable_diffusion)/2.0d0

  if (isroot) then
     write(6,'(A,E19.12,A)') "Timestep for phase-field integration is ",dt, " seconds."
  end if

end subroutine estimate_timestep
