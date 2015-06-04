subroutine estimate_timestep
  use thermo_constants
  use commondata
  use diffusion_constants
  implicit none

  real*8 :: dt_stable_diffusion    ! Maximum stable forward-euler timestep for integrating the mu field
  real*8 :: dt_stable_phase_field  ! Maximum stable forward-euler timestep for integrating the phase field
  real*8 :: max_M_sigma            ! Maximum of Mobility*sigma for all the phases

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  dt_stable_diffusion = (dpf*dpf)/(2*max(D_S_pht, D_S_pyr, D_S_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_mkw)) ! Standard Forward-Euler criterion

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  max_M_sigma = max(M_met_mkw*sigma_met_mkw,M_met_pht*sigma_met_pht,M_met_pyr*sigma_met_pyr_0,M_met_env*sigma_met_env)
  max_M_sigma = max(max_M_sigma,M_mkw_pht*sigma_mkw_pht,M_mkw_pyr*sigma_mkw_pyr_0,M_mkw_env*sigma_mkw_env)
  max_M_sigma = max(max_M_sigma,M_pht_pyr*sigma_pht_pyr_0,M_pht_env*sigma_pht_env)
  max_M_sigma = max(max_M_sigma,M_pyr_env*sigma_pyr_env_0)

  dt_stable_phase_field = (dpf*dpf)/max_M_sigma ! The prefactor before the laplacian for PF evolution is the sum of M*sigma for each phase.

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  dt = min(dt_stable_phase_field,dt_stable_diffusion)/0.08d0

  if (isroot) write(6,'(A,E19.12,A)') "Timestep for phase-field integration is ",dt, " seconds."

end subroutine estimate_timestep
