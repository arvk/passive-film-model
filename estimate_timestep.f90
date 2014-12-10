subroutine estimate_timestep
  use commondata
  use diffusion_constants

  real*8 :: dt_stable       ! Maximum stable forward-euler timestep

  dt_stable = (dpf*dpf)/(2*dmax1(D_S_met, D_S_pht, D_S_pyr,D_Fe_met, D_Fe_pht, D_Fe_pyr))

  dt = dt_stable/2.0d0

end subroutine estimate_timestep