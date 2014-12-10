subroutine estimate_timestep
  use fields
  use diffusion_constants

  real*8 :: dt_stable       ! Maximum stable forward-euler timestep

  dt_stable = (dpf*dpf)/(2*dmax1(D_S_pht, D_S_pyr, D_Fe_pht, D_Fe_pyr))

  dt = dt_stable/5000.0d0
  write(6,*) dt

end subroutine estimate_timestep
