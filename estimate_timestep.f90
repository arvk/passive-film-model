subroutine estimate_timestep
  use thermo_constants
  use commondata
  use diffusion_constants
  implicit none
#include <finclude/petscsys.h>
  !! **Estimate a stable timestep for phase field and chemical potential field evolution**

  PetscScalar :: dt_stable_diffusion    !! Maximum stable forward-euler timestep for integrating the chemical-potential field
  PetscScalar :: dt_stable_phase_field  !! Maximum stable forward-euler timestep for integrating the phase field
  PetscScalar :: max_M_sigma            !! Maximum of Mobility*sigma for all the phases
  PetscScalar :: fesphase, fesphase2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  dt_stable_diffusion = (dpf*dpf)/(2*max(D_S_pht, D_S_pyr, D_S_mkw, D_Fe_pht, D_Fe_pyr, D_Fe_mkw)) ! Standard Forward-Euler criterion

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  max_M_sigma = 0.0d0
  do fesphase = 0,(nphases-1)
     do fesphase2 = 0,(nphases-1)
        if (fesphase.ne.fesphase2) then
           max_M_sigma = max(max_M_sigma,Mob_pf(fesphase,fesphase2)*sigma(fesphase,fesphase2))
        end if
     end do
  end do

  dt_stable_phase_field = (dpf*dpf)/max_M_sigma ! The prefactor before the laplacian for PF evolution is the sum of M*sigma for each phase.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  dt = 0.1d0

  if (isroot) write(6,'(A,F9.5,A)') " INFO: Timestep for phase-field integration is ",dt, " seconds."

end subroutine estimate_timestep
