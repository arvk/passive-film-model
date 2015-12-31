subroutine allocate_matrices()
  use commondata
  use fields
  use thermo_constants
  implicit none
#include <petsc/finclude/petscsys.h>
  !!####Allocate all matrices responsible for storing the field variables

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate global matrices
  allocate(met_g(psx_g,psy_g,psz_g))
  allocate(mkw_g(psx_g,psy_g,psz_g))
  allocate(pht_g(psx_g,psy_g,psz_g))
  allocate(pyr_g(psx_g,psy_g,psz_g))
  allocate(env_g(psx_g,psy_g,psz_g))
  allocate(mu_g(psx_g,psy_g,psz_g))
  allocate(elpot_g(psx_g,psy_g,psz_g))

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  ! Allocate field mobilities
  allocate(Mob_pf(0:(nphases-1),0:(nphases-1)))

  ! Allocate surface energies
  allocate(sigma(0:(nphases-1),0:(nphases-1)))

  ! Allocate realtive phase stabilities
  allocate(w_pf(0:(nphases-1)))

  ! Allocate realtive phase permittivities
  allocate(permittivity(0:(nphases-1)))

  ! Derivative of sulfur concentration with chemical potential
  allocate(drho_dmu(0:(nphases-1)))

  ! Sulfur density in different phases
  allocate(rhoS(0:(nphases-1)))

  ! Collected sulfidation rates for different FeS phases and environments
  allocate(sulf_rate_gas(0:(nphases-1)))
  allocate(sulf_rate_liq(0:(nphases-1)))
  allocate(sulf_rate(0:(nphases-1)))

end subroutine allocate_matrices
