subroutine calculate_sulfidation_rates
  use thermo_constants
  use commondata
  use fields
  implicit none
  !! **Rates of surface sulfidation for different \(FeS\) phases and environments**
  !--------------------------------------------------------------------------------

  !!#### Gaseous

  !!$$k_{metal}^{gas} = 10^{(0.00473 T) - 5.645 + \frac{0.4 \times (\mu_{env}+63562)}{R T}} \ nm/s $$
  !!Reference: Assessing Corrosion in Oil Refining and Petrochemical Processing, Materials Research, Vol 7, No 1, pp. 163-173, 2004
  !!---
  sulf_rate_gas(nmet) = 10**((0.00473*T)-5.645+(0.4*(avg_mu_env+63562)/(R*T)))
  sulf_rate_gas(nmet) = max(sulf_rate_gas(nmet)*1E-9,0.0d0)

  !! $$k_{mackinawite}^{gas} = 0.01372 + 0.04356 \times e^{\frac{\mu_{env}}{R T}} \ nm/s $$
  !! Reference: Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
  !!---
  sulf_rate_gas(nmkw) = 0.01372 + 0.04356*(exp(avg_mu_env/(R*T)))
  sulf_rate_gas(nmkw) = max(sulf_rate_gas(nmkw)*1E-9,0.0d0)

  !! $$k_{pyrrhotite}^{gas} = 10^9 \times e^{\frac{-11766}{T}-0.6478} \ nm/s $$
  !! Reference: Mechanisms Governing the Growth, Reactivity and Stability of Iron Sulfides, Ph.D Thesis William Herbert, MIT
  !!---
  sulf_rate_gas(npht) = exp(-(11766/T)-0.6478)*1E9
  sulf_rate_gas(npht) = max(sulf_rate_gas(npht)*1E-9,0.0d0)

  !! $$k_{pyrite}^{gas} = 7.45 \times 10^8 \times e^{\frac{-98400}{R T}} \ nm/s $$
  !! Reference: Kinetics of sulfidation of chalcopyrite with gaseous sulfur. Padilla R. et. al., Met Trans B, Vol 34B, Feb 2003, 61-68
  !!---
  sulf_rate_gas(npyr) = 7.45E8 * exp(-(98400/(R*T)))
  sulf_rate_gas(npyr) = max(sulf_rate_gas(npyr)*1E-9,0.0d0)

  sulf_rate_gas(nenv) = 0.0d0 ! DUMMY VALUES


  !!#### Aqueous

  !!$$k_{metal}^{liq} = 0.0666 \times \left( 1 + \frac{\mu_{env} + 63562}{R T} \right) $$
  !!Reference: Kinetics of iron sulfide and mixed iron sulfide/carbonate scale precipitation in CO2/H2S corrosion, Corrosion 2006, Paper 06644
  !!---
  sulf_rate_liq(nmet) = 0.0666*(1 + ((avg_mu_env+63562)/(R*T)))
  sulf_rate_liq(nmet) = max(sulf_rate_liq(nmet)*1E-9,0.0d0)

  !!$$k_{mackinawite}^{liq} = 0.1332 \times \left(1 + \frac{2 \times (\mu_{env} + 63562)}{R T} \right) $$
  !!Reference: Mechanistic model of H2S corrosion of mild steel
  !!---
  sulf_rate_liq(nmkw) = 0.1332*(1 + (2*(avg_mu_env+63562)/(R*T)))
  sulf_rate_liq(nmkw) = max(sulf_rate_liq(nmkw)*1E-9,0.0d0)

  !!$$k_{pyrrhotite}^{liq} = 2.41628 \ nm/s $$
  !!Reference: Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
  !!---
  sulf_rate_liq(npht) = 2.41628
  sulf_rate_liq(npht) = max(sulf_rate_liq(npht)*1E-9,0.0d0)

  !!$$k_{pyrite}^{liq} = 0.003543 \ nm/s $$
  !!Reference: Crystal growth of pyrite in Aqueous solutions. Inhibition by organophosphorous compounds, Harmandas NG. et. al., Langmuir 14, 1250-1255, 1998
  !!---
  sulf_rate_liq(npyr) = 0.003543
  sulf_rate_liq(npyr) = max(sulf_rate_liq(npyr)*1E-9,0.0d0)

  sulf_rate_liq(nenv) = 0.0d0 ! DUMMY VALUES

end subroutine calculate_sulfidation_rates

