subroutine calculate_sulfidation_rates
  use thermo_constants
  use commondata
  use fields
  implicit none

  !! Collected sulfidation rates for different FeS phases and environments
  sulf_rate_gas(nmet) = 10**((0.00473*T)-5.645+(0.4*(avg_mu_env+63562)/(R*T))) !! Ref = Assessing Corrosion in Oil Refining and Petrochemical Processing, Materials Research, Vol 7, No 1, pp. 163-173, 2004
  sulf_rate_gas(nmet) = max(sulf_rate_gas(nmet)*1E-9,0.0d0)

  sulf_rate_gas(nmkw) = 0.01372 + 0.04356*(exp(avg_mu_env/(R*T))) !! Ref = Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
  sulf_rate_gas(nmkw) = max(sulf_rate_gas(nmkw)*1E-9,0.0d0)

  sulf_rate_gas(npht) = exp(-(11766/T)-0.6478)*1E9 !! Ref = Mechanisms Governing the Growth, Reactivity and Stability of Iron Sulfides, Ph.D Thesis William Herbert, MIT
  sulf_rate_gas(npht) = max(sulf_rate_gas(npht)*1E-9,0.0d0)

  sulf_rate_gas(npyr) = 7.45E8 * exp(-(98400/(R*T))) !! Ref = Kinetics of sulfidation of chalcopyrite with gaseous sulfur. Padilla R. et. al., Met Trans B, Vol 34B, Feb 2003, 61-68,
  sulf_rate_gas(npyr) = max(sulf_rate_gas(npyr)*1E-9,0.0d0)

  sulf_rate_gas(nenv) = 0.0d0 !! DUMMY VALUES

  sulf_rate_liq(nmet) = 0.0666*(1 + ((avg_mu_env+63562)/(R*T))) !! Ref = Kinetics of iron sulfiede and mixed iron sulfide/carbonate scale precipitation in CO2/H2S corrosion, Corrosion 2006, Paper 06644
  sulf_rate_liq(nmet) = max(sulf_rate_liq(nmet)*1E-9,0.0d0)

  sulf_rate_liq(nmkw) = 0.1332*(1 + (2*(avg_mu_env+63562)/(R*T))) !! Ref = Mechanistic model of H2S corrosion of mild steel
  sulf_rate_liq(nmkw) = max(sulf_rate_liq(nmkw)*1E-9,0.0d0)

  sulf_rate_liq(npht) = 2.41628 !! Ref = Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
  sulf_rate_liq(npht) = max(sulf_rate_liq(npht)*1E-9,0.0d0)

  sulf_rate_liq(npyr) = 0.003543 !! Ref = Crystal growth of pyrite in Aqueous solutions. Inhibition by organophosphorous compounds, Harmandas NG. et. al., Langmuir 14, 1250-1255, 1998.
  sulf_rate_liq(npyr) = max(sulf_rate_liq(npyr)*1E-9,0.0d0)

  sulf_rate_liq(nenv) = 0.0d0 !! DUMMY VALUES

end subroutine calculate_sulfidation_rates

