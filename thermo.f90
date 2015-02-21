subroutine thermo
  use thermo_constants
  use commondata
  use fields
  implicit none

  integer :: i              ! Loop index
  real*8 :: my_mu           ! Chemical potential
  real*8 :: my_met, my_mkw, my_pht, my_pyr  ! Free energy of each phase
  real*8 :: min_met_mkw, min_mkw_pht, min_pht_pyr


  !! Estimate the phase-stability region

  min_met_mkw = 50000.0d0; min_mkw_pht = 50000.0d0; min_pht_pyr = 50000.0d0

  do i = -30000,30000
     my_mu = R*T*(dble(i)/1000.0d0)

     my_met = 0.0d0 - (my_mu*0.0015)
     my_mkw = (1E-6*my_mu*my_mu) + (20.53*T) - 70060 - (my_mu*0.95d0)
     my_pht = (1E-6*my_mu*my_mu) + (20.53*T) - 72050 - (my_mu*1.0d0)
     my_pyr = (1E-9*my_mu*my_mu) + (50.355*T) - 98710 - (my_mu*2)

     if (abs(my_met-my_mkw) .lt. min_met_mkw) then
        mus_met_mkw_eqb = my_mu
        min_met_mkw = abs(my_met-my_mkw)
     end if

     if (abs(my_mkw-my_pht) .lt. min_mkw_pht) then
        mus_mkw_pht_eqb = my_mu
        min_mkw_pht = abs(my_mkw-my_pht)
     end if

     if (abs(my_pht-my_pyr) .lt. min_pht_pyr) then
        mus_pht_pyr_eqb = my_mu
        min_pht_pyr = abs(my_pht-my_pyr)
     end if
     
  end do


  !! Density of sulfur in each FeS_x phase (number of moles of sulfur / m^3 of phase)

  rho_mkw = 0.95d0*48683.0d0 !! Mkw density data from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683 
  rho_met = 0.0015d0*140401
  rho_pht = 52275.0d0
  rho_pyr = 2.0d0*41667.0d0
  rho_env = 0.0015d0*(13303/T)
           

end subroutine thermo

