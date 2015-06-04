subroutine thermo
  use thermo_constants
  use commondata
  use fields
  implicit none

  integer :: mu_x_10k                              ! Integer equal to 10000*mu (Loop)
  real*8 :: my_mu                                  ! Chemical potential
  real*8 :: my_met, my_mkw, my_pht, my_pyr         ! Free energy of each phase
  real*8 :: min_met_mkw, min_mkw_pht, min_pht_pyr  ! Difference in free energy between two phases

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  !! Estimate the phase-stability region
  min_met_mkw = 50000.0d0; min_mkw_pht = 50000.0d0; min_pht_pyr = 50000.0d0

  do mu_x_10k = -300000,300000
     my_mu = R*T*(dble(mu_x_10k)/10000.0d0)   ! my_mu = -30RT to 30 RT

     my_met = 0.0d0 - (my_mu*0.0015)
     my_mkw = (1E-6*my_mu*my_mu) + (20.53*T) - 65060 - (my_mu*0.80d0)
     my_pht = (1E-6*my_mu*my_mu) + (20.53*T) - 72050 - (my_mu*1.0d0)
     my_pyr = (1E-9*my_mu*my_mu) + (50.355*T) - 98710 - (my_mu*2.0d0)

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

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  !! Density of sulfur in each phase (no. of moles of S / m^3 )
  rho_mkw = 48683.0d0           ! REF= Mkw density data from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
  rho_met = 0.0015d0*140401     ! REF= 0.15% assumed
  rho_pht = 52275.0d0           ! REF= de Villiers et. al. American Minerologist, 95(1): 148-152
  rho_pyr = 2.0d0*41667.0d0     ! REF= Muscat et. al. PRB 65(5):054107, 2002
  rho_env = 0.0015d0*(13303/T)  ! REF= 0.15% assumed


end subroutine thermo

