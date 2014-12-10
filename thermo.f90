subroutine thermo
  use thermo_constants
  use commondata
  use fields
  implicit none

  integer :: i              ! Loop index
  real*8 :: my_mu           ! Chemical potential
  real*8 :: my_pht, my_met, my_pyr  ! Free energy of each phase
  real*8 :: min_met_pht, min_pht_pyr

  min_met_pht = 50000.0d0; min_pht_pyr = 50000.0d0

  do i = -30000,30000
     my_mu = R*T*(dble(i)/1000.0d0)

     my_pht = (1E-6*my_mu*my_mu) + (20.53*T) - 72050 - (my_mu*1.0d0)
     my_met = 0.0d0 - (my_mu*0.0015)
     my_pyr = (1E-9*my_mu*my_mu) + (50.355*T) - 98710 - (my_mu*2)

     if (abs(my_met-my_pht) .lt. min_met_pht) then
        mus_met_pht_eqb = my_mu
        min_met_pht = abs(my_met-my_pht)
     end if

     if (abs(my_pht-my_pyr) .lt. min_pht_pyr) then
        mus_pht_pyr_eqb = my_mu
        min_pht_pyr = abs(my_pht-my_pyr)
     end if
     
  end do

end subroutine thermo

