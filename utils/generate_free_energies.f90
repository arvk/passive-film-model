program generate_free_energies
  implicit none

  real*8 :: T  !! Temperature (in Kelvin)
  real*8, parameter :: R = 8.314  !! Gas constant (in Joule/mol-K)

  real*8 :: mu !! Sulfur chemical potential
  real*8 :: g_met, g_mkw, g_pht, g_pyr, g_env

  integer :: i

  T = 500 ! (in K)

  call system("rm -rf free-energies.txt")
  open (unit = 132, file = "free-energies.txt", status = "new")

  do i = -20000,6000
     mu = R*T*(dble(i)/1000.0d0)  !! mu ranged from -30RT to 0

     g_met = 0.0d0 - (mu*0.0015)
     g_mkw = (1E-6*mu*mu) + (20.53*T) - 65060 - (mu*0.80d0)
     g_pht = (1E-6*mu*mu) + (20.53*T) - 72050 - (mu*1.0d0)
     g_pyr = (1E-9*mu*mu) + (50.355*T) - 98710 - (mu*2.0d0)
     g_env = 0.0d0 - (mu*0.0015)

     write(132,*) mu,g_met,g_mkw,g_pht,g_pyr,g_env

  end do

  close(132)

end program generate_free_energies

