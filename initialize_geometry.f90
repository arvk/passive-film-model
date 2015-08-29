subroutine initialize_geometry()
  use commondata
  use fields
  use thermo_constants
  use kmc_data
  implicit none

  integer :: x,y,z                               ! Index for x-, y- and z-direction (Loop)
  integer :: met_z_end                           ! Location of boundary between metal and environment
  real*8 :: avg_mu_met                           ! Chemical potential in the metal region
  real*8, parameter :: infinitesimal = 0.00001d0 ! A hard-coded 'small' number

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  met_z_end = 15*(psz_g/16)

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end         ! Initialize metal region
           met_g(x,y,z) = 1.0d0-infinitesimal
           mkw_g(x,y,z) = 0.0d0+infinitesimal
           pht_g(x,y,z) = 0.0d0+infinitesimal
           pyr_g(x,y,z) = 0.0d0+infinitesimal
           env_g(x,y,z) = 0.0d0+infinitesimal
        end do

        do z = met_z_end+1,psz_g   ! Initialize environmental region
           met_g(x,y,z) = 0.0d0+infinitesimal
           mkw_g(x,y,z) = 0.0d0+infinitesimal
           pht_g(x,y,z) = 0.0d0+infinitesimal
           pyr_g(x,y,z) = 0.0d0+infinitesimal
           env_g(x,y,z) = 1.0d0-infinitesimal
        end do

     end do
  end do

  metal_amount = sum(met_g)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  avg_mu_met = mus_met_mkw_eqb - (R*T*0.5d0)
  avg_mu_env = mus_mkw_pht_eqb - (R*T*2.5d0)

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end
           mu_g(x,y,z) = avg_mu_met    ! Initialize metal region
        end do

        do z = met_z_end+1,psz_g
           mu_g(x,y,z) = avg_mu_env    ! Initialize environmental region
        end do

     end do
  end do

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  elpot_g = metal_potential            ! Initialize global electrical potential

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  pH_g = pH_in                         ! Initialize global pH field

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call random_seed(put=seed)
  call random_number(opyr_g)
  opyr_g = opyr_g * (3.14159265d0/2.0d0)   ! Populate the global orientation field with random numbers

end subroutine initialize_geometry

