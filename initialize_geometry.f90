subroutine initialize_geometry()
  use commondata
  use fields
  use thermo_constants
  use kmc_data
  implicit none

  integer :: x,y,z
  integer :: met_z_end
  integer :: rank_loop
  real*8 :: avg_mu_met
  real*8, parameter :: epsilon = 0.00001d0 ! A small number

  met_z_end = 15*(psz_g/16)

  avg_mu_met = mus_met_mkw_eqb - (R*T*0.5d0) 
  avg_mu_env = mus_mkw_pht_eqb - (R*T*2.5d0) 

!! Initialize global phase-fraction fields
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end
           met_g(x,y,z) = 1.0d0-epsilon; mkw_g(x,y,z) = 0.0d0+epsilon; pht_g(x,y,z) = 0.0d0+epsilon; pyr_g(x,y,z) = 0.0d0+epsilon; env_g(x,y,z) = 0.0d0+epsilon
        end do

        do z = met_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0+epsilon; mkw_g(x,y,z) = 0.0d0+epsilon; pht_g(x,y,z) = 0.0d0+epsilon; pyr_g(x,y,z) = 0.0d0+epsilon; env_g(x,y,z) = 1.0d0-epsilon
        end do

     end do
  end do



!! Initialize global chemical-potential field
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end
           mu_g(x,y,z) = avg_mu_met
        end do

        do z = met_z_end+1,psz_g
           mu_g(x,y,z) = avg_mu_env
        end do

     end do
  end do


!! Initialize global electrical potential
  elpot_g = metal_potential


!! Initialize global pyrite orientation field

!! Populate the global orientation field with random numbers
  call random_seed(put=seed)
  call random_number(opyr_g)

  opyr_g = opyr_g * (3.14159265d0/2.0d0)

end subroutine initialize_geometry

