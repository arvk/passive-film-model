subroutine initialize_geometry()
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use kmc_data
  implicit none

  integer :: x,y,z
  integer :: met_z_end
  integer :: rank_loop
  real*8 :: avg_mu_met

  met_z_end = 15*(psz_g/16)

  avg_mu_met = mus_met_mkw_eqb - (R*T*0.5d0) 
  avg_mu_env = mus_mkw_pht_eqb - (R*T*2.5d0) 

  min_mu = R*T*(0-50.0d0)
  max_mu = avg_mu_env+(R*T*0.5)

!! Initialize global phase-fraction fields
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end
           met_g(x,y,z) = 1.0d0; mkw_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 0.0d0
        end do

        do z = met_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0; mkw_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 1.0d0
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

