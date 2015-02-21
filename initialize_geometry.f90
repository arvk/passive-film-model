subroutine initialize_geometry()
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use kmc_data
  implicit none

  integer :: x,y,z
  integer :: pht_z_beg, pht_z_end
  integer :: rank_loop

  pht_z_beg = 9*(psz_g/16)
  pht_z_end = 11*(psz_g/16)

  avg_mu_met = mus_met_mkw_eqb - (R*T*0.1d0) 
  avg_mu_mkw = mus_mkw_pht_eqb - (R*T*0.1d0) 
  avg_mu_env = (R*T*(0.0d0-0.995d0))
  avg_mu_pht = 0.5d0*(avg_mu_mkw+avg_mu_env)

  min_mu = R*T*(0-50.0d0)
  max_mu = avg_mu_env+(R*T*3.005)


!! Initialize global phase-fraction fields
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,pht_z_beg
           met_g(x,y,z) = 1.0d0; mkw_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 0.0d0
        end do

        do z = pht_z_beg+1, pht_z_end
           met_g(x,y,z) = 0.0d0; mkw_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 1.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 0.0d0
        end do

        do z = pht_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0; mkw_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 1.0d0
        end do

     end do
  end do



!! Initialize global chemical-potential field
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,pht_z_beg
           mu_g(x,y,z) = avg_mu_met
        end do

        do z = pht_z_beg+1, pht_z_end
           mu_g(x,y,z) = avg_mu_met + (z-pht_z_beg)*150
        end do

        do z = pht_z_end+1,psz_g
           mu_g(x,y,z) = avg_mu_env
        end do

        do z = -2,2
           mu_g(x,y,pht_z_beg+z) = mu_g(x,y,pht_z_beg-3)+((mu_g(x,y,pht_z_beg+3)-mu_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
        end do
        
     end do
  end do


!! Initialize global pyrite orientation field

!! Populate the global orientation field with random numbers
  call random_seed(put=seed)
  call random_number(opyr_g)

  opyr_g = opyr_g * (3.14159265d0/2.0d0)

end subroutine initialize_geometry

