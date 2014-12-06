subroutine initialize_geometry()
  use commondata
  use laplacians
  use thermo_constants
  use kmc_data

  implicit none

  integer :: x,y,z
  integer :: pht_z_beg, pht_z_end
  integer :: rank_loop


  pht_z_beg = 8*(psz_g/20)
  pht_z_end = 16*(psz_g/20)

  max_mu = R*T*5.0d0
  min_mu = R*T*(0-50.0d0)

  avg_mu_met = mus_met_pht_eqb - (R*T*0.50) 
  avg_mu_env = mus_met_pht_eqb + (R*T*0.75) 
  avg_mu_pht = mus_met_pht_eqb + (R*T*0.50) 

  write(6,*) avg_mu_met, avg_mu_pht, avg_mu_env

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,pht_z_beg
           met_g(x,y,z) = 1.0d0      ! Initialize metal
           pht_g(x,y,z) = 0.0d0
           env_g(x,y,z) = 0.0d0
           mu_g(x,y,z) = avg_mu_met
           ph_g(x,y,z) = pH_in
        end do

        do z = pht_z_beg+1, pht_z_end
           met_g(x,y,z) = 0.0d0
           pht_g(x,y,z) = 1.0d0      ! Initialize pyrrhotite
           env_g(x,y,z) = 0.0d0
           mu_g(x,y,z) = avg_mu_met + (z-pht_z_beg)*((avg_mu_env-avg_mu_met)/(pht_z_end-pht_z_beg))
           ph_g(x,y,z) = pH_in
        end do

        do z = pht_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0
           pht_g(x,y,z) = 0.0d0
           env_g(x,y,z) = 1.0d0      ! Initialize environment
           mu_g(x,y,z) = avg_mu_env
           ph_g(x,y,z) = pH_in
        end do

        do z = -2,2
           met_g(x,y,pht_z_beg+z) = met_g(x,y,pht_z_beg-3)+((met_g(x,y,pht_z_beg+3)-met_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           pht_g(x,y,pht_z_beg+z) = pht_g(x,y,pht_z_beg-3)+((pht_g(x,y,pht_z_beg+3)-pht_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           env_g(x,y,pht_z_beg+z) = env_g(x,y,pht_z_beg-3)+((env_g(x,y,pht_z_beg+3)-env_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           mu_g(x,y,pht_z_beg+z) = mu_g(x,y,pht_z_beg-3)+((mu_g(x,y,pht_z_beg+3)-mu_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           ph_g(x,y,pht_z_beg+z) = ph_g(x,y,pht_z_beg-3)+((ph_g(x,y,pht_z_beg+3)-ph_g(x,y,pht_z_beg-3))*(0.2*(z+2)))

           met_g(x,y,pht_z_end+z) = met_g(x,y,pht_z_end-3)+((met_g(x,y,pht_z_end+3)-met_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           pht_g(x,y,pht_z_end+z) = pht_g(x,y,pht_z_end-3)+((pht_g(x,y,pht_z_end+3)-pht_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           env_g(x,y,pht_z_end+z) = env_g(x,y,pht_z_end-3)+((env_g(x,y,pht_z_end+3)-env_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           mu_g(x,y,pht_z_end+z) = mu_g(x,y,pht_z_end-3)+((mu_g(x,y,pht_z_end+3)-mu_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           ph_g(x,y,pht_z_end+z) = ph_g(x,y,pht_z_end-3)+((ph_g(x,y,pht_z_end+3)-ph_g(x,y,pht_z_end-3))*(0.2*(z+2)))
        end do
        
     end do
  end do



end subroutine initialize_geometry


