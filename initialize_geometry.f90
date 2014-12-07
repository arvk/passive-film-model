subroutine initialize_geometry()
  use commondata
  use laplacians
  use thermo_constants
  use kmc_data

  implicit none

  integer :: x,y,z
  integer :: pht_z_beg, pht_z_end
  integer :: rank_loop
  real*8 :: dt = 5E-2

  pht_z_beg = 543*(psz_g/576)
  pht_z_end = 556*(psz_g/576)

  avg_mu_met = mus_met_pht_eqb + (R*T*0.1d0) 
!  avg_mu_env = (R*T*(0.0d0-0.99d0))
  avg_mu_env = (R*T*(0.0d0-0.995d0))
  avg_mu_pht = 0.5d0*(avg_mu_met+avg_mu_env)

  min_mu = R*T*(0-50.0d0)
  max_mu = avg_mu_env+(R*T*3.005)

  write(6,*) avg_mu_met, avg_mu_pht, avg_mu_env, mus_met_pht_eqb, mus_pht_pyr_eqb

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,pht_z_beg
           met_g(x,y,z) = 1.0d0      ! Initialize metal
           pht_g(x,y,z) = 0.0d0
           env_g(x,y,z) = 0.0d0
           pyr_g(x,y,z) = 0.0d0
           mu_g(x,y,z) = avg_mu_met
           ph_g(x,y,z) = pH_in
        end do

        do z = pht_z_beg+1, pht_z_end
           met_g(x,y,z) = 0.0d0
           pht_g(x,y,z) = 1.0d0      ! Initialize pyrrhotite
           env_g(x,y,z) = 0.0d0
           pyr_g(x,y,z) = 0.0d0
           mu_g(x,y,z) = avg_mu_met + (z-pht_z_beg)*10
           ph_g(x,y,z) = pH_in
        end do

        do z = pht_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0
           pht_g(x,y,z) = 0.0d0
           env_g(x,y,z) = 1.0d0      ! Initialize environment
           pyr_g(x,y,z) = 0.0d0
           mu_g(x,y,z) = avg_mu_env
           ph_g(x,y,z) = pH_in
        end do

        do z = -2,2
           met_g(x,y,pht_z_beg+z) = met_g(x,y,pht_z_beg-3)+((met_g(x,y,pht_z_beg+3)-met_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           pht_g(x,y,pht_z_beg+z) = pht_g(x,y,pht_z_beg-3)+((pht_g(x,y,pht_z_beg+3)-pht_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           env_g(x,y,pht_z_beg+z) = env_g(x,y,pht_z_beg-3)+((env_g(x,y,pht_z_beg+3)-env_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           mu_g(x,y,pht_z_beg+z) = mu_g(x,y,pht_z_beg-3)+((mu_g(x,y,pht_z_beg+3)-mu_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           ph_g(x,y,pht_z_beg+z) = ph_g(x,y,pht_z_beg-3)+((ph_g(x,y,pht_z_beg+3)-ph_g(x,y,pht_z_beg-3))*(0.2*(z+2)))

           ! met_g(x,y,pht_z_end+z) = met_g(x,y,pht_z_end-3)+((met_g(x,y,pht_z_end+3)-met_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           ! pht_g(x,y,pht_z_end+z) = pht_g(x,y,pht_z_end-3)+((pht_g(x,y,pht_z_end+3)-pht_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           ! env_g(x,y,pht_z_end+z) = env_g(x,y,pht_z_end-3)+((env_g(x,y,pht_z_end+3)-env_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           ! mu_g(x,y,pht_z_end+z) = mu_g(x,y,pht_z_end-3)+((mu_g(x,y,pht_z_end+3)-mu_g(x,y,pht_z_end-3))*(0.2*(z+2)))
           ! ph_g(x,y,pht_z_end+z) = ph_g(x,y,pht_z_end-3)+((ph_g(x,y,pht_z_end+3)-ph_g(x,y,pht_z_end-3))*(0.2*(z+2)))
        end do
        
     end do
  end do



end subroutine initialize_geometry


