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

  !! Random number generation for initializing the orientation field
  integer, dimension(:), allocatable :: seed
  integer :: n,istat
  integer, dimension(8) :: datetime


  pht_z_beg = 2*(psz_g/4)
  pht_z_end = 3*(psz_g/4)

  avg_mu_met = mus_met_pht_eqb + (R*T*0.1d0) 
  avg_mu_env = (R*T*(0.0d0-0.995d0))
  avg_mu_pht = 0.5d0*(avg_mu_met+avg_mu_env)

  min_mu = R*T*(0-50.0d0)
  max_mu = avg_mu_env+(R*T*3.005)


!! Initialize global phase-fraction fields
  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,pht_z_beg
           met_g(x,y,z) = 1.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 0.0d0
        end do

        do z = pht_z_beg+1, pht_z_end
           met_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 1.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 0.0d0
        end do

        do z = pht_z_end+1,psz_g
           met_g(x,y,z) = 0.0d0; pht_g(x,y,z) = 0.0d0; pyr_g(x,y,z) = 0.0d0; env_g(x,y,z) = 1.0d0
        end do

        do z = -2,2
           met_g(x,y,pht_z_beg+z) = met_g(x,y,pht_z_beg-3)+((met_g(x,y,pht_z_beg+3)-met_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           pht_g(x,y,pht_z_beg+z) = pht_g(x,y,pht_z_beg-3)+((pht_g(x,y,pht_z_beg+3)-pht_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
           env_g(x,y,pht_z_beg+z) = env_g(x,y,pht_z_beg-3)+((env_g(x,y,pht_z_beg+3)-env_g(x,y,pht_z_beg-3))*(0.2*(z+2)))
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

!! Initialize the random seed from /dev/random
  call random_seed(size=n)
  allocate(seed(n))

  call random_seed(get=seed)
  
  open(89, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(89) seed
     close(89)
  else
     call date_and_time(values=datetime)
     seed(n) = datetime(8); seed(1) = datetime(8)*datetime(7)*datetime(6)
  end if

  call random_seed(put=seed)

!! Populate the global orientation field with random numbers
  call random_number(opyr_g)

  opyr_g = opyr_g * (3.14d0/4.0d0)

end subroutine initialize_geometry

