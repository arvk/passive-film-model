subroutine read_geometry
  use commondata
  use fields
  use kmc_data
  use thermo_constants
  use laplacians

  implicit none

  integer :: x, y, z   ! Loop variables
  character*1 :: hash

  open (unit=101, file="PHT.out", status="old")
  read(101,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(101,*)  pht_g(x,y,z)
        end do
     end do
  end do
  close(101)  

  open (unit=102, file="ENV.out", status="old")
  read(102,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(102,*)  env_g(x,y,z)
        end do
     end do
  end do
  close(102)  

  open (unit=103, file="MET.out", status="old")
  read(103,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(103,*)  met_g(x,y,z)
        end do
     end do
  end do
  close(103)  

  open (unit=105, file="MUS.out", status="old")
  read(105,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           read(105,*)  mu_g(x,y,z)
        end do
     end do
  end do
  close(105)  

  open (unit=108, file="PYR.out", status="old")
  read(108,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           read(108,*)  pyr_g(x,y,z)
        end do
     end do
  end do
  close(108)  

  avg_mu_met = mus_met_pht_eqb + (R*T*0.1d0)
  avg_mu_env = (R*T*(0.0d0-0.995d0))
  avg_mu_pht = 0.5d0*(avg_mu_met+avg_mu_env)

  min_mu = R*T*(0-50.0d0)
  max_mu = avg_mu_env+(R*T*3.005)

end subroutine read_geometry


