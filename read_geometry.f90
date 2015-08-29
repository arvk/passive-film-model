subroutine read_geometry
  use commondata
  use fields
  use kmc_data
  use thermo_constants
  implicit none

  integer :: x, y, z   ! Index along x-, y- and z-directions (Loop)
  character*1 :: hash  ! Hash character

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Read phase-fraction fields; Units == 1XX
  open (unit=101, file="MET.out", status="old")
  read(101,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(101,*)  met_g(x,y,z)          ! Read metal phase fraction
        end do
     end do
  end do
  close(101)

  metal_amount = sum(met_g)

  open (unit=102, file="MKW.out", status="old")
  read(102,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(102,*)  mkw_g(x,y,z)         ! Read mackinawite phase fraction
        end do
     end do
  end do
  close(102)

  open (unit=103, file="PHT.out", status="old")
  read(103,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(103,*)  pht_g(x,y,z)         ! Read pyrrhotite phase fraction
        end do
     end do
  end do
  close(103)

  open (unit=104, file="PYR.out", status="old")
  read(104,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(104,*)  pyr_g(x,y,z)         ! Read pyrite phase fraction
        end do
     end do
  end do
  close(104)

  open (unit=105, file="ENV.out", status="old")
  read(105,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(105,*)  env_g(x,y,z)        ! Read environment phase fraction
        end do
     end do
  end do
  close(105)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Read chemical-potential field; Units == 2XX
  open (unit=200, file="MUS.out", status="old")
  read(200,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(200,*)  mu_g(x,y,z)         ! Read chemical potential field
        end do
     end do
  end do
  close(200)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Read orientation field; Units == 3XX
  open (unit=303, file="OPYR.out", status="old")
  read(303,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(303,*)  opyr_g(x,y,z)       ! Read orientation field
        end do
     end do
  end do
  close(303)

  avg_mu_env = mus_mkw_pht_eqb - (R*T*2.5d0)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Read pH field; Units == 4XX
  open (unit=400, file="MUS.out", status="old")
  read(400,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           read(400,*)  pH_g(x,y,z)         ! Read pH field
        end do
     end do
  end do
  close(400)

end subroutine read_geometry


