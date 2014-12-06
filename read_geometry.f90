subroutine read_geometry
  use commondata
  use laplacians

  implicit none

  integer :: x, y, z   ! Loop variables
  character*1 :: hash

  open (unit=101, file="PHT.out", status="old")
  read(101,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
    do y = 1,psy_g
       read(101,*)  (pht(x,y,z),z=1,psz_g) 
    end do
 end do
  close(101)  

  open (unit=102, file="ENV.out", status="old")
  read(102,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
    do y = 1,psy_g
       read(102,*)  (env(x,y,z),z=1,psz_g) 
    end do
 end do
  close(102)  

  open (unit=103, file="MET.out", status="old")
  read(103,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
    do y = 1,psy_g
       read(103,*)  (met(x,y,z),z=1,psz_g) 
    end do
 end do
 close(103)  

  open (unit=105, file="MUS.out", status="old")
  read(105,*) hash, psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
   do y = 1,psy_g
      read(105,*)  (mu(x,y,z),z=1,psz_g) 
    end do
 end do
 close(105)  

end subroutine read_geometry


