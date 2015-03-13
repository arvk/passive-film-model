subroutine voids_create()
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none

  integer :: x, y, z   ! Loop variables

    do z = 2,psz+1
     do y = 1,psy
        do x = 1,psx
           if ((mkw(x,y,z).lt.0.3).and.(mkw(x,y,z).gt.0.02).and.(met(x,y,z).gt.0.7).and.(met(x,y,z).lt.0.98)) then
              voids(x,y,z) = 1.0d0
              write(6,*) "Made 1 void @ ",x,y,z
           end if
        end do
     end do
  end do

end subroutine voids_create
