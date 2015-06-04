subroutine calc_grad_pf
  use commondata
  use fields
  use gradients
  implicit none

  integer :: x, y, z ! Index for x-, y- and z-direction (Loop)
  integer :: wrap    ! Wrap index along the x-, y-, and z-direction

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  delypyr = 0.0d0 ! Initialize y-gradient to zero
  delzpyr = 0.0d0 ! Initialize z-gradient to zero

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+(2*ghost_width)
           delypyr(x,y,z) = (pyr(x,wrap(y+1,psy),z)-pyr(x,wrap(y-1,psy),z))/2.0d0   ! Calculating y-gradient in the bulk
        end do
     end do
  end do

  do x = 1,psx
     do y = 1,psy
        delzpyr(x,y,1) = (pyr(x,y,2)-pyr(x,y,1))
        do z = 2,psz+(2*ghost_width)-2
           delzpyr(x,y,z) = (pyr(x,y,z+1)-pyr(x,y,z-1))/2.0d0                       ! Calculating z-gradient in the bulk
        end do
        delzpyr(x,y,psz+(2*ghost_width)) = (pyr(x,y,psz+(2*ghost_width))-pyr(x,y,psz+(2*ghost_width)-1))
     end do
  end do

end subroutine calc_grad_pf
