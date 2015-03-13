subroutine calc_grad_pf
  use commondata
  use fields
  use gradients
  implicit none

  integer :: x, y, z   ! Loop variables
  integer :: wrap

!!!! Initialize all variables to zero
  delypyr = 0.0d0 
  delzpyr = 0.0d0 

!!!!! Calculating y-gradient in the bulk

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           delypyr(x,y,z) = (pyr(x,wrap(y+1,psy),z)-pyr(x,wrap(y-1,psy),z))/2.0d0
        end do
     end do
  end do

!!!!! Calculating z-gradient in the bulk
  do x = 1,psx
     do y = 1,psy
        delzpyr(x,y,1) = (pyr(x,y,2)-pyr(x,y,1))
        do z = 2,psz+1
           delzpyr(x,y,z) = (pyr(x,y,z+1)-pyr(x,y,z-1))/2.0d0
        end do
        delzpyr(x,y,psz+2) = (pyr(x,y,psz+2)-pyr(x,y,psz+1))
     end do
  end do

end subroutine calc_grad_pf
