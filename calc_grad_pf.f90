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

!!!!! Calculating laplacian in the bulk

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           delypyr(x,y,z) = (pyr(x,wrap(y+1,psy),z)-pyr(x,wrap(y-1,psy),z))/(2*dpf)
           delzpyr(x,y,z) = (pyr(x,y,z+1)-pyr(x,y,z-1))/(2*dpf)

        end do
     end do
  end do


!!! Impose boundary conditions on the Laplacian

  do x = 1,psx
     do y = 1,psy
        delypyr(x,y,1) = 0.0d0 ; delypyr(x,y,psz+2) = 0.0d0
     end do
  end do

end subroutine calc_grad_pf
