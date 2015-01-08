subroutine calc_lap_mu
  use commondata
  use fields
  use laplacians
  implicit none

  integer :: x, y, z   ! Loop variables
  integer :: wrap

!!!! Initialize all variables to zero
  del2mu = 0.0d0 

!!!!! Calculating laplacian in the bulk

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           del2mu(x,y,z) = 14*(mu(wrap(x+1,psx),y,z)+mu(wrap(x-1,psx),y,z)+mu(x,wrap(y+1,psy),z)+&
&mu(x,wrap(y-1,psy),z)+mu(x,y,z+1)+mu(x,y,z-1))
           del2mu(x,y,z) = del2mu(x,y,z) + 3*(mu(wrap(x+1,psx),wrap(y+1,psy),z)+mu(wrap(x-1,psx),wrap(y+1,psy),z)+&
&mu(wrap(x+1,psx),wrap(y-1,psy),z)+mu(wrap(x-1,psx),wrap(y-1,psy),z))
           del2mu(x,y,z) = del2mu(x,y,z) + 3*(mu(x,wrap(y+1,psy),z+1)+mu(x,wrap(y+1,psy),z-1)+&
&mu(x,wrap(y-1,psy),z+1)+mu(x,wrap(y-1,psy),z-1))
           del2mu(x,y,z) = del2mu(x,y,z) + 3*(mu(wrap(x+1,psx),y,z+1)+mu(wrap(x-1,psx),y,z+1)+&
&mu(wrap(x+1,psx),y,z-1)+mu(wrap(x-1,psx),y,z-1))
           del2mu(x,y,z) = del2mu(x,y,z) + 1*(mu(wrap(x+1,psx),wrap(y+1,psy),z+1)+mu(wrap(x-1,psx),&
&wrap(y+1,psy),z+1)+mu(wrap(x+1,psx),wrap(y-1,psy),z+1)+mu(wrap(x-1,psx),wrap(y-1,psy),z+1))
           del2mu(x,y,z) = del2mu(x,y,z) + 1*(mu(wrap(x+1,psx),wrap(y+1,psy),z-1)+mu(wrap(x-1,psx),&
&wrap(y+1,psy),z-1)+mu(wrap(x+1,psx),wrap(y-1,psy),z-1)+mu(wrap(x-1,psx),wrap(y-1,psy),z-1))
           del2mu(x,y,z) = (del2mu(x,y,z) - 128*mu(x,y,z))/(30*dpf*dpf)

        end do
     end do
  end do


!!! Impose boundary conditions on the Laplacian

  do x = 1,psx
     do y = 1,psy
        del2mu(x,y,1) = 0.0d0 ; del2mu(x,y,psz+2) = 0.0d0
     end do
  end do

end subroutine calc_lap_mu
