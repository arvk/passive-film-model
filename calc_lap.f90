subroutine calc_lap
  use commondata
  use laplacians
  implicit none

  integer :: x, y, z   ! Loop variables
  integer :: wrap

!!!! Initialize all variables to zero
  del2pht = 0.0d0 
  del2env = 0.0d0 
  del2met = 0.0d0 
  del2mu = 0.0d0 
  del2ph = 0.0d0 


!!!!! Calculating laplacian in the bulk

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           del2met(x,y,z) = 14*(met(wrap(x+1,psx),y,z)+met(wrap(x-1,psx),y,z)+met(x,wrap(y+1,psy),z)&
&+met(x,wrap(y-1,psy),z)+met(x,y,z+1)+met(x,y,z-1))
           del2met(x,y,z) = del2met(x,y,z) + 3*(met(wrap(x+1,psx),wrap(y+1,psy),z)+met(wrap(x-1,psx),wrap(y+1,psy),z)+&
&met(wrap(x+1,psx),wrap(y-1,psy),z)+met(wrap(x-1,psx),wrap(y-1,psy),z))
           del2met(x,y,z) = del2met(x,y,z) + 3*(met(x,wrap(y+1,psy),z+1)+met(x,wrap(y+1,psy),z-1)+&
&met(x,wrap(y-1,psy),z+1)+met(x,wrap(y-1,psy),z-1))
           del2met(x,y,z) = del2met(x,y,z) + 3*(met(wrap(x+1,psx),y,z+1)+met(wrap(x-1,psx),y,z+1)+&
&met(wrap(x+1,psx),y,z-1)+met(wrap(x-1,psx),y,z-1))
           del2met(x,y,z) = del2met(x,y,z) + 1*(met(wrap(x+1,psx),wrap(y+1,psy),z+1)+met(wrap(x-1,psx),&
wrap(y+1,psy),z+1)+met(wrap(x+1,psx),wrap(y-1,psy),z+1)+met(wrap(x-1,psx),wrap(y-1,psy),z+1))
           del2met(x,y,z) = del2met(x,y,z) + 1*(met(wrap(x+1,psx),wrap(y+1,psy),z-1)+met(wrap(x-1,psx),&
wrap(y+1,psy),z-1)+met(wrap(x+1,psx),wrap(y-1,psy),z-1)+met(wrap(x-1,psx),wrap(y-1,psy),z-1))
           del2met(x,y,z) = (del2met(x,y,z) - 128*met(x,y,z))/(30*dpf*dpf)

           del2pht(x,y,z) = 14*(pht(wrap(x+1,psx),y,z)+pht(wrap(x-1,psx),y,z)+pht(x,wrap(y+1,psy),z)+&
&pht(x,wrap(y-1,psy),z)+pht(x,y,z+1)+pht(x,y,z-1))
           del2pht(x,y,z) = del2pht(x,y,z) + 3*(pht(wrap(x+1,psx),wrap(y+1,psy),z)+pht(wrap(x-1,psx),wrap(y+1,psy),z)+&
&pht(wrap(x+1,psx),wrap(y-1,psy),z)+pht(wrap(x-1,psx),wrap(y-1,psy),z))
           del2pht(x,y,z) = del2pht(x,y,z) + 3*(pht(x,wrap(y+1,psy),z+1)+pht(x,wrap(y+1,psy),z-1)+&
&pht(x,wrap(y-1,psy),z+1)+pht(x,wrap(y-1,psy),z-1))
           del2pht(x,y,z) = del2pht(x,y,z) + 3*(pht(wrap(x+1,psx),y,z+1)+pht(wrap(x-1,psx),y,z+1)+&
&pht(wrap(x+1,psx),y,z-1)+pht(wrap(x-1,psx),y,z-1))
           del2pht(x,y,z) = del2pht(x,y,z) + 1*(pht(wrap(x+1,psx),wrap(y+1,psy),z+1)+pht(wrap(x-1,psx),&
&wrap(y+1,psy),z+1)+pht(wrap(x+1,psx),wrap(y-1,psy),z+1)+pht(wrap(x-1,psx),wrap(y-1,psy),z+1))
           del2pht(x,y,z) = del2pht(x,y,z) + 1*(pht(wrap(x+1,psx),wrap(y+1,psy),z-1)+pht(wrap(x-1,psx),&
&wrap(y+1,psy),z-1)+pht(wrap(x+1,psx),wrap(y-1,psy),z-1)+pht(wrap(x-1,psx),wrap(y-1,psy),z-1))
           del2pht(x,y,z) = (del2pht(x,y,z) - 128*pht(x,y,z))/(30*dpf*dpf)

           del2env(x,y,z) = 14*(env(wrap(x+1,psx),y,z)+env(wrap(x-1,psx),y,z)+env(x,wrap(y+1,psy),z)+&
&env(x,wrap(y-1,psy),z)+env(x,y,z+1)+env(x,y,z-1))
           del2env(x,y,z) = del2env(x,y,z) + 3*(env(wrap(x+1,psx),wrap(y+1,psy),z)+env(wrap(x-1,psx),wrap(y+1,psy),z)+&
&env(wrap(x+1,psx),wrap(y-1,psy),z)+env(wrap(x-1,psx),wrap(y-1,psy),z))
           del2env(x,y,z) = del2env(x,y,z) + 3*(env(x,wrap(y+1,psy),z+1)+env(x,wrap(y+1,psy),z-1)+&
&env(x,wrap(y-1,psy),z+1)+env(x,wrap(y-1,psy),z-1))
           del2env(x,y,z) = del2env(x,y,z) + 3*(env(wrap(x+1,psx),y,z+1)+env(wrap(x-1,psx),y,z+1)+&
&env(wrap(x+1,psx),y,z-1)+env(wrap(x-1,psx),y,z-1))
           del2env(x,y,z) = del2env(x,y,z) + 1*(env(wrap(x+1,psx),wrap(y+1,psy),z+1)+env(wrap(x-1,psx),&
&wrap(y+1,psy),z+1)+env(wrap(x+1,psx),wrap(y-1,psy),z+1)+env(wrap(x-1,psx),wrap(y-1,psy),z+1))
           del2env(x,y,z) = del2env(x,y,z) + 1*(env(wrap(x+1,psx),wrap(y+1,psy),z-1)+env(wrap(x-1,psx),&
&wrap(y+1,psy),z-1)+env(wrap(x+1,psx),wrap(y-1,psy),z-1)+env(wrap(x-1,psx),wrap(y-1,psy),z-1))
           del2env(x,y,z) = (del2env(x,y,z) - 128*env(x,y,z))/(30*dpf*dpf)

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

           del2ph(x,y,z) = 14*(ph(wrap(x+1,psx),y,z)+ph(wrap(x-1,psx),y,z)+ph(x,wrap(y+1,psy),z)+&
&ph(x,wrap(y-1,psy),z)+ph(x,y,z+1)+ph(x,y,z-1))
           del2ph(x,y,z) = del2ph(x,y,z) + 3*(ph(wrap(x+1,psx),wrap(y+1,psy),z)+ph(wrap(x-1,psx),wrap(y+1,psy),z)+&
&ph(wrap(x+1,psx),wrap(y-1,psy),z)+ph(wrap(x-1,psx),wrap(y-1,psy),z))
           del2ph(x,y,z) = del2ph(x,y,z) + 3*(ph(x,wrap(y+1,psy),z+1)+ph(x,wrap(y+1,psy),z-1)+&
&ph(x,wrap(y-1,psy),z+1)+ph(x,wrap(y-1,psy),z-1))
           del2ph(x,y,z) = del2ph(x,y,z) + 3*(ph(wrap(x+1,psx),y,z+1)+ph(wrap(x-1,psx),y,z+1)+&
&ph(wrap(x+1,psx),y,z-1)+ph(wrap(x-1,psx),y,z-1))
           del2ph(x,y,z) = del2ph(x,y,z) + 1*(ph(wrap(x+1,psx),wrap(y+1,psy),z+1)+ph(wrap(x-1,psx),&
&wrap(y+1,psy),z+1)+ph(wrap(x+1,psx),wrap(y-1,psy),z+1)+ph(wrap(x-1,psx),wrap(y-1,psy),z+1))
           del2ph(x,y,z) = del2ph(x,y,z) + 1*(ph(wrap(x+1,psx),wrap(y+1,psy),z-1)+ph(wrap(x-1,psx),&
&wrap(y+1,psy),z-1)+ph(wrap(x+1,psx),wrap(y-1,psy),z-1)+ph(wrap(x-1,psx),wrap(y-1,psy),z-1))
           del2ph(x,y,z) = (del2ph(x,y,z) - 128*ph(x,y,z))/(30*dpf*dpf)

        end do
     end do
  end do


!!! Impose boundary conditions on the Laplacian

  do x = 1,psx
     do y = 1,psy
        del2pht(x,y,1) = 0.0d0 ; del2pht(x,y,psz+2) = 0.0d0
        del2env(x,y,1) = 0.0d0 ; del2env(x,y,psz+2) = 0.0d0
        del2met(x,y,1) = 0.0d0 ; del2met(x,y,psz+2) = 0.0d0
        del2mu(x,y,1) = 0.0d0 ; del2mu(x,y,psz+2) = 0.0d0
        del2ph(x,y,1) = 0.0d0 ; del2ph(x,y,psz+2) = 0.0d0
     end do
  end do

end subroutine calc_lap




integer function wrap(a,lim)
  implicit none
  integer, intent(in) :: a,lim
  wrap = modulo(a-1,lim)+1
  return
end function wrap
