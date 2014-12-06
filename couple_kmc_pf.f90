subroutine couple_kmc_pf
  use commondata
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z  ! Loop variables
  integer :: i, j     ! Loop variables
  real*8 :: sum_kg
  real*8 :: oldpht, newpht, oldenv, newenv
  real*8 :: coupling_const
  integer, dimension(psx_g,psy_g) :: interface_loc

  interface_loc = 0

  do x = 1,psx_g
     do y = 1,psy_g
        do z = psz_g-1,1,-1
           if ((pht_g(x,y,z) .gt. 1.0E-3).and.(pht_g(x,y,z+1) .lt. 1.0E-3)) then
              interface_loc(x,y) = z
              exit
           end if
        end do
     end do
  end do

  do x = 1,psx_g
     do y = 1,psy_g

        sum_kg = 0.0d0

        do i = 1,kg_scale
           do j = 1,kg_scale
              sum_kg = sum_kg + kg_g(((x-1)*kg_scale)+i,((y-1)*kg_scale)+j)
           end do
        end do

        oldpht = pht_g(x,y,interface_loc(x,y))
        oldenv = env_g(x,y,interface_loc(x,y))

        newpht = pht_g(x,y,interface_loc(x,y))+(0.2d0*((sum_kg/(kg_scale*kg_scale))-pht_g(x,y,interface_loc(x,y))))
        newpht = max(min(newpht,1.0d0),0.0d0)

        env_g(x,y,interface_loc(x,y)) = env_g(x,y,interface_loc(x,y)) - (newpht-pht_g(x,y,interface_loc(x,y)))
        pht_g(x,y,interface_loc(x,y)) = pht_g(x,y,interface_loc(x,y)) + (newpht-pht_g(x,y,interface_loc(x,y)))

        newpht = (pht_g(x,y,interface_loc(x,y))/(pht_g(x,y,interface_loc(x,y))+env_g(x,y,interface_loc(x,y))))*(oldpht+oldenv)
        newenv = (env_g(x,y,interface_loc(x,y))/(pht_g(x,y,interface_loc(x,y))+env_g(x,y,interface_loc(x,y))))*(oldpht+oldenv)

        pht_g(x,y,interface_loc(x,y)) = newpht
        env_g(x,y,interface_loc(x,y)) = newenv

        ph_g(x,y,interface_loc(x,y)) = ph_g(x,y,interface_loc(x,y)) - ((newpht-oldpht)*7.0d0)
        ph_g(x,y,interface_loc(x,y)) = max(min(ph_g(x,y,interface_loc(x,y)),14.0d0),0.0d0)

        if (newpht.gt.oldpht) then
           mu_g(x,y,interface_loc(x,y)) = mu_g(x,y,interface_loc(x,y)) + ((newpht-oldpht)*(avg_mu_pht-mu_g(x,y,interface_loc(x,y))))
        else
           mu_g(x,y,interface_loc(x,y)) = mu_g(x,y,interface_loc(x,y)) - ((newpht-oldpht)*(avg_mu_env-mu_g(x,y,interface_loc(x,y))))
        end if

     end do
  end do



end subroutine couple_kmc_pf
