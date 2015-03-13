subroutine couple_kmc_pf
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z  ! Loop variables
  integer :: i, j     ! Loop variables
  real*8 :: sum_kg
  real*8 :: pht_bc, env_bc ! Local phase fractions at (x,y,z) before coupling. This is read from pht(x,y,z) and env(x,y,z)
  real*8 :: pht_ac, env_ac ! Local phase fractions at (x,y,z) after coupling. This is written to pht(x,y,z) and env(x,y,z)
  real*8 :: coupling_const
  real*8 :: avg_mu_pht
  integer, dimension(psx_g,psy_g) :: interface_loc

  interface_loc = 0

  do x = 1,psx_g
     do y = 1,psy_g
        do z = psz_g-1,1,-1
           if ((pht_g(x,y,z) .gt. 5.0E-1).and.(pht_g(x,y,z+1) .lt. 5.0E-1)) then
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

        pht_bc = pht_g(x,y,interface_loc(x,y))
        env_bc = env_g(x,y,interface_loc(x,y))

        pht_ac = pht_g(x,y,interface_loc(x,y))+(0.2d0*((sum_kg/(kg_scale*kg_scale))-pht_g(x,y,interface_loc(x,y))))
        pht_ac = max(min(pht_ac,1.0d0),0.0d0)

        env_g(x,y,interface_loc(x,y)) = env_g(x,y,interface_loc(x,y)) - (pht_ac-pht_g(x,y,interface_loc(x,y)))
        pht_g(x,y,interface_loc(x,y)) = pht_g(x,y,interface_loc(x,y)) + (pht_ac-pht_g(x,y,interface_loc(x,y)))

        pht_ac = (pht_g(x,y,interface_loc(x,y))/(pht_g(x,y,interface_loc(x,y))+env_g(x,y,interface_loc(x,y))))*(pht_bc+env_bc)
        env_ac = (env_g(x,y,interface_loc(x,y))/(pht_g(x,y,interface_loc(x,y))+env_g(x,y,interface_loc(x,y))))*(pht_bc+env_bc)

        pht_g(x,y,interface_loc(x,y)) = pht_ac
        env_g(x,y,interface_loc(x,y)) = env_ac

        if (pht_ac.gt.pht_bc) then
           mu_g(x,y,interface_loc(x,y)) = mu_g(x,y,interface_loc(x,y)) + ((pht_ac-pht_bc)*(avg_mu_pht-mu_g(x,y,interface_loc(x,y))))
        else
           mu_g(x,y,interface_loc(x,y)) = mu_g(x,y,interface_loc(x,y)) - ((pht_ac-pht_bc)*(avg_mu_env-mu_g(x,y,interface_loc(x,y))))
        end if

     end do
  end do



end subroutine couple_kmc_pf
