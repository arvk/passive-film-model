subroutine pre_equilibrate()
  use commondata
  use laplacians
  use field_evolution
  use thermo_constants
  implicit none

  integer :: x, y, z   ! Loop variables
  real*8 :: ms_diff
  real*8 :: dt
  integer :: i
  real*8, dimension(psx,psy,psz+2) :: my_newmu, my_oldmu
  integer :: wrap

  !! Diffusivities
  real*8 :: D_S_met, D_S_env, D_S_pht
  real*8 :: D_Fe_met, D_Fe_env, D_Fe_pht
  real*8 :: D_inter_met, D_inter_env, D_inter_pht, D
  real*8 :: D_H_met, D_H_env, D_H_pht, D_H

  dt = 0.00000001d0

  my_newmu = 0.0d0
  my_oldmu = 0.0d0

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           my_oldmu(x,y,z) = mu(x,y,z)
           my_newmu(x,y,z) = mu(x,y,z)
        end do
     end do
  end do

  ms_diff = 200.0d0

!do while(ms_diff .gt. 10.00d0)

do i = 1,10

     call swap_mu()

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

              D = 1E-10

              dmu_dt(x,y,z) = (D*del2mu(x,y,z))

           end do
        end do
     end do

     do x = 1,psx
        do y = 1,psy
           do z = 2,psz+1
              my_newmu(x,y,z) = my_newmu(x,y,z) + (dt*dmu_dt(x,y,z))
              my_newmu(x,y,z) = max(min(my_newmu(x,y,z),max_mu),min_mu)
           end do
        end do
     end do


     do x = 1,psx
        do y = 1,psy
           do z = 2,psz+1
              if (env(x,y,z).gt.0.97) then
                 env(x,y,z) = 1.0d0
                 met(x,y,z) = 0.0d0
                 pht(x,y,z) = 0.0d0
                 my_newmu(x,y,z) = avg_mu_env
              end if
              if (met(x,y,z).gt.0.97) then
                 my_newmu(x,y,z) = avg_mu_met
              end if
           end do
        end do
     end do


     ms_diff = 0.0d0

     do x = 1,psx
        do y = 1,psy
           do z = 2,psz+1
              ms_diff = ms_diff + abs(my_newmu(x,y,z)-my_oldmu(x,y,z))
           end do
        end do
     end do

     ms_diff = ms_diff/(psx*psy*psz)

     do x = 1,psx
        do y = 1,psy
           do z = 2,psz+1
              my_oldmu(x,y,z)=my_newmu(x,y,z)
           end do
        end do
     end do

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           mu(x,y,z) = my_newmu(x,y,z)
        end do
     end do
  end do

  end do

  ! call gather_pf()

  ! if (rank.eq.0) then
  !    call system("rm -rf MYMU.out")
  !    open (unit=1011, file="MYMU.out", status="new")
  !    write(1011,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  !    do z = 1,psz_g
  !       write(1011,*)  mu_g(3,3,z)
  !    end do
  !    close(1011)  
  ! end if

  !  write(6,*) 'Done Pre-Eqbing'




end subroutine pre_equilibrate

