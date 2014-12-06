subroutine pre_equilibrate()
  use commondata
  use field_evolution
  use thermo_constants
  implicit none

  integer :: x, y, z   ! Loop variables
  real*8 :: ms_diff
  integer :: i
  real*8, dimension(psx,psy,psz+2) :: newmu

  newmu = 0.0d0


  ms_diff = 2.0d0

  write(6,*) 'Pre-equilibrating.'

  do while(ms_diff .gt. 0.25d0)

     call calc_lap()
     call musolve(0.17d0,0.00d0)

     do x = 1,psx_g
        do y = 1,psy_g
           do z = 1,psz_g
              if (pht(x,y,z).gt.0.97) then
                 mu(x,y,z) = max(min(mu(x,y,z),mus_pht_pyr_eqb),mus_met_pht_eqb)
              end if          
           end do
        end do
     end do

     ms_diff = 0.0d0

     do x = 1,psx_g
        do y = 1,psy_g
           do z = 1,psz_g
              ms_diff = ms_diff + abs(newmu(x,y,z)-oldmu(x,y,z))
           end do
        end do
     end do

     ms_diff = ms_diff/sum(pht)

  end do

  write(6,*) 'Done pre-equilibrating.'

end subroutine pre_equilibrate
