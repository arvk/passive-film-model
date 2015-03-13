subroutine dissolve_film()
  use commondata
  use fields
  use thermo_constants
  implicit none


  integer :: x, y, z   ! Loop variables
  real*8 :: diss_rate_met,diss_rate_mkw,diss_rate_pht,diss_rate_pyr,diss_rate_env


!!!! Define dissolution rates (in nm/s)
  diss_rate_met = 1.5d0
  diss_rate_mkw = 0.1d0
  diss_rate_pht = 0.1d0
  diss_rate_pyr = 0.1d0
  diss_rate_env = 0.1d0



     do x = 1,psx
        do y = 1,psy
           do z = psz+1,2,-1
              if ((env(x,y,z) .le. 9.9E-1).and.(env(x,y,z+1) .gt. 9.9E-1)) then

                 met(x,y,z) = max(met(x,y,z) - ((diss_rate_met*1E-9*dt)/dpf),0.0d0)
                 mkw(x,y,z) = max(mkw(x,y,z) - ((diss_rate_mkw*1E-9*dt)/dpf),0.0d0)
                 pht(x,y,z) = max(pht(x,y,z) - ((diss_rate_pht*1E-9*dt)/dpf),0.0d0)
                 pyr(x,y,z) = max(pyr(x,y,z) - ((diss_rate_pyr*1E-9*dt)/dpf),0.0d0)
                 env(x,y,z) = 1.0d0 - (met(x,y,z) + mkw(x,y,z) + pht(x,y,z) + pyr(x,y,z))

                 exit
              end if
           end do
        end do
     end do




end subroutine dissolve_film

