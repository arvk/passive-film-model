subroutine dissolve_film()
  use commondata
  use fields
  use thermo_constants
  implicit none


  integer :: x, y, z   ! Loop variables
  real*8 :: diss_rate_met,diss_rate_mkw,diss_rate_pht,diss_rate_pyr,diss_rate_env


!!!! Define dissolution rates (in nm/s)
  diss_rate_met = 0.257d0 !!! REF = Electrodissolution Kinetics of Iron in Chloride Solutions by Robert J. Chin* and Ken Nobe, JECS Vol 119, No. 11, Nov. 1972
                          !!! Full explression is diss_rate_met = 0.257d0*((0.0d0-log10(14-pH))**0.6)*exp((0.85*F/RT)*(V+0.45))
  diss_rate_mkw = 0.015d0 !!! REF = CORROSION MECHANISMS AND MATERIAL PERFORMANCE IN ENVIRONMENTS CONTAINING HYDROGEN SULFIDE AND ELEMENTAL SULFUR, Liane Smith* and Bruce Craig, SACNUC Workshop 22nd and 23rd October, 2008, Brussels and References therein
  diss_rate_pht = 0.015d0
  diss_rate_pyr = 0.015d0
  diss_rate_env = 0.015d0



     do x = 1,psx
        do y = 1,psy
           do z = psz+(2*ghost_width)-1,2,-1
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


!! Void filling
     do x = 1,psx
        do y = 1,psy
           do z = psz+(2*ghost_width)-1,2,-1
              if ((voids(x,y,z).gt.0.0d0).and.(env(x,y,z+1).eq.1.0d0)) then
                 met(x,y,z) = 0.0d0
                 met(x,y,z) = 0.0d0
                 met(x,y,z) = 0.0d0
                 met(x,y,z) = 0.0d0
                 env(x,y,z) = 1.0d0
                 mu(x,y,z) = avg_mu_env
                 voids(x,y,z) = 0.0d0
              end if
           end do
        end do
     end do





end subroutine dissolve_film

