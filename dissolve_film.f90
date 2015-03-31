subroutine dissolve_film()
  use commondata
  use fields
  use thermo_constants
  implicit none


  integer :: x, y, z   ! Loop variables
  real*8 :: diss_rate_met,diss_rate_mkw,diss_rate_pht,diss_rate_pyr,diss_rate_env
  real*8 :: dissolved_met,dissolved_mkw,dissolved_pht,dissolved_pyr
  real*8 :: pf_after_dissolve


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

                 pf_after_dissolve = met(x,y,z) - ((diss_rate_met*1E-9*dt)/dpf)
                 dissolved_met = met(x,y,z) - pf_after_dissolve
                 met(x,y,z) = pf_after_dissolve

                 pf_after_dissolve = mkw(x,y,z) - ((diss_rate_mkw*1E-9*dt)/dpf)
                 dissolved_mkw = mkw(x,y,z) - pf_after_dissolve
                 mkw(x,y,z) = pf_after_dissolve

                 pf_after_dissolve = pht(x,y,z) - ((diss_rate_pht*1E-9*dt)/dpf)
                 dissolved_pht = pht(x,y,z) - pf_after_dissolve
                 pht(x,y,z) = pf_after_dissolve

                 pf_after_dissolve = pyr(x,y,z) - ((diss_rate_pyr*1E-9*dt)/dpf)
                 dissolved_pyr = pyr(x,y,z) - pf_after_dissolve
                 pyr(x,y,z) = pf_after_dissolve

                 env(x,y,z) = env(x,y,z) + dissolved_met + dissolved_mkw + dissolved_pht + dissolved_pyr

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

