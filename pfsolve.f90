subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: ierr,status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter

  !!---------------PF evolution-----------------!!

  !! Bulk free energy 
  real*8 :: f_met, f_mkw, f_pht, f_pyr, f_env
  real*8 :: w_met, w_mkw, w_pht, w_pyr, w_env

  !! Derivative of bulk free energy with phase field
  real*8 :: dF_dmet_mkw, dF_dmet_pht, dF_dmet_pyr, dF_dmet_env
  real*8 :: dF_dmkw_met, dF_dmkw_pht, dF_dmkw_pyr, dF_dmkw_env
  real*8 :: dF_dpht_met, dF_dpht_mkw, dF_dpht_pyr, dF_dpht_env
  real*8 :: dF_dpyr_met, dF_dpyr_mkw, dF_dpyr_pht, dF_dpyr_env
  real*8 :: dF_denv_met, dF_denv_mkw, dF_denv_pht, dF_denv_pyr

  real*8 :: sigma_pyr_met, sigma_pyr_mkw, sigma_pyr_pht, sigma_pyr_env
  real*8 :: sigma_met_pyr, sigma_mkw_pyr, sigma_pht_pyr, sigma_env_pyr

  integer :: int_count

  real*8, dimension(psx,psy,psz+2) :: del_opyr
  real*8 :: odiff
  integer :: wrap
  real*8 :: delx,dely,delz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize field variables to 0
 dmet_dt = 0.0d0 ; dmkw_dt = 0.0d0 ; dpht_dt = 0.0d0 ; dpyr_dt = 0.0d0 ; denv_dt = 0.0d0

     if (mod(iter,swap_freq_pf).eq.1) then
        call swap_pf()
     end if

  !! Calculate laplacian of phase/composition fields
  call calc_lap_pf()
  call calc_grad_pf()



  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           delx = odiff(opyr(wrap(x+1,psx),y,z),opyr(wrap(x-1,psx),y,z))
           dely = odiff(opyr(x,wrap(y+1,psy),z),opyr(x,wrap(y-1,psy),z))
           delz = odiff(opyr(x,y,z+1),opyr(x,y,z-1))

           del_opyr(x,y,z) = sqrt((delx*delx)+(dely*dely)+(delz*delz))/dpf

        end do
     end do
  end do



  !################################################################
  !#######################--PF EVOLUTION--#########################
  !################################################################

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0



           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_mkw = sigma_pyr_mkw_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_mkw_pyr = sigma_pyr_mkw
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env




!~~~~~~~~MET~~~~~~~~!

           dF_dmet_mkw = (sigma_met_mkw*((met(x,y,z)*del2mkw(x,y,z))-(mkw(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_mkw*met(x,y,z)*mkw(x,y,z)*(met(x,y,z)-mkw(x,y,z))) - &
                & (2*(met(x,y,z)+mkw(x,y,z))*(met(x,y,z)+mkw(x,y,z))*(w_met-w_mkw))

           dF_dmet_pht = (sigma_met_pht*((met(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_pht*met(x,y,z)*pht(x,y,z)*(met(x,y,z)-pht(x,y,z))) - &
                & (2*(met(x,y,z)+pht(x,y,z))*(met(x,y,z)+pht(x,y,z))*(w_met-w_pht))

           dF_dmet_env = (sigma_met_env*((met(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_env*met(x,y,z)*env(x,y,z)*(met(x,y,z)-env(x,y,z))) - &
                & (2*(met(x,y,z)+env(x,y,z))*(met(x,y,z)+env(x,y,z))*(w_met-w_env))

           dF_dmet_pyr = (sigma_met_pyr*((met(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_pyr*met(x,y,z)*pyr(x,y,z)*(met(x,y,z)-pyr(x,y,z))) - &
                & (2*(met(x,y,z)+pyr(x,y,z))*(met(x,y,z)+pyr(x,y,z))*(w_met-w_pyr))

!~~~~~~~~MKW~~~~~~~~!

           dF_dmkw_met = (sigma_mkw_met*((mkw(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2mkw(x,y,z)))) + &
                & (2*hill_mkw_met*mkw(x,y,z)*met(x,y,z)*(mkw(x,y,z)-met(x,y,z))) - &
                & (2*(mkw(x,y,z)+met(x,y,z))*(mkw(x,y,z)+met(x,y,z))*(w_mkw-w_met))

           dF_dmkw_pht = (sigma_mkw_pht*((mkw(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2mkw(x,y,z)))) + &
                & (2*hill_mkw_pht*mkw(x,y,z)*pht(x,y,z)*(mkw(x,y,z)-pht(x,y,z))) - &
                & (2*(mkw(x,y,z)+pht(x,y,z))*(mkw(x,y,z)+pht(x,y,z))*(w_mkw-w_pht))

           dF_dmkw_env = (sigma_mkw_env*((mkw(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2mkw(x,y,z)))) + &
                & (2*hill_mkw_env*mkw(x,y,z)*env(x,y,z)*(mkw(x,y,z)-env(x,y,z))) - &
                & (2*(mkw(x,y,z)+env(x,y,z))*(mkw(x,y,z)+env(x,y,z))*(w_mkw-w_env))

           dF_dmkw_pyr = (sigma_mkw_pyr*((mkw(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2mkw(x,y,z)))) + &
                & (2*hill_mkw_pyr*mkw(x,y,z)*pyr(x,y,z)*(mkw(x,y,z)-pyr(x,y,z))) - &
                & (2*(mkw(x,y,z)+pyr(x,y,z))*(mkw(x,y,z)+pyr(x,y,z))*(w_mkw-w_pyr))

!~~~~~~~~PHT~~~~~~~~!

           dF_dpht_met = (sigma_pht_met*((pht(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_met*pht(x,y,z)*met(x,y,z)*(pht(x,y,z)-met(x,y,z))) - &
                & (2*(pht(x,y,z)+met(x,y,z))*(pht(x,y,z)+met(x,y,z))*(w_pht-w_met))

           dF_dpht_mkw = (sigma_pht_mkw*((pht(x,y,z)*del2mkw(x,y,z))-(mkw(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_mkw*pht(x,y,z)*mkw(x,y,z)*(pht(x,y,z)-mkw(x,y,z))) - &
                & (2*(pht(x,y,z)+mkw(x,y,z))*(pht(x,y,z)+mkw(x,y,z))*(w_pht-w_mkw))

           dF_dpht_env = (sigma_pht_env*((pht(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_env*pht(x,y,z)*env(x,y,z)*(pht(x,y,z)-env(x,y,z))) - &
                & (2*(pht(x,y,z)+env(x,y,z))*(pht(x,y,z)+env(x,y,z))*(w_pht-w_env))

           dF_dpht_pyr = (sigma_pht_pyr*((pht(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_pyr*pht(x,y,z)*pyr(x,y,z)*(pht(x,y,z)-pyr(x,y,z))) - &
                & (2*(pht(x,y,z)+pyr(x,y,z))*(pht(x,y,z)+pyr(x,y,z))*(w_pht-w_pyr))

!~~~~~~~~PYR~~~~~~~~!

           dF_dpyr_met = (sigma_pyr_met*((pyr(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_met*pyr(x,y,z)*met(x,y,z)*(pyr(x,y,z)-met(x,y,z))) - &
                & (2*(pyr(x,y,z)+met(x,y,z))*(pyr(x,y,z)+met(x,y,z))*(w_pyr-w_met))

           dF_dpyr_mkw = (sigma_pyr_mkw*((pyr(x,y,z)*del2mkw(x,y,z))-(mkw(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_mkw*pyr(x,y,z)*mkw(x,y,z)*(pyr(x,y,z)-mkw(x,y,z))) - &
                & (2*(pyr(x,y,z)+mkw(x,y,z))*(pyr(x,y,z)+mkw(x,y,z))*(w_pyr-w_mkw))

           dF_dpyr_pht = (sigma_pyr_pht*((pyr(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_pht*pyr(x,y,z)*pht(x,y,z)*(pyr(x,y,z)-pht(x,y,z))) - &
                & (2*(pyr(x,y,z)+pht(x,y,z))*(pyr(x,y,z)+pht(x,y,z))*(w_pyr-w_pht))

           dF_dpyr_env = (sigma_pyr_env*((pyr(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_env*pyr(x,y,z)*env(x,y,z)*(pyr(x,y,z)-env(x,y,z))) - &
                & (2*(pyr(x,y,z)+env(x,y,z))*(pyr(x,y,z)+env(x,y,z))*(w_pyr-w_env))

!~~~~~~~~PYR~~~~~~~~!

           dF_denv_met = (sigma_env_met*((env(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_met*env(x,y,z)*met(x,y,z)*(env(x,y,z)-met(x,y,z))) - &
                & (2*(env(x,y,z)+met(x,y,z))*(env(x,y,z)+met(x,y,z))*(w_env-w_met))

           dF_denv_mkw = (sigma_env_mkw*((env(x,y,z)*del2mkw(x,y,z))-(mkw(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_mkw*env(x,y,z)*mkw(x,y,z)*(env(x,y,z)-mkw(x,y,z))) - &
                & (2*(env(x,y,z)+mkw(x,y,z))*(env(x,y,z)+mkw(x,y,z))*(w_env-w_mkw))

           dF_denv_pht = (sigma_env_pht*((env(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_pht*env(x,y,z)*pht(x,y,z)*(env(x,y,z)-pht(x,y,z))) - &
                & (2*(env(x,y,z)+pht(x,y,z))*(env(x,y,z)+pht(x,y,z))*(w_env-w_pht))

           dF_denv_pyr = (sigma_env_pyr*((env(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_pyr*env(x,y,z)*pyr(x,y,z)*(env(x,y,z)-pyr(x,y,z))) - &
                & (2*(env(x,y,z)+pyr(x,y,z))*(env(x,y,z)+pyr(x,y,z))*(w_env-w_pyr))













           !! Error correction        
           correct_rounding = 0.50d0*(dF_dmet_mkw - dF_dmkw_met)
           dF_dmet_mkw = correct_rounding ; dF_dmkw_met = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmet_pht - dF_dpht_met)
           dF_dmet_pht = correct_rounding ; dF_dpht_met = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmet_pyr - dF_dpyr_met)
           dF_dmet_pyr = correct_rounding ; dF_dpyr_met = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmet_env - dF_denv_met)
           dF_dmet_env = correct_rounding ; dF_denv_met = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmkw_pht - dF_dpht_mkw)
           dF_dmkw_pht = correct_rounding ; dF_dpht_mkw = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmkw_pyr - dF_dpyr_mkw)
           dF_dmkw_pyr = correct_rounding ; dF_dpyr_mkw = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmkw_env - dF_denv_mkw)
           dF_dmkw_env = correct_rounding ; dF_denv_mkw = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpht_pyr - dF_dpyr_pht)
           dF_dpht_pyr = correct_rounding ; dF_dpyr_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpht_env - dF_denv_pht)
           dF_dpht_env = correct_rounding ; dF_denv_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpyr_env - dF_denv_pyr)
           dF_dpyr_env = correct_rounding ; dF_denv_pyr = 0.0d0 - correct_rounding







           !! Calculate phase field(s) update
           dmet_dt(x,y,z) = (M_met_mkw*dF_dmet_mkw)+(M_met_pht*dF_dmet_pht)+(M_met_pyr*dF_dmet_pyr)+(M_met_env*dF_dmet_env)
           dmkw_dt(x,y,z) = (M_mkw_met*dF_dmkw_met)+(M_mkw_pht*dF_dmkw_pht)+(M_mkw_pyr*dF_dmkw_pyr)+(M_mkw_env*dF_dmkw_env)
           dpht_dt(x,y,z) = (M_pht_met*dF_dpht_met)+(M_pht_mkw*dF_dpht_mkw)+(M_pht_pyr*dF_dpht_pyr)+(M_pht_env*dF_dpht_env)
           dpyr_dt(x,y,z) = (M_pyr_met*dF_dpyr_met)+(M_pyr_mkw*dF_dpyr_mkw)+(M_pyr_pht*dF_dpyr_pht)+(M_pyr_env*dF_dpyr_env)
           denv_dt(x,y,z) = (M_env_met*dF_denv_met)+(M_env_mkw*dF_denv_mkw)+(M_env_pht*dF_denv_pht)+(M_env_pyr*dF_denv_pyr)
           dpyr_dt(x,y,z) = dpyr_dt(x,y,z) - (25.0d0*pyr(x,y,z)*((M_pyr_met+M_pyr_mkw+M_pyr_pht+M_pyr_env)*del_opyr(x,y,z)))
        end do
     end do
  end do



  !! Apply boundary conditions to phase field(s) update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           dmet_dt(x,y,2) = 0.0d0
           dmkw_dt(x,y,2) = 0.0d0
           dpht_dt(x,y,2) = 0.0d0
           dpyr_dt(x,y,2) = 0.0d0
           denv_dt(x,y,2) = 0.0d0
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           dmet_dt(x,y,psz+1) = 0.0d0
           dmkw_dt(x,y,psz+1) = 0.0d0
           dpht_dt(x,y,psz+1) = 0.0d0
           dpyr_dt(x,y,psz+1) = 0.0d0
           denv_dt(x,y,psz+1) = 0.0d0
        end do
     end do
  end if


  !! Update phase fields
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmet(x,y,z) = max(min((met(x,y,z) + (dt*dmet_dt(x,y,z))),1.0d0),0.0d0)
           newmkw(x,y,z) = max(min((mkw(x,y,z) + (dt*dmkw_dt(x,y,z))),1.0d0),0.0d0)
           newpht(x,y,z) = max(min((pht(x,y,z) + (dt*dpht_dt(x,y,z))),1.0d0),0.0d0)
           newenv(x,y,z) = max(min((env(x,y,z) + (dt*denv_dt(x,y,z))),1.0d0),0.0d0)
           newpyr(x,y,z) = max(min((pyr(x,y,z) + (dt*dpyr_dt(x,y,z))),1.0d0),0.0d0)

           correct_rounding = (newmet(x,y,z)+newmkw(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))

           newmet(x,y,z) = newmet(x,y,z)/correct_rounding
           newmkw(x,y,z) = newmkw(x,y,z)/correct_rounding
           newpht(x,y,z) = newpht(x,y,z)/correct_rounding
           newenv(x,y,z) = newenv(x,y,z)/correct_rounding
           newpyr(x,y,z) = newpyr(x,y,z)/correct_rounding
        end do
     end do
  end do


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (env(x,y,z).gt.0.99) then
              newmet(x,y,z) = 0.0d0
              newmkw(x,y,z) = 0.0d0
              newpht(x,y,z) = 0.0d0
              newpyr(x,y,z) = 0.0d0
              newenv(x,y,z) = 1.0d0
           end if
        end do
     end do
  end do


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           dmet_dt(x,y,z) = (newmet(x,y,z) - met(x,y,z))/dt
           dmkw_dt(x,y,z) = (newmkw(x,y,z) - mkw(x,y,z))/dt
           dpht_dt(x,y,z) = (newpht(x,y,z) - pht(x,y,z))/dt
           dpyr_dt(x,y,z) = (newpyr(x,y,z) - pyr(x,y,z))/dt
           denv_dt(x,y,z) = (newenv(x,y,z) - env(x,y,z))/dt
        end do
     end do
  end do


  sulfidation_rate = 0.0d0
  int_count = 1

  do x = 1,psx
     do y = 1,psy
        do z = psz+1,2,-1
           if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then
              int_count = int_count + 1
              sulfidation_rate = sulfidation_rate + 0.04356E-9*(exp(avg_mu_env/(R*T))) !! Ref = Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
           end if
        end do
     end do
  end do

  sulfidation_rate = max((sulfidation_rate/int_count) + 0.01372E-9,0.0d0)

  !###########################################################
  !##################--UPDATE ALL FIELDS--####################
  !###########################################################

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           met(x,y,z) = newmet(x,y,z)
           mkw(x,y,z) = newmkw(x,y,z)
           pht(x,y,z) = newpht(x,y,z)
           pyr(x,y,z) = newpyr(x,y,z)
           env(x,y,z) = newenv(x,y,z)
        end do
     end do
  end do

end subroutine pfsolve

