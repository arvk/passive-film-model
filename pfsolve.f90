subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use field_evolution
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: ierr,status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter

  !!---------------PF evolution-----------------!!

  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr

  !! Derivative of bulk free energy with phase field
  real*8 :: dF_dpht_met, dF_dpht_env, dF_dpht_pyr
  real*8 :: dF_dmet_pht, dF_dmet_env, dF_dmet_pyr
  real*8 :: dF_denv_met, dF_denv_pht, dF_denv_pyr
  real*8 :: dF_dpyr_met, dF_dpyr_pht, dF_dpyr_env

  integer :: int_count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize field variables to 0
  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dpyr_dt = 0.0d0 ; dmu_dt = 0.0d0 

     if (mod(iter,swap_freq_pf).eq.1) then
        call swap_pf()
     end if

  !! Calculate laplacian of phase/composition fields
  call calc_lap_pf()


  !################################################################
  !#######################--PF EVOLUTION--#########################
  !################################################################

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)


           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401

           ! if (exp(mu(x,y,z)/(R*T)) .lt. 0.0015d0) then
           !    w_met = f_met - (mu(x,y,z)*exp(mu(x,y,z)/(R*T))*140401)
           ! else
           !    w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           ! end if

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)


           dF_dpht_met = (sigma_pht_met*((pht(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_met*pht(x,y,z)*met(x,y,z)*(pht(x,y,z)-met(x,y,z))) - &
                & (2*(pht(x,y,z)+met(x,y,z))*(pht(x,y,z)+met(x,y,z))*(w_pht-w_met))

           dF_dpht_env = (sigma_pht_env*((pht(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_env*pht(x,y,z)*env(x,y,z)*(pht(x,y,z)-env(x,y,z))) - &
                & (2*(pht(x,y,z)+env(x,y,z))*(pht(x,y,z)+env(x,y,z))*(w_pht-w_env))

           dF_dpht_pyr = (sigma_pht_pyr*((pht(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_pyr*pht(x,y,z)*pyr(x,y,z)*(pht(x,y,z)-pyr(x,y,z))) - &
                & (2*(pht(x,y,z)+pyr(x,y,z))*(pht(x,y,z)+pyr(x,y,z))*(w_pht-w_pyr))


           dF_denv_met = (sigma_env_met*((env(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_met*env(x,y,z)*met(x,y,z)*(env(x,y,z)-met(x,y,z))) - &
                & (2*(env(x,y,z)+met(x,y,z))*(env(x,y,z)+met(x,y,z))*(w_env-w_met))

           dF_denv_pht = (sigma_env_pht*((env(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_pht*env(x,y,z)*pht(x,y,z)*(env(x,y,z)-pht(x,y,z))) - &
                & (2*(env(x,y,z)+pht(x,y,z))*(env(x,y,z)+pht(x,y,z))*(w_env-w_pht))

           dF_denv_pyr = (sigma_env_pyr*((env(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_pyr*env(x,y,z)*pyr(x,y,z)*(env(x,y,z)-pyr(x,y,z))) - &
                & (2*(env(x,y,z)+pyr(x,y,z))*(env(x,y,z)+pyr(x,y,z))*(w_env-w_pyr))


           dF_dmet_pht = (sigma_met_pht*((met(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_pht*met(x,y,z)*pht(x,y,z)*(met(x,y,z)-pht(x,y,z))) - &
                & (2*(met(x,y,z)+pht(x,y,z))*(met(x,y,z)+pht(x,y,z))*(w_met-w_pht))

           dF_dmet_env = (sigma_met_env*((met(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_env*met(x,y,z)*env(x,y,z)*(met(x,y,z)-env(x,y,z))) - &
                & (2*(met(x,y,z)+env(x,y,z))*(met(x,y,z)+env(x,y,z))*(w_met-w_env))

           dF_dmet_pyr = (sigma_met_pyr*((met(x,y,z)*del2pyr(x,y,z))-(pyr(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_pyr*met(x,y,z)*pyr(x,y,z)*(met(x,y,z)-pyr(x,y,z))) - &
                & (2*(met(x,y,z)+pyr(x,y,z))*(met(x,y,z)+pyr(x,y,z))*(w_met-w_pyr))


           dF_dpyr_met = (sigma_pyr_met*((pyr(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_met*pyr(x,y,z)*met(x,y,z)*(pyr(x,y,z)-met(x,y,z))) - &
                & (2*(pyr(x,y,z)+met(x,y,z))*(pyr(x,y,z)+met(x,y,z))*(w_pyr-w_met))

           dF_dpyr_pht = (sigma_pyr_pht*((pyr(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_pht*pyr(x,y,z)*pht(x,y,z)*(pyr(x,y,z)-pht(x,y,z))) - &
                & (2*(pyr(x,y,z)+pht(x,y,z))*(pyr(x,y,z)+pht(x,y,z))*(w_pyr-w_pht))

           dF_dpyr_env = (sigma_pyr_env*((pyr(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2pyr(x,y,z)))) + &
                & (2*hill_pyr_env*pyr(x,y,z)*env(x,y,z)*(pyr(x,y,z)-env(x,y,z))) - &
                & (2*(pyr(x,y,z)+env(x,y,z))*(pyr(x,y,z)+env(x,y,z))*(w_pyr-w_env))



           !! Error correction        
           correct_rounding = 0.50d0*(dF_dpht_met - dF_dmet_pht)
           dF_dpht_met = correct_rounding ; dF_dmet_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpht_env - dF_denv_pht)
           dF_dpht_env = correct_rounding ; dF_denv_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpht_pyr - dF_dpyr_pht)
           dF_dpht_pyr = correct_rounding ; dF_dpyr_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_denv_met - dF_dmet_env)
           dF_denv_met = correct_rounding ; dF_dmet_env = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_denv_pyr - dF_dpyr_env)
           dF_denv_pyr = correct_rounding ; dF_dpyr_env = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dmet_pyr - dF_dpyr_met)
           dF_dmet_pyr = correct_rounding ; dF_dpyr_met = 0.0d0 - correct_rounding

           !! Calculate phase field(s) update
           dpht_dt(x,y,z) = (M_pht_met*dF_dpht_met)+(M_pht_env*dF_dpht_env)+(M_pht_pyr*dF_dpht_pyr)
           denv_dt(x,y,z) = (M_env_met*dF_denv_met)+(M_env_pht*dF_denv_pht)+(M_env_pyr*dF_denv_pyr)
           dmet_dt(x,y,z) = (M_met_pht*dF_dmet_pht)+(M_met_env*dF_dmet_env)+(M_met_pyr*dF_dmet_pyr)
           dpyr_dt(x,y,z) = (M_pyr_pht*dF_dpyr_pht)+(M_pyr_env*dF_dpyr_env)+(M_pyr_met*dF_dpyr_met)
        end do
     end do
  end do



  !! Apply boundary conditions to phase field(s) update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           dpht_dt(x,y,2) = 0.0d0
           denv_dt(x,y,2) = 0.0d0
           dmet_dt(x,y,2) = 0.0d0
           dpyr_dt(x,y,2) = 0.0d0
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           dpht_dt(x,y,psz+1) = 0.0d0
           denv_dt(x,y,psz+1) = 0.0d0
           dmet_dt(x,y,psz+1) = 0.0d0
           dpyr_dt(x,y,psz+1) = 0.0d0
        end do
     end do
  end if


  !! Update phase fields
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmet(x,y,z) = max(min((met(x,y,z) + (dt*dmet_dt(x,y,z))),1.0d0),0.0d0)
           newpht(x,y,z) = max(min((pht(x,y,z) + (dt*dpht_dt(x,y,z))),1.0d0),0.0d0)
           newenv(x,y,z) = max(min((env(x,y,z) + (dt*denv_dt(x,y,z))),1.0d0),0.0d0)
           newpyr(x,y,z) = max(min((pyr(x,y,z) + (dt*dpyr_dt(x,y,z))),1.0d0),0.0d0)

           newmet(x,y,z) = newmet(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newpht(x,y,z) = newpht(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newenv(x,y,z) = newenv(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newpyr(x,y,z) = newpyr(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
        end do
     end do
  end do


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (env(x,y,z).gt.0.99) then
              newenv(x,y,z) = 1.0d0
              newmet(x,y,z) = 0.0d0
              newpht(x,y,z) = 0.0d0
              newpyr(x,y,z) = 0.0d0
           end if
        end do
     end do
  end do


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           dpht_dt(x,y,z) = (newpht(x,y,z) - pht(x,y,z))/dt
           denv_dt(x,y,z) = (newenv(x,y,z) - env(x,y,z))/dt
           dmet_dt(x,y,z) = (newmet(x,y,z) - met(x,y,z))/dt
           dpyr_dt(x,y,z) = (newpyr(x,y,z) - pyr(x,y,z))/dt
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
              sulfidation_rate = sulfidation_rate + 0.04356E-9*(exp(avg_mu_env/(R*T)) - exp(-2.3025/(R*T)))  !! Ref = Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
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
           pht(x,y,z) = newpht(x,y,z)
           env(x,y,z) = newenv(x,y,z)
           met(x,y,z) = newmet(x,y,z)
           pyr(x,y,z) = newpyr(x,y,z)
        end do
     end do
  end do

end subroutine pfsolve
