subroutine pfsolve(iter)
  use commondata
  use laplacians
  use thermo_constants
  use field_evolution
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  real*8 :: sulfidation_rate     ! Sulfidation rate / Film growth rate in m/s
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: ierr,status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter

  !!---------------PF evolution-----------------!!

  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr

  !! Thermodynamic parameters
  real*8 :: hill_pht_env, hill_pht_met, hill_pht_pyr
  real*8 :: hill_met_env, hill_met_pht, hill_met_pyr
  real*8 :: hill_env_pht, hill_env_met, hill_env_pyr
  real*8 :: hill_pyr_env, hill_pyr_met, hill_pyr_pht

  !! Interfacial energy
  real*8 :: sigma_pht_env, sigma_pht_met, sigma_pht_pyr
  real*8 :: sigma_met_env, sigma_met_pht, sigma_met_pyr
  real*8 :: sigma_env_pht, sigma_env_met, sigma_env_pyr
  real*8 :: sigma_pyr_env, sigma_pyr_met, sigma_pyr_pht

  !! Kinetic parameters (Interface mobility)
  real*8 :: M_pht_met, M_met_pht, M_met_pyr
  real*8 :: M_pht_env, M_env_pht, M_env_pyr
  real*8 :: M_met_env, M_env_met, M_pht_pyr
  real*8 :: M_pyr_env, M_pyr_met, M_pyr_pht

  !! Derivative of bulk free energy with phase field
  real*8 :: dF_dpht_met, dF_dpht_env, dF_dpht_pyr
  real*8 :: dF_dmet_pht, dF_dmet_env, dF_dmet_pyr
  real*8 :: dF_denv_met, dF_denv_pht, dF_denv_pyr
  real*8 :: dF_dpyr_met, dF_dpyr_pht, dF_dpyr_env

  !! Timestep for evolution
  real*8 :: dt = 0.75E-3

  integer :: int_count
  integer :: muloop,muloop_max

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Thermodynamic parameters
  hill_pht_met = 8E10 ; hill_met_pht = 8E10 ; hill_met_pyr = 8E10
  hill_pht_env = 5E11 ; hill_env_pht = 5E11 ; hill_env_pyr = 5E7
  hill_met_env = 2E14 ; hill_env_met = 2E14 ; hill_pht_pyr = 5E7
  hill_pyr_env = 5E7 ; hill_pyr_met = 8E10 ; hill_pyr_pht = 5E7

  sigma_pht_env= -2E2*dpf ; sigma_pht_met= -2E2*dpf ; sigma_pht_pyr= -5E0*dpf
  sigma_met_env= -2E3*dpf ; sigma_met_pht= -2E2*dpf ; sigma_met_pyr= -2E2*dpf
  sigma_env_pht= -2E2*dpf ; sigma_env_met= -2E3*dpf ; sigma_env_pyr= -5E0*dpf
  sigma_pyr_pht= -5E0*dpf ; sigma_pyr_met= -2E2*dpf ; sigma_pyr_env= -5E0*dpf

  !! Kinetic parameters  
  M_pht_met = 1.7E-09 ; M_met_pht = 1.7E-09 ; M_met_pyr = 1.7E-09
  M_pht_env = 1.7E-12 ; M_env_pht = 1.7E-12 ; M_env_pyr = 8.7E-09
  M_met_env = 1.7E-26 ; M_env_met = 1.7E-26 ; M_pht_pyr = 8.7E-09
  M_pyr_met = 1.7E-09 ; M_pyr_env = 8.7E-09 ; M_pyr_pht = 8.7E-09

  !! Initialize field variables to 0
  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dpyr_dt = 0.0d0 ; dmu_dt = 0.0d0 

     if (mod(iter,swap_freq_pf).eq.1) then
        call swap_pf()
     end if


  !! Calculate laplacian of phase/composition fields
  call calc_lap()



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


  !!  sulfidation_rate = 0.01372E-9 for the linear rate

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

!  write(6,*) 'Sulf rate is ', sulfidation_rate

muloop_max = 1

do muloop = 1,muloop_max
  call swap_mu()
  call calc_lap()
  call musolve(dt/muloop_max, sulfidation_rate,iter)
end do

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







subroutine musolve(dt, sulfidation_rate,iter)
  use commondata
  use laplacians
  use thermo_constants
  use field_evolution
  use diffusion_constants
  implicit none

  integer, intent(in) :: iter

  integer :: x, y, z   ! Loop variables
  real*8 :: asd ! Anti-Surface-Diffusion

  !! Sulfur density in different phases
  real*8 :: rho_pht, rho_env, rho_met, rho_pyr

  !! Diffusivities
  real*8 :: D_inter_met, D_inter_env, D_inter_pht, D_inter_pyr, D
  real*8 :: D_H_met, D_H_env, D_H_pht, D_H_pyr, D_H

  !! Derivative of sulfur density with chemical potential
  real*8 :: drho_dmu_pht, drho_dmu_env, drho_dmu_met, drho_dmu_pyr, Chi

  !! Timesteps
  real*8, intent(in) :: dt  ! Actual timestep for evolution
  real*8 :: dt_stable       ! Maximum stable forward-euler timestep

  real*8, intent(in) :: sulfidation_rate     ! Sulfidation rate / Film growth rate in m/s

  integer, dimension(psx,psy) :: interface_loc

  real*8, dimension(psx,psy,psz+2) :: newmu





!! Verify that the timestep is stable for forward-Euler integration
  dt_stable = (dpf*dpf)/(2*dmax1(D_S_met, D_S_pht, D_S_pyr,D_Fe_met, D_Fe_pht, D_Fe_pyr))
  write(6,*) "The largest stable timestep for forward-Euler integration of the concentration field is ", dt_stable
  write(6,*) "The timestep currently used is ", dt

  newmu = 0.0d0

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Calculate sulfur density in each phase
           rho_pht = 52275
           rho_met = 0.0015d0*140401

           if (iter.lt.nomc/100) then
              rho_env = rho_met
           else      
              rho_env = 0.0015d0*(13303/T)
           end if

           rho_pyr = 2.0d0*41667.0d0

!              rho_env = 0.0015d0*(13303/T)
        
!           rho_met = exp(avg_mu_met/(R*T))*140401
!           rho_env = exp(avg_mu_env/(R*T))*(13303/T)

           !! Calculate inter-diffusivities in each phase 
!           D_Fe_pht = 10**((0-7056/T)-2.679)    ! Ref = William's data


           ! if (abs(mu(x,y,z)/(2*250896)) .lt. 0.05) then
           !    D_inter_pht = (D_S_pht*(1+(mu(x,y,z)/(2*250896)))) + (D_Fe_pht*(1-(mu(x,y,z)/(2*250896))))
           ! else
           !    D_inter_pht = D_Fe_pht
           ! end if

              D_inter_pht = D_Fe_pht

!           D_Fe_met = 7.87E-7 * exp((96500*(0.0-0.60d0))/(R*T))   ! Ref = PHYSICAL REVIEW B 80, 144111-4 2009
           D_inter_met = D_Fe_met

!           D_S_env = 1.70E-9  ! Ref = J. Chem. Eng. Data 1994, 39, 330-332 and http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
           D_inter_env = D_S_env

!           D_inter_pyr = D_Fe_pyr                               ! S diffusivity should be negligible

           D_inter_pyr = D_inter_pht/10

           D = ((pht(x,y,z)*D_inter_pht)+(env(x,y,z)*D_inter_env)+(met(x,y,z)*D_inter_met)+(pyr(x,y,z)*D_inter_pyr))

           !! Hydrogen diffusivity

!           D_H_env = 9.31E-9                      ! Ref = http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
!           D_H_met = 1.5E-7 * exp(0-(8492/(R*T))) ! Ref = PHYSICAL REVIEW B 70 , 064102-7 (2004)
!           D_H_pht = 1.75E-9 ! Ref = Int. J. Electrochem. Sci.,8, (2013) 2880-2891
!           D_H_pyr = D_H_met ! Assumed. Cannot find reference

!           D_H = (met(x,y,z)*D_H_met)+(pht(x,y,z)*D_H_pht)+(env(x,y,z)*D_H_env)

           !! Calculate derivative of sulfur concentration with chemical potential
           drho_dmu_pht = 52275.0d0/(2*250896.0d0)
           drho_dmu_met = (0.0015d0*140401)/(R*T)

           if (iter.lt.nomc/100) then
              drho_dmu_env = (0.0015d0*140401)/(R*T)
           else
              drho_dmu_env = (0.0015d0*(13303/T))/(R*T)
           end if

!           drho_dmu_pyr = 2*41667.0d0/(2*25089600.0d0)
           drho_dmu_pyr = (2*41667.0d0)/(2*250896.0d0)

           !! Calculate chemical 'specific heat'
           Chi = (pht(x,y,z)*drho_dmu_pht) + (met(x,y,z)*drho_dmu_met) + (env(x,y,z)*drho_dmu_env) + (pyr(x,y,z)*drho_dmu_pyr)

           !! Apply chemical-potential-field evolution equation
           dmu_dt(x,y,z) = (D*del2mu(x,y,z)) - &
                &       (((dpht_dt(x,y,z)*rho_pht) + (dmet_dt(x,y,z)*rho_met) + (denv_dt(x,y,z)*rho_env) + (dpyr_dt(x,y,z)*rho_pyr))/Chi)
!           dph_dt(x,y,z) = (D_H*del2ph(x,y,z))

        end do
     end do
  end do



  !! Apply boundary conditions to chemical-potential-field update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           dmu_dt(x,y,2) = 0.0d0
           dph_dt(x,y,2) = 0.0d0
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           dmu_dt(x,y,psz+1) = 0.0d0
           dph_dt(x,y,psz+1) = 0.0d0
        end do
     end do
  end if


  !!! Update chemical potential field
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmu(x,y,z) = mu(x,y,z) + (dt*dmu_dt(x,y,z))

           ! if (iter.lt.nomc/100) then
           !    newmu(x,y,z) = max(min(newmu(x,y,z),mus_pht_pyr_eqb),mus_met_pht_eqb)
           ! else
           !    newmu(x,y,z) = max(min(newmu(x,y,z),max_mu),min_mu)
           ! end if

              newmu(x,y,z) = max(min(newmu(x,y,z),max_mu),min_mu)
       
           newph(x,y,z) = ph(x,y,z)! + (dt*dph_dt(x,y,z))
           newph(x,y,z) = max(min(newph(x,y,z),14.0d0),0.0d0)
        end do
     end do
  end do


  !! Normalize chemical potential in bulk phases
  ! do x = 1,psx
  !    do y = 1,psy
  !       do z = 2,psz+1
  !          if (env(x,y,z).gt.0.97) then
  !             newmu(x,y,z) = avg_mu_env
  !          end if
  !       end do
  !    end do
  ! end do


  ! do x = 1,psx
  !    do y = 1,psy
  !       do z = 2,psz+1
  !          if (env(x,y,z).gt.0.97) then
  !             newph(x,y,z) = pH_in
  !          end if
  !       end do
  !    end do
  ! end do


  !! Impose boundary counditions on the composition field (sulfidation rate)
  
  if (dt.gt.0) then

     interface_loc = 0

     rho_pht = 52275.0d0
     rho_met = 0.0015d0*140401       

     do x = 1,psx
        do y = 1,psy
           do z = psz+1,2,-1
              if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then
                 interface_loc(x,y) = z
!                 newmu(x,y,interface_loc(x,y)) = newmu(x,y,interface_loc(x,y)) - (0.5d0*dt*dmu_dt(x,y,interface_loc(x,y)))
!                 newmu(x,y,interface_loc(x,y)) = newmu(x,y,interface_loc(x,y)) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt*20)/(dpf))
                 newmu(x,y,interface_loc(x,y)) = mu(x,y,interface_loc(x,y)) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt*20)/(dpf))
                 newmu(x,y,interface_loc(x,y)) = min(newmu(x,y,interface_loc(x,y)),max_mu)
                 newmu(x,y,interface_loc(x,y)+1) = newmu(x,y,interface_loc(x,y)-1)
                 exit
              end if
           end do
        end do
     end do

  end if

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           mu(x,y,z) = newmu(x,y,z)
        end do
     end do
  end do

end subroutine musolve
