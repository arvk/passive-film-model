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
  real*8 :: f_pht, f_env, f_met
  real*8 :: w_pht, w_env, w_met

  !! Thermodynamic parameters
  real*8 :: hill_pht_env, hill_pht_met
  real*8 :: hill_met_env, hill_met_pht
  real*8 :: hill_env_pht, hill_env_met

  !! Interfacial energy
  real*8 :: sigma_pht_env, sigma_pht_met
  real*8 :: sigma_met_env, sigma_met_pht
  real*8 :: sigma_env_pht, sigma_env_met

  !! Kinetic parameters (Interface mobility)
  real*8 :: M_pht_met, M_met_pht
  real*8 :: M_pht_env, M_env_pht
  real*8 :: M_met_env, M_env_met

  !! Derivative of bulk free energy with phase field
  real*8 :: dF_dpht_met, dF_dpht_env
  real*8 :: dF_dmet_pht, dF_dmet_env
  real*8 :: dF_denv_met, dF_denv_pht

  !! Timestep for evolution
  real*8 :: dt = 0.000005d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Thermodynamic parameters
  hill_pht_met = 2E4*85700.0d0 ; hill_met_pht = 2E4*85700.0d0
  hill_pht_env = 7E4*26370.0d0 ; hill_env_pht = 7E4*26370.0d0
  hill_met_env = 1E0*95000.0d0 ; hill_env_met = 1E0*95000.0d0

  sigma_pht_env= -7E1*dpf ; sigma_pht_met= -2E2*dpf
  sigma_met_env= -9E1*dpf ; sigma_met_pht= -2E2*dpf
  sigma_env_pht= -7E1*dpf ; sigma_env_met= -9E1*dpf


  !! Kinetic parameters  
  M_pht_met = 1.7E-09 ; M_met_pht = 1.7E-09
  M_pht_env = 1.7E-10 ; M_env_pht = 1.7E-10
  M_met_env = 1.7E-26 ; M_env_met = 1.7E-26

  !! Initialize field variables to 0
  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dmu_dt = 0.0d0 


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

           if (exp(mu(x,y,z)/(R*T)) .lt. 0.0015d0) then
              w_met = f_met - (mu(x,y,z)*exp(mu(x,y,z)/(R*T))*140401)
           else
              w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           end if


           dF_dpht_met = (sigma_pht_met*((pht(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_met*pht(x,y,z)*met(x,y,z)*(pht(x,y,z)-met(x,y,z))) - &
                & (2*(pht(x,y,z)+met(x,y,z))*(pht(x,y,z)+met(x,y,z))*(w_pht-w_met))

           dF_dpht_env = (sigma_pht_env*((pht(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2pht(x,y,z)))) + &
                & (2*hill_pht_env*pht(x,y,z)*env(x,y,z)*(pht(x,y,z)-env(x,y,z))) - &
                & (2*(pht(x,y,z)+env(x,y,z))*(pht(x,y,z)+env(x,y,z))*(w_pht-w_env))


           dF_denv_met = (sigma_env_met*((env(x,y,z)*del2met(x,y,z))-(met(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_met*env(x,y,z)*met(x,y,z)*(env(x,y,z)-met(x,y,z))) - &
                & (2*(env(x,y,z)+met(x,y,z))*(env(x,y,z)+met(x,y,z))*(w_env-w_met))

           dF_denv_pht = (sigma_env_pht*((env(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2env(x,y,z)))) + &
                & (2*hill_env_pht*env(x,y,z)*pht(x,y,z)*(env(x,y,z)-pht(x,y,z))) - &
                & (2*(env(x,y,z)+pht(x,y,z))*(env(x,y,z)+pht(x,y,z))*(w_env-w_pht))


           dF_dmet_pht = (sigma_met_pht*((met(x,y,z)*del2pht(x,y,z))-(pht(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_pht*met(x,y,z)*pht(x,y,z)*(met(x,y,z)-pht(x,y,z))) - &
                & (2*(met(x,y,z)+pht(x,y,z))*(met(x,y,z)+pht(x,y,z))*(w_met-w_pht))

           dF_dmet_env = (sigma_met_env*((met(x,y,z)*del2env(x,y,z))-(env(x,y,z)*del2met(x,y,z)))) + &
                & (2*hill_met_env*met(x,y,z)*env(x,y,z)*(met(x,y,z)-env(x,y,z))) - &
                & (2*(met(x,y,z)+env(x,y,z))*(met(x,y,z)+env(x,y,z))*(w_met-w_env))



           !! Error correction        
           correct_rounding = 0.50d0*(dF_dpht_met - dF_dmet_pht)
           dF_dpht_met = correct_rounding ; dF_dmet_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_dpht_env - dF_denv_pht)
           dF_dpht_env = correct_rounding ; dF_denv_pht = 0.0d0 - correct_rounding

           correct_rounding = 0.50d0*(dF_denv_met - dF_dmet_env)
           dF_denv_met = correct_rounding ; dF_dmet_env = 0.0d0 - correct_rounding


           !! Calculate phase field(s) update
           dpht_dt(x,y,z) = (M_pht_met*dF_dpht_met)+(M_pht_env*dF_dpht_env)
           denv_dt(x,y,z) = (M_env_met*dF_denv_met)+(M_env_pht*dF_denv_pht)
           dmet_dt(x,y,z) = (M_met_pht*dF_dmet_pht)+(M_met_env*dF_dmet_env)
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
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           dpht_dt(x,y,psz+1) = 0.0d0
           denv_dt(x,y,psz+1) = 0.0d0
           dmet_dt(x,y,psz+1) = 0.0d0
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

           newmet(x,y,z) = newmet(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z))
           newpht(x,y,z) = newpht(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z))
           newenv(x,y,z) = newenv(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z))
        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           dpht_dt(x,y,z) = (newpht(x,y,z) - pht(x,y,z))/dt
           denv_dt(x,y,z) = (newenv(x,y,z) - env(x,y,z))/dt
           dmet_dt(x,y,z) = (newmet(x,y,z) - met(x,y,z))/dt
        end do
     end do
  end do


  !!  sulfidation_rate = 0.01372E-9 for the linear rate
  call musolve(dt, 0.01372E-9)


  !###########################################################
  !##################--UPDATE ALL FIELDS--####################
  !###########################################################

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           pht(x,y,z) = newpht(x,y,z)
           env(x,y,z) = newenv(x,y,z)
           met(x,y,z) = newmet(x,y,z)
        end do
     end do
  end do




end subroutine pfsolve







subroutine musolve(dt, sulfidation_rate)
  use commondata
  use laplacians
  use thermo_constants
  use field_evolution
  implicit none

  integer :: x, y, z   ! Loop variables
  real*8 :: asd ! Anti-Surface-Diffusion

  !! Sulfur density in different phases
  real*8 :: rho_pht, rho_env, rho_met

  !! Diffusivities
  real*8 :: D_S_met, D_S_env, D_S_pht
  real*8 :: D_Fe_met, D_Fe_env, D_Fe_pht
  real*8 :: D_inter_met, D_inter_env, D_inter_pht, D
  real*8 :: D_H_met, D_H_env, D_H_pht, D_H

  !! Derivative of sulfur density with chemical potential
  real*8 :: drho_dmu_pht, drho_dmu_env, drho_dmu_met, Chi

  !! Timestep for evolution
  real*8, intent(in) :: dt
  real*8, intent(in) :: sulfidation_rate     ! Sulfidation rate / Film growth rate in m/s

  integer, dimension(psx,psy) :: interface_loc


  real*8, dimension(psx,psy,psz+2) :: newmu

  newmu = 0.0d0

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Calculate sulfur density in each phase
           rho_pht = 52275*(1+(mu(x,y,z)/(2*250896)))
           rho_met = max(exp(mu(x,y,z)/(R*T)),0.0015d0)*140401
           rho_env = exp(mu(x,y,z)/(R*T))*(13303/T)

           !! Calculate inter-diffusivities in each phase 
           D_S_pht = 0.1818E-20 ; D_Fe_pht = 7.2727E-20
           if (abs(mu(x,y,z)/(2*250896)) .lt. 0.05) then
              D_inter_pht = (D_S_pht*(1+(mu(x,y,z)/(2*250896)))) + (D_Fe_pht*(1-(mu(x,y,z)/(2*250896))))
           else
              D_inter_pht = D_Fe_pht
           end if

           D_S_met = 0.1818E-30 ; D_Fe_met = 0.1818E-29
           D_inter_met = 0.20E-18

           D_S_env = 909.0E-20 ; D_Fe_env = 0.0E-20
           D_inter_env = 909.0E-20

           D = 1000.0d0*((pht(x,y,z)*D_inter_pht)+(env(x,y,z)*D_inter_env)+(met(x,y,z)*D_inter_met))

           D_H_met = 0.1818E-30 ; D_H_pht = 0.1818E-29 ; D_H_env = 900.0E-20
           D_H = (met(x,y,z)*D_H_met)+(pht(x,y,z)*D_H_pht)+(env(x,y,z)*D_H_env)

           !! Calculate derivative of sulfur concentration with chemical potential
           drho_dmu_pht = 52275.0d0/(2*250896.0d0)

           if (exp(mu(x,y,z)/(R*T)) .lt. 0.0015) then
              drho_dmu_met = (exp(mu(x,y,z)/(R*T))*140401)/(R*T)
           else
              drho_dmu_met = (0.0015d0*140401)/(R*T)
           end if

           drho_dmu_env = (exp(mu(x,y,z)/(R*T))*(13303/T))/(R*T)

           !! Calculate chemical 'specific heat'
           Chi = (pht(x,y,z)*drho_dmu_pht) + (met(x,y,z)*drho_dmu_met) + (env(x,y,z)*drho_dmu_env)

           !! Apply chemical-potential-field evolution equation
           dmu_dt(x,y,z) = (D*del2mu(x,y,z)) - &
                &       (((dpht_dt(x,y,z)*rho_pht) + (dmet_dt(x,y,z)*rho_met) + (denv_dt(x,y,z)*rho_env))/Chi)
           dph_dt(x,y,z) = (D_H*del2ph(x,y,z))

        end do
     end do
  end do



  !! Apply boundary conditions to chemical-potential-field update
  ! if (rank.eq.0) then
  !    do x = 1,psx
  !       do y = 1,psy
  !          dmu_dt(x,y,2) = 0.0d0
  !          dph_dt(x,y,2) = 0.0d0
  !       end do
  !    end do
  ! elseif(rank.eq.procs-1) then
  !    do x = 1,psx
  !       do y = 1,psy
  !          dmu_dt(x,y,psz+1) = 0.0d0
  !          dph_dt(x,y,psz+1) = 0.0d0
  !       end do
  !    end do
  ! end if


  !!! Update chemical potential field
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmu(x,y,z) = mu(x,y,z) + (dt*dmu_dt(x,y,z))
           newmu(x,y,z) = max(min(newmu(x,y,z),max_mu),min_mu)

           newph(x,y,z) = ph(x,y,z) + (dt*dph_dt(x,y,z))
           newph(x,y,z) = max(min(newph(x,y,z),14.0d0),0.0d0)
        end do
     end do
  end do


  !! Normalize chemical potential in bulk phases
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (env(x,y,z).gt.0.97) then
              newmu(x,y,z) = avg_mu_env
           end if
        end do
     end do
  end do


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (env(x,y,z).gt.0.97) then
              newph(x,y,z) = pH_in
           end if
        end do
     end do
  end do


  !! Impose boundary counditions on the composition field (sulfidation rate)

  interface_loc = 0

  rho_pht = 52275.0d0
  rho_met = 0.0015d0*140401       

  do x = 1,psx
     do y = 1,psy
        do z = psz+1,2,-1
           if ((pht(x,y,z) .gt. 1.0E-3).and.(pht(x,y,z+1) .lt. 1.0E-3)) then
              interface_loc(x,y) = z
              newmu(x,y,interface_loc(x,y)) = newmu(x,y,interface_loc(x,y)) - (0.5d0*dt*dmu_dt(x,y,interface_loc(x,y)))
              newmu(x,y,interface_loc(x,y)) = newmu(x,y,interface_loc(x,y)) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt)/(dpf))
              newmu(x,y,interface_loc(x,y)) = min(newmu(x,y,interface_loc(x,y)),max_mu)
              newmu(x,y,interface_loc(x,y)+1) = newmu(x,y,interface_loc(x,y)-1)
              exit
           end if
        end do
     end do
  end do

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           mu(x,y,z) = newmu(x,y,z)
        end do
     end do
  end do


end subroutine musolve
