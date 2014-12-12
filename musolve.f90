subroutine musolve(iter)
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use field_evolution
  use diffusion_constants
  implicit none

  integer, intent(in) :: iter

  integer :: x, y, z   ! Loop variables
  real*8 :: asd ! Anti-Surface-Diffusion

  !! Diffusivities
  real*8 :: D_inter_met, D_inter_env, D_inter_pht, D_inter_pyr, D

  !! Derivative of sulfur density with chemical potential
  real*8 :: drho_dmu_pht, drho_dmu_env, drho_dmu_met, drho_dmu_pyr, Chi

!  real*8, intent(in) :: sulfidation_rate     ! Sulfidation rate / Film growth rate in m/s

  integer, dimension(psx,psy) :: interface_loc

  real*8, dimension(psx,psy,psz+2) :: newmu

  newmu = 0.0d0


  if (mod(iter,swap_freq_pf).eq.1) then
     call swap_mu()
  end if
  call calc_lap_mu()


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Calculate inter-diffusivities in each phase 

           D_inter_pht = D_Fe_pht
           D_inter_met = D_Fe_met
           D_inter_env = D_S_env
           D_inter_pyr = D_inter_pht/10

           D = ((pht(x,y,z)*D_inter_pht)+(env(x,y,z)*D_inter_env)+(met(x,y,z)*D_inter_met)+(pyr(x,y,z)*D_inter_pyr))

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

        end do
     end do
  end do



  !! Apply boundary conditions to chemical-potential-field update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           dmu_dt(x,y,2) = 0.0d0
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           dmu_dt(x,y,psz+1) = 0.0d0
        end do
     end do
  end if


  !!! Update chemical potential field
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmu(x,y,z) = mu(x,y,z) + (dt*dmu_dt(x,y,z))
           newmu(x,y,z) = max(min(newmu(x,y,z),max_mu),min_mu)     
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
