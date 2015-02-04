subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  external PhtFunction, PhtJacobian

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  PetscErrorCode ierr
  Vec pht_vec,rhs_pht_vec,ret_pht_vec
  Mat jac_pht
  PetscScalar, pointer :: point_pht_vec(:)
  SNES snes_pht

  integer :: x, y, z   ! Loop variables
  integer :: linindex
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter
  integer :: snesiter

  !!---------------PF evolution-----------------!!

  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr

  !! Derivative of bulk free energy with phase field
  real*8 :: dF_dpht_met, dF_dpht_env, dF_dpht_pyr
  real*8 :: dF_dmet_pht, dF_dmet_env, dF_dmet_pyr
  real*8 :: dF_denv_met, dF_denv_pht, dF_denv_pyr
  real*8 :: dF_dpyr_met, dF_dpyr_pht, dF_dpyr_env



  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  real*8 :: D_local

  integer :: int_count

  real*8, dimension(psx,psy,psz+2) :: del_opyr
  real*8 :: odiff
  integer :: wrap
  real*8 :: delx,dely,delz


  ! Vectors for implicit solver
  real*8, dimension(psx*psy*psz) :: approxsol
  real*8, dimension(psx*psy*psz) :: B
  integer, dimension(psx*psy*psz) :: vector_locator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! Initialize field variables to 0
  dpht_dt = 0.0d0 ; denv_dt = 0.0d0 ; dmet_dt = 0.0d0 ; dpyr_dt = 0.0d0 ; dmu_dt = 0.0d0 

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







  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (pht(x,y,z+1)/dt) + &
                & 2*(w_pht-w_met)*(met(x,y,z+1)*met(x,y,z+1)) + &
                & 2*(w_pht-w_pyr)*(pyr(x,y,z+1)*pyr(x,y,z+1)) + &
                & 2*(w_pht-w_env)*(env(x,y,z+1)*env(x,y,z+1)) 

           if (z.eq.1) then

              D_local = 0.5d0*(M_pht_met*sigma_pht_met*(met(x,y,1)+met(x,y,2)) + &
                   & M_pht_pyr*sigma_pht_pyr*(pyr(x,y,1)+pyr(x,y,2)) + &
                   & M_pht_env*sigma_pht_env*(env(x,y,1)+env(x,y,2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pht(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pht_met*sigma_pht_met*(met(x,y,psz+1)+met(x,y,psz+2)) + &
                   & M_pht_pyr*sigma_pht_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)) + &
                   & M_pht_env*sigma_pht_env*(env(x,y,psz+1)+env(x,y,psz+2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pht(x,y,z+1+1))
           end if

           approxsol(linindex) = pht(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do

  call VecCreate(PETSC_COMM_SELF,pht_vec,ierr)
  call VecSetSizes(pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(pht_vec,ierr)
  call VecSetUp(pht_vec,ierr)


  call VecCreate(PETSC_COMM_SELF,ret_pht_vec,ierr)
  call VecSetSizes(ret_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pht_vec,ierr)
  call VecSetUp(ret_pht_vec,ierr)


  call VecCreate(PETSC_COMM_SELF,rhs_pht_vec,ierr)
  call VecSetSizes(rhs_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pht_vec,ierr)
  call VecSetUp(rhs_pht_vec,ierr)
  call VecSetValues(rhs_pht_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)


  call MatCreate(PETSC_COMM_SELF,jac_pht,ierr)
  call MatSetSizes(jac_pht,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pht,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pht,ierr)
  call SNESSetFunction(snes_pht,ret_pht_vec,PhtFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pht,jac_pht,jac_pht,PhtJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pht,ierr)
  call SNESSolve(snes_pht,rhs_pht_vec,pht_vec,ierr)

  call SNESGetIterationNumber(snes_pht,snesiter,ierr)
  write(6,*) "Rank,iter", rank, snesiter



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


           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env


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
           dpyr_dt(x,y,z) = dpyr_dt(x,y,z) - (25.0d0*pyr(x,y,z)*((M_pyr_pht+M_pyr_pht+M_pyr_pht)*del_opyr(x,y,z)))
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
           pht(x,y,z) = newpht(x,y,z)
           env(x,y,z) = newenv(x,y,z)
           met(x,y,z) = newmet(x,y,z)
           pyr(x,y,z) = newpyr(x,y,z)
        end do
     end do
  end do



  call VecDestroy(pht_vec,ierr)
  call VecDestroy(rhs_pht_vec,ierr)
  call MatDestroy(jac_pht,ierr)
  call SNESDestroy(snes_pht,ierr)


end subroutine pfsolve








subroutine PhtFunction(snes,pht_vec,ret_vec,dummy,ierr)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  SNES snes
  Vec pht_vec,ret_vec
  Mat lhs_pht_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pht_vec(:)

  real*8, dimension(psx,psy,psz+2) :: D_pht

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8, dimension(psx*psy*psz) :: vecread
  real*8, dimension(psx*psy*psz) :: lnr,sqr

  real*8, dimension(:), allocatable :: A, LU
  integer, dimension(:), allocatable :: JA, IA

  integer :: x, y, z   ! Loop variables
  integer :: linindex, contindex
  integer :: wrap

  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap


  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2

           D_pht(x,y,z) = M_pht_met*sigma_pht_met*met(x,y,z) + M_pht_pyr*sigma_pht_pyr*pyr(x,y,z) + M_pht_env*sigma_pht_env*env(x,y,z) 

        end do
     end do
  end do


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

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)





           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do




  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(LU((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0


  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)-1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y-1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x-1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_pht(x,y,(z+1)+1)+D_pht(x,y,(z+1)-1)+&
                &D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,wrap(y-1,psy),z+1)+&
                &D_pht(wrap(x+1,psx),y,z+1)+D_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x+1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)+1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6








  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

        do rowindex = IA(linindex),IA(linindex+1)-2

           is_sorted = .TRUE.

           JAleft = JA(rowindex)
           JAright = JA(rowindex+1)

           if (JAleft .gt. JAright) then

              JAswap = JA(rowindex)
              JA(rowindex) = JA(rowindex+1)
              JA(rowindex+1) = JAswap

              Aswap = A(rowindex)
              A(rowindex) = A(rowindex+1)
              A(rowindex+1) = Aswap

              is_sorted = .FALSE.
              exit
           end if

        end do


     end do

  end do


  IA = IA - 1
  JA = JA - 1


  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_pht_mat,ierr)

  call MatMult(lhs_pht_mat,pht_vec,ret_vec,ierr)


  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

!!! Calculate linear prefactor
           lnr(linindex) = (M_pht_met*sigma_pht_met*del2met(x,y,z+1)) + (2*hill_pht_met*met(x,y,z+1)*met(x,y,z+1)) + (2*(w_pht-w_met)) + &
& (M_pht_pyr*sigma_pht_pyr*del2pyr(x,y,z+1)) + (2*hill_pht_pyr*pyr(x,y,z+1)*pyr(x,y,z+1)) + (2*(w_pht-w_pyr)) + &
& (M_pht_env*sigma_pht_env*del2env(x,y,z+1)) + (2*hill_pht_env*env(x,y,z+1)*env(x,y,z+1)) + (2*(w_pht-w_env)) 

           lnr(linindex) = lnr(linindex)*point_pht_vec(linindex)

!!! Calculate quadratic prefactor
           sqr(linindex) = ((2*(w_pht-w_met)) - 2*hill_pht_met) + ((2*(w_pht-w_pyr)) - 2*hill_pht_pyr) + ((2*(w_pht-w_env)) - 2*hill_pht_env) 

           sqr(linindex) = sqr(linindex)*point_pht_vec(linindex)*point_pht_vec(linindex)
        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)


  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call MatDestroy(lhs_pht_mat,ierr)



  return

end subroutine PhtFunction













subroutine PhtJacobian(snes,pht_vec,pht_jacob,pht_precond,dummy,ierr)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  SNES snes
  Vec pht_vec
  Mat pht_jacob,pht_precond
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pht_vec(:)


  real*8, dimension(psx,psy,psz+2) :: D_pht

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8, dimension(psx*psy*psz) :: vecread
  real*8, dimension(psx*psy*psz) :: lnr,sqr

  real*8, dimension(:), allocatable :: A, LU
  integer, dimension(:), allocatable :: JA, IA

  integer :: x, y, z   ! Loop variables
  integer :: linindex, contindex
  integer :: wrap

  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap


  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env

  PetscInt rowval,colval
  PetscScalar val
  integer :: io,jo


  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2

           D_pht(x,y,z) = M_pht_met*sigma_pht_met*met(x,y,z) + M_pht_pyr*sigma_pht_pyr*pyr(x,y,z) + M_pht_env*sigma_pht_env*env(x,y,z) 

        end do
     end do
  end do


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

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)



           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do


  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(LU((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0

  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

      !!! Calculate linear prefactor
           lnr(linindex) = (M_pht_met*sigma_pht_met*del2met(x,y,z+1)) + (2*hill_pht_met*met(x,y,z+1)*met(x,y,z+1)) + (2*(w_pht-w_met)) + &
& (M_pht_pyr*sigma_pht_pyr*del2pyr(x,y,z+1)) + (2*hill_pht_pyr*pyr(x,y,z+1)*pyr(x,y,z+1)) + (2*(w_pht-w_pyr)) + &
& (M_pht_env*sigma_pht_env*del2env(x,y,z+1)) + (2*hill_pht_env*env(x,y,z+1)*env(x,y,z+1)) + (2*(w_pht-w_env)) 

      !!! Calculate quadratic prefactor
           sqr(linindex) = ((2*(w_pht-w_met)) - 2*hill_pht_met) + ((2*(w_pht-w_pyr)) - 2*hill_pht_pyr) + ((2*(w_pht-w_env)) - 2*hill_pht_env) 

           sqr(linindex) = sqr(linindex)*2.0d0*point_pht_vec(linindex)


           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)-1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y-1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x-1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_pht(x,y,(z+1)+1)+D_pht(x,y,(z+1)-1)+&
                &D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,wrap(y-1,psy),z+1)+&
                &D_pht(wrap(x+1,psx),y,z+1)+D_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           A(contindex) = A(contindex) + lnr(linindex) + sqr(linindex)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x+1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)+1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)


  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6



  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

        do rowindex = IA(linindex),IA(linindex+1)-2

           is_sorted = .TRUE.

           JAleft = JA(rowindex)
           JAright = JA(rowindex+1)

           if (JAleft .gt. JAright) then

              JAswap = JA(rowindex)
              JA(rowindex) = JA(rowindex+1)
              JA(rowindex+1) = JAswap

              Aswap = A(rowindex)
              A(rowindex) = A(rowindex+1)
              A(rowindex+1) = Aswap

              is_sorted = .FALSE.
              exit
           end if

        end do


     end do

  end do


  IA = IA - 1
  JA = JA - 1







  do io = 1,psx*psy*psz ! Nrows

     rowval = io-1

     do jo = 1,IA(io+1)-IA(io)

        colval = JA(IA(io) + jo)

        val = A(IA(io) + jo)

        call MatSetValue(pht_precond,rowval,colval,val,INSERT_VALUES,ierr)

     end do
     

  end do


  call MatAssemblyBegin(pht_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pht_precond,MAT_FINAL_ASSEMBLY,ierr)


  return

end subroutine PhtJacobian
