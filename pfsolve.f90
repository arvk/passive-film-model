subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  external MetFunction, MetJacobian, PhtFunction, PhtJacobian, PyrFunction, PyrJacobian, EnvFunction, EnvJacobian

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  PetscErrorCode ierr

  Vec met_vec,rhs_met_vec,ret_met_vec
  Vec pht_vec,rhs_pht_vec,ret_pht_vec
  Vec pyr_vec,rhs_pyr_vec,ret_pyr_vec
  Vec env_vec,rhs_env_vec,ret_env_vec

  Mat jac_met, jac_pht, jac_pyr, jac_env

  SNES snes_met, snes_pht, snes_pyr, snes_env

  PetscScalar, pointer :: point_met_vec(:)
  PetscScalar, pointer :: point_pht_vec(:)
  PetscScalar, pointer :: point_pyr_vec(:)
  PetscScalar, pointer :: point_env_vec(:)


  integer :: x, y, z   ! Loop variables
  integer :: linindex
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter
  integer :: its

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
  real*8, dimension(psx*psy*psz) :: B
  real*8, dimension(psx*psy*psz) :: approxsol
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



  call VecCreate(PETSC_COMM_SELF,ret_met_vec,ierr)
  call VecSetSizes(ret_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_met_vec,ierr)
  call VecSetUp(ret_met_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,ret_pht_vec,ierr)
  call VecSetSizes(ret_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pht_vec,ierr)
  call VecSetUp(ret_pht_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,ret_pyr_vec,ierr)
  call VecSetSizes(ret_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pyr_vec,ierr)
  call VecSetUp(ret_pyr_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,ret_env_vec,ierr)
  call VecSetSizes(ret_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_env_vec,ierr)
  call VecSetUp(ret_env_vec,ierr)


  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (met(x,y,z+1)/dt) - &
                & 2*M_met_pht*(w_met-w_pht)*(pht(x,y,z+1)*pht(x,y,z+1)) - &
                & 2*M_met_pyr*(w_met-w_pyr)*(pyr(x,y,z+1)*pyr(x,y,z+1)) - &
                & 2*M_met_env*(w_met-w_env)*(env(x,y,z+1)*env(x,y,z+1)) 

           if (z.eq.1) then

              D_local = 0.5d0*(M_met_pht*sigma_met_pht*(pht(x,y,1)+pht(x,y,2)) + &
                   & M_met_pyr*sigma_met_pyr*(pyr(x,y,1)+pyr(x,y,2)) + &
                   & M_met_env*sigma_met_env*(env(x,y,1)+env(x,y,2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*met(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_met_pht*sigma_met_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)) + &
                   & M_met_pyr*sigma_met_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)) + &
                   & M_met_env*sigma_met_env*(env(x,y,psz+1)+env(x,y,psz+2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*met(x,y,z+1+1))

           end if

           approxsol(linindex) = met(x,y,z+1)
           vector_locator(linindex) = linindex-1

        end do
     end do
  end do

  call VecCreate(PETSC_COMM_SELF,met_vec,ierr)
  call VecSetSizes(met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(met_vec,ierr)
  call VecSetUp(met_vec,ierr)
  call VecSetValues(met_vec,psx*psy*psz,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_met_vec,ierr)
  call VecSetSizes(rhs_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_met_vec,ierr)
  call VecSetUp(rhs_met_vec,ierr)
  call VecSetValues(rhs_met_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_met,ierr)
  call MatSetSizes(jac_met,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_met,ierr)









  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (pht(x,y,z+1)/dt) - &
                & 2*M_pht_met*(w_pht-w_met)*(met(x,y,z+1)*met(x,y,z+1)) - &
                & 2*M_pht_pyr*(w_pht-w_pyr)*(pyr(x,y,z+1)*pyr(x,y,z+1)) - &
                & 2*M_pht_env*(w_pht-w_env)*(env(x,y,z+1)*env(x,y,z+1)) 

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
  call VecSetValues(pht_vec,psx*psy*psz,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pht_vec,ierr)
  call VecSetSizes(rhs_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pht_vec,ierr)
  call VecSetUp(rhs_pht_vec,ierr)
  call VecSetValues(rhs_pht_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pht,ierr)
  call MatSetSizes(jac_pht,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pht,ierr)























  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (pyr(x,y,z+1)/dt) - &
                & 2*M_pyr_met*(w_pyr-w_met)*(met(x,y,z+1)*met(x,y,z+1)) - &
                & 2*M_pyr_pht*(w_pyr-w_pht)*(pht(x,y,z+1)*pht(x,y,z+1)) - &
                & 2*M_pyr_env*(w_pyr-w_env)*(env(x,y,z+1)*env(x,y,z+1)) 

           if (z.eq.1) then

              D_local = 0.5d0*(M_pyr_met*sigma_pyr_met*(met(x,y,1)+met(x,y,2)) + &
                   & M_pyr_pht*sigma_pyr_pht*(pht(x,y,1)+pht(x,y,2)) + &
                   & M_pyr_env*sigma_pyr_env*(env(x,y,1)+env(x,y,2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pyr(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pyr_met*sigma_pyr_met*(met(x,y,psz+1)+met(x,y,psz+2)) + &
                   & M_pyr_pht*sigma_pyr_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)) + &
                   & M_pyr_env*sigma_pyr_env*(env(x,y,psz+1)+env(x,y,psz+2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pyr(x,y,z+1+1))

           end if

           approxsol(linindex) = pyr(x,y,z+1)
           vector_locator(linindex) = linindex-1

        end do
     end do
  end do

  call VecCreate(PETSC_COMM_SELF,pyr_vec,ierr)
  call VecSetSizes(pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(pyr_vec,ierr)
  call VecSetUp(pyr_vec,ierr)
  call VecSetValues(pyr_vec,psx*psy*psz,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pyr_vec,ierr)
  call VecSetSizes(rhs_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pyr_vec,ierr)
  call VecSetUp(rhs_pyr_vec,ierr)
  call VecSetValues(rhs_pyr_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pyr,ierr)
  call MatSetSizes(jac_pyr,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pyr,ierr)


















  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (env(x,y,z+1)/dt) - &
                & 2*M_env_met*(w_env-w_met)*(met(x,y,z+1)*met(x,y,z+1)) - &
                & 2*M_env_pyr*(w_env-w_pyr)*(pyr(x,y,z+1)*pyr(x,y,z+1)) - &
                & 2*M_env_pht*(w_env-w_pht)*(pht(x,y,z+1)*pht(x,y,z+1)) 

           if (z.eq.1) then

              D_local = 0.5d0*(M_env_met*sigma_env_met*(met(x,y,1)+met(x,y,2)) + &
                   & M_env_pyr*sigma_env_pyr*(pyr(x,y,1)+pyr(x,y,2)) + &
                   & M_env_pht*sigma_env_pht*(pht(x,y,1)+pht(x,y,2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*env(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_env_met*sigma_env_met*(met(x,y,psz+1)+met(x,y,psz+2)) + &
                   & M_env_pyr*sigma_env_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)) + &
                   & M_env_pht*sigma_env_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)))

              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*env(x,y,z+1+1))

           end if

           approxsol(linindex) = env(x,y,z+1)
           vector_locator(linindex) = linindex-1

        end do
     end do
  end do

  call VecCreate(PETSC_COMM_SELF,env_vec,ierr)
  call VecSetSizes(env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(env_vec,ierr)
  call VecSetUp(env_vec,ierr)
  call VecSetValues(env_vec,psx*psy*psz,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_env_vec,ierr)
  call VecSetSizes(rhs_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_env_vec,ierr)
  call VecSetUp(rhs_env_vec,ierr)
  call VecSetValues(rhs_env_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_env,ierr)
  call MatSetSizes(jac_env,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_env,ierr)














  call SNESCreate(PETSC_COMM_SELF,snes_met,ierr)
  call SNESSetFunction(snes_met,ret_met_vec,MetFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_met,jac_met,jac_met,MetJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_met,ierr)
  call SNESSolve(snes_met,rhs_met_vec,met_vec,ierr)


  call SNESCreate(PETSC_COMM_SELF,snes_pht,ierr)
  call SNESSetFunction(snes_pht,ret_pht_vec,PhtFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pht,jac_pht,jac_pht,PhtJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pht,ierr)
  call SNESSolve(snes_pht,rhs_pht_vec,pht_vec,ierr)


  call SNESCreate(PETSC_COMM_SELF,snes_pyr,ierr)
  call SNESSetFunction(snes_pyr,ret_pyr_vec,PyrFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pyr,jac_pyr,jac_pyr,PyrJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pyr,ierr)
  call SNESSolve(snes_pyr,rhs_pyr_vec,pyr_vec,ierr)


  call SNESCreate(PETSC_COMM_SELF,snes_env,ierr)
  call SNESSetFunction(snes_env,ret_env_vec,EnvFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_env,jac_env,jac_env,EnvJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_env,ierr)
  call SNESSolve(snes_env,rhs_env_vec,env_vec,ierr)





  call VecGetArrayF90(met_vec,point_met_vec,ierr)
  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)
  call VecGetArrayF90(pyr_vec,point_pyr_vec,ierr)
  call VecGetArrayF90(env_vec,point_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

              if (point_met_vec(linindex).eq.point_met_vec(linindex)) then
                 newmet(x,y,z+1) =  max(min(point_met_vec(linindex),1.0d0),0.0d0)
              end if

              if (point_pht_vec(linindex).eq.point_pht_vec(linindex)) then
                 newpht(x,y,z+1) =  max(min(point_pht_vec(linindex),1.0d0),0.0d0)
              end if

              if (point_pyr_vec(linindex).eq.point_pyr_vec(linindex)) then
                 newpyr(x,y,z+1) =  max(min(point_pyr_vec(linindex),1.0d0),0.0d0)
              end if

              if (point_env_vec(linindex).eq.point_env_vec(linindex)) then
                 newenv(x,y,z+1) =  max(min(point_env_vec(linindex),1.0d0),0.0d0)
              end if

        end do
     end do
  end do

  call VecRestoreArrayF90(met_vec,point_met_vec,ierr)
  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)
  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)
  call VecRestoreArrayF90(env_vec,point_env_vec,ierr)


  !! Apply boundary conditions to phase field(s) update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           newpht(x,y,2) = pht(x,y,2) 
           newenv(x,y,2) = env(x,y,2)
           newmet(x,y,2) = met(x,y,2)
           newpyr(x,y,2) = pyr(x,y,2)
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           newpht(x,y,psz+1) = pht(x,y,psz+1) 
           newenv(x,y,psz+1) = env(x,y,psz+1)
           newmet(x,y,psz+1) = met(x,y,psz+1)
           newpyr(x,y,psz+1) = pyr(x,y,psz+1)
        end do
     end do
  end if




  !! Update phase fields
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           if (newenv(x,y,z).gt.0.97) then
              newmet(x,y,z) = met(x,y,z)
              newpht(x,y,z) = pht(x,y,z)
              newpyr(x,y,z) = pyr(x,y,z)
              newenv(x,y,z) = env(x,y,z)
           end if
           
        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           newmet(x,y,z) = newmet(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newpht(x,y,z) = newpht(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newenv(x,y,z) = newenv(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
           newpyr(x,y,z) = newpyr(x,y,z)/(newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z))
        end do
     end do
  end do




  ! do x = 1,psx
  !    do y = 1,psy
  !       do z = 2,psz+1
  !          if (env(x,y,z).gt.0.99) then
  !             newenv(x,y,z) = 1.0d0
  !             newmet(x,y,z) = 0.0d0
  !             newpht(x,y,z) = 0.0d0
  !             newpyr(x,y,z) = 0.0d0
  !          end if
  !       end do
  !    end do
  ! end do


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


  call VecDestroy(met_vec,ierr)
  call VecDestroy(rhs_met_vec,ierr)
  call MatDestroy(jac_met,ierr)
  call SNESDestroy(snes_met,ierr)

  call VecDestroy(pht_vec,ierr)
  call VecDestroy(rhs_pht_vec,ierr)
  call MatDestroy(jac_pht,ierr)
  call SNESDestroy(snes_pht,ierr)

  call VecDestroy(pyr_vec,ierr)
  call VecDestroy(rhs_pyr_vec,ierr)
  call MatDestroy(jac_pyr,ierr)
  call SNESDestroy(snes_pyr,ierr)

  call VecDestroy(env_vec,ierr)
  call VecDestroy(rhs_env_vec,ierr)
  call MatDestroy(jac_env,ierr)
  call SNESDestroy(snes_env,ierr)

end subroutine pfsolve
