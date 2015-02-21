subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  external F_Met_Pht, F_Met_Pyr, F_Met_Env
  external F_Pht_Met, F_Pht_Pyr, F_Pht_Env
  external F_Pyr_Met, F_Pyr_Pht, F_Pyr_Env
  external F_Env_Met, F_Env_Pht, F_Env_Pyr
  external J_Met_Pht, J_Met_Pyr, J_Met_Env
  external J_Pht_Met, J_Pht_Pyr, J_Pht_Env
  external J_Pyr_Met, J_Pyr_Pht, J_Pyr_Env
  external J_Env_Met, J_Env_Pht, J_Env_Pyr

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  PetscErrorCode ierr
  Vec met_vec, pht_vec, pyr_vec, env_vec
  Vec rhs_met_vec, rhs_pht_vec, rhs_pyr_vec, rhs_env_vec
  Vec ret_met_vec, ret_pht_vec, ret_pyr_vec, ret_env_vec
  Mat jac_met, jac_pht, jac_pyr, jac_env
  SNES snes_met, snes_pht, snes_pyr, snes_env
  KSP ksp_met, ksp_pht, ksp_pyr, ksp_env
  PC pc_met, pc_pht, pc_pyr, pc_env
  KSPType ksptype
  PCType pctype
  PetscScalar, pointer :: point_met_vec(:)
  PetscScalar, pointer :: point_pht_vec(:)
  PetscScalar, pointer :: point_pyr_vec(:)
  PetscScalar, pointer :: point_env_vec(:)

  integer :: x, y, z   ! Loop variables
  integer :: linindex
  integer :: int_count
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter

  real*8 :: D_local
  integer :: wrap
  real*8 :: delx,dely,delz

  !! Bulk free energy 
  real*8 :: f_met, f_pyr, f_pht, f_env
  real*8 :: w_met, w_pyr, w_pht, w_env

  ! Vectors for implicit solver
  real*8, dimension(psx*psy*psz) :: B, approxsol
  real*8, dimension(psx,psy,psz+2) :: dmet_pht, dmet_pyr, dmet_env
  real*8, dimension(psx,psy,psz+2) :: dpht_met, dpht_pyr, dpht_env
  real*8, dimension(psx,psy,psz+2) :: dpyr_pht, dpyr_met, dpyr_env
  real*8, dimension(psx,psy,psz+2) :: denv_pht, denv_met, denv_pyr
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8 :: sum_fields

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (mod(iter,swap_freq_pf).eq.1) then
     call swap_pf()
  end if




!--------------------------------------------------------------------!
!-------------------------MET-PHT------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_met_vec,ierr)
  call VecSetSizes(ret_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_met_vec,ierr)
  call VecSetUp(ret_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (met(x,y,z+1)/dt) - (2*M_met_pht*w_met*(pht(x,y,z+1)*pht(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_met_pht*sigma_met_pht*(pht(x,y,1)+pht(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*met(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_met_pht*sigma_met_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)))
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

  call VecAssemblyBegin(met_vec,ierr)
  call VecAssemblyEnd(met_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_met_vec,ierr)
  call VecSetSizes(rhs_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_met_vec,ierr)
  call VecSetUp(rhs_met_vec,ierr)
  call VecSetValues(rhs_met_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_met_vec,ierr)
  call VecAssemblyEnd(rhs_met_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_met,ierr)
  call MatSetSizes(jac_met,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_met,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_met,ierr)
  call SNESSetFunction(snes_met,ret_met_vec,F_Met_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_met,jac_met,jac_met,J_Met_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_met,ierr)
  call SNESGetKSP(snes_met,ksp_met,ierr)
  call KSPSetType(ksp_met,KSPCG,ierr)
  call KSPGetPC(ksp_met,pc_met,ierr)
  call PCSetType(pc_met,PCGAMG,ierr)
  call SNESSolve(snes_met,rhs_met_vec,met_vec,ierr)


  call VecGetArrayF90(met_vec,point_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dmet_pht(x,y,z+1) =  point_met_vec(linindex)-met(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(met_vec,point_met_vec,ierr)

  call VecDestroy(met_vec,ierr)
  call VecDestroy(rhs_met_vec,ierr)
  call MatDestroy(jac_met,ierr)
  call SNESDestroy(snes_met,ierr)


























!--------------------------------------------------------------------!
!-------------------------MET-PYR------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_met_vec,ierr)
  call VecSetSizes(ret_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_met_vec,ierr)
  call VecSetUp(ret_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (met(x,y,z+1)/dt) - (2*M_met_pyr*w_met*(pyr(x,y,z+1)*pyr(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_met_pyr*sigma_met_pyr*(pyr(x,y,1)+pyr(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*met(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_met_pyr*sigma_met_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)))
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

  call VecAssemblyBegin(met_vec,ierr)
  call VecAssemblyEnd(met_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_met_vec,ierr)
  call VecSetSizes(rhs_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_met_vec,ierr)
  call VecSetUp(rhs_met_vec,ierr)
  call VecSetValues(rhs_met_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_met_vec,ierr)
  call VecAssemblyEnd(rhs_met_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_met,ierr)
  call MatSetSizes(jac_met,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_met,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_met,ierr)
  call SNESSetFunction(snes_met,ret_met_vec,F_Met_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_met,jac_met,jac_met,J_Met_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_met,ierr)
  call SNESGetKSP(snes_met,ksp_met,ierr)
  call KSPSetType(ksp_met,KSPCG,ierr)
  call KSPGetPC(ksp_met,pc_met,ierr)
  call PCSetType(pc_met,PCGAMG,ierr)
  call SNESSolve(snes_met,rhs_met_vec,met_vec,ierr)


  call VecGetArrayF90(met_vec,point_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           dmet_pyr(x,y,z+1) =  point_met_vec(linindex)-met(x,y,z+1)
        end do
     end do
  end do

  call VecRestoreArrayF90(met_vec,point_met_vec,ierr)

  call VecDestroy(met_vec,ierr)
  call VecDestroy(rhs_met_vec,ierr)
  call MatDestroy(jac_met,ierr)
  call SNESDestroy(snes_met,ierr)



























!--------------------------------------------------------------------!
!-------------------------MET-ENV------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_met_vec,ierr)
  call VecSetSizes(ret_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_met_vec,ierr)
  call VecSetUp(ret_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (met(x,y,z+1)/dt) - (2*M_met_env*w_met*(env(x,y,z+1)*env(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_met_env*sigma_met_env*(env(x,y,1)+env(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*met(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_met_env*sigma_met_env*(env(x,y,psz+1)+env(x,y,psz+2)))
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

  call VecAssemblyBegin(met_vec,ierr)
  call VecAssemblyEnd(met_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_met_vec,ierr)
  call VecSetSizes(rhs_met_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_met_vec,ierr)
  call VecSetUp(rhs_met_vec,ierr)
  call VecSetValues(rhs_met_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_met_vec,ierr)
  call VecAssemblyEnd(rhs_met_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_met,ierr)
  call MatSetSizes(jac_met,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_met,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_met,ierr)
  call SNESSetFunction(snes_met,ret_met_vec,F_Met_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_met,jac_met,jac_met,J_Met_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_met,ierr)
  call SNESGetKSP(snes_met,ksp_met,ierr)
  call KSPSetType(ksp_met,KSPCG,ierr)
  call KSPGetPC(ksp_met,pc_met,ierr)
  call PCSetType(pc_met,PCGAMG,ierr)
  call SNESSolve(snes_met,rhs_met_vec,met_vec,ierr)


  call VecGetArrayF90(met_vec,point_met_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dmet_env(x,y,z+1) =  point_met_vec(linindex)-met(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(met_vec,point_met_vec,ierr)

  call VecDestroy(met_vec,ierr)
  call VecDestroy(rhs_met_vec,ierr)
  call MatDestroy(jac_met,ierr)
  call SNESDestroy(snes_met,ierr)























!--------------------------------------------------------------------!
!-------------------------PHT-MET------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pht_vec,ierr)
  call VecSetSizes(ret_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pht_vec,ierr)
  call VecSetUp(ret_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pht(x,y,z+1)/dt) - (2*M_pht_met*w_pht*(met(x,y,z+1)*met(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pht_met*sigma_pht_met*(met(x,y,1)+met(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pht(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pht_met*sigma_pht_met*(met(x,y,psz+1)+met(x,y,psz+2)))
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

  call VecAssemblyBegin(pht_vec,ierr)
  call VecAssemblyEnd(pht_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pht_vec,ierr)
  call VecSetSizes(rhs_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pht_vec,ierr)
  call VecSetUp(rhs_pht_vec,ierr)
  call VecSetValues(rhs_pht_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pht_vec,ierr)
  call VecAssemblyEnd(rhs_pht_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pht,ierr)
  call MatSetSizes(jac_pht,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pht,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pht,ierr)
  call SNESSetFunction(snes_pht,ret_pht_vec,F_Pht_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pht,jac_pht,jac_pht,J_Pht_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pht,ierr)
  call SNESGetKSP(snes_pht,ksp_pht,ierr)
  call KSPSetType(ksp_pht,KSPCG,ierr)
  call KSPGetPC(ksp_pht,pc_pht,ierr)
  call PCSetType(pc_pht,PCGAMG,ierr)
  call SNESSolve(snes_pht,rhs_pht_vec,pht_vec,ierr)


  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpht_met(x,y,z+1) =  point_pht_vec(linindex)-pht(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)

  call VecDestroy(pht_vec,ierr)
  call VecDestroy(rhs_pht_vec,ierr)
  call MatDestroy(jac_pht,ierr)
  call SNESDestroy(snes_pht,ierr)
















!--------------------------------------------------------------------!
!-------------------------PHT-PYR------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pht_vec,ierr)
  call VecSetSizes(ret_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pht_vec,ierr)
  call VecSetUp(ret_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pht(x,y,z+1)/dt) - (2*M_pht_pyr*w_pht*(pyr(x,y,z+1)*pyr(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pht_pyr*sigma_pht_pyr*(pyr(x,y,1)+pyr(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pht(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pht_pyr*sigma_pht_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)))
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

  call VecAssemblyBegin(pht_vec,ierr)
  call VecAssemblyEnd(pht_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pht_vec,ierr)
  call VecSetSizes(rhs_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pht_vec,ierr)
  call VecSetUp(rhs_pht_vec,ierr)
  call VecSetValues(rhs_pht_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pht_vec,ierr)
  call VecAssemblyEnd(rhs_pht_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pht,ierr)
  call MatSetSizes(jac_pht,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pht,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pht,ierr)
  call SNESSetFunction(snes_pht,ret_pht_vec,F_Pht_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pht,jac_pht,jac_pht,J_Pht_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pht,ierr)
  call SNESGetKSP(snes_pht,ksp_pht,ierr)
  call KSPSetType(ksp_pht,KSPCG,ierr)
  call KSPGetPC(ksp_pht,pc_pht,ierr)
  call PCSetType(pc_pht,PCGAMG,ierr)
  call SNESSolve(snes_pht,rhs_pht_vec,pht_vec,ierr)


  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpht_pyr(x,y,z+1) =  point_pht_vec(linindex)-pht(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)

  call VecDestroy(pht_vec,ierr)
  call VecDestroy(rhs_pht_vec,ierr)
  call MatDestroy(jac_pht,ierr)
  call SNESDestroy(snes_pht,ierr)




























!--------------------------------------------------------------------!
!-------------------------PHT-ENV------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pht_vec,ierr)
  call VecSetSizes(ret_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pht_vec,ierr)
  call VecSetUp(ret_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pht(x,y,z+1)/dt) - (2*M_pht_env*w_pht*(env(x,y,z+1)*env(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pht_env*sigma_pht_env*(env(x,y,1)+env(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pht(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pht_env*sigma_pht_env*(env(x,y,psz+1)+env(x,y,psz+2)))
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

  call VecAssemblyBegin(pht_vec,ierr)
  call VecAssemblyEnd(pht_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pht_vec,ierr)
  call VecSetSizes(rhs_pht_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pht_vec,ierr)
  call VecSetUp(rhs_pht_vec,ierr)
  call VecSetValues(rhs_pht_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pht_vec,ierr)
  call VecAssemblyEnd(rhs_pht_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pht,ierr)
  call MatSetSizes(jac_pht,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pht,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pht,ierr)
  call SNESSetFunction(snes_pht,ret_pht_vec,F_Pht_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pht,jac_pht,jac_pht,J_Pht_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pht,ierr)
  call SNESGetKSP(snes_pht,ksp_pht,ierr)
  call KSPSetType(ksp_pht,KSPCG,ierr)
  call KSPGetPC(ksp_pht,pc_pht,ierr)
  call PCSetType(pc_pht,PCGAMG,ierr)
  call SNESSolve(snes_pht,rhs_pht_vec,pht_vec,ierr)


  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpht_env(x,y,z+1) =  point_pht_vec(linindex)-pht(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)

  call VecDestroy(pht_vec,ierr)
  call VecDestroy(rhs_pht_vec,ierr)
  call MatDestroy(jac_pht,ierr)
  call SNESDestroy(snes_pht,ierr)















!--------------------------------------------------------------------!
!-------------------------PYR-MET------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pyr_vec,ierr)
  call VecSetSizes(ret_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pyr_vec,ierr)
  call VecSetUp(ret_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pyr(x,y,z+1)/dt) - (2*M_pyr_met*w_pyr*(met(x,y,z+1)*met(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pyr_met*sigma_pyr_met*(met(x,y,1)+met(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pyr(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pyr_met*sigma_pyr_met*(met(x,y,psz+1)+met(x,y,psz+2)))
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

  call VecAssemblyBegin(pyr_vec,ierr)
  call VecAssemblyEnd(pyr_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pyr_vec,ierr)
  call VecSetSizes(rhs_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pyr_vec,ierr)
  call VecSetUp(rhs_pyr_vec,ierr)
  call VecSetValues(rhs_pyr_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pyr_vec,ierr)
  call VecAssemblyEnd(rhs_pyr_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pyr,ierr)
  call MatSetSizes(jac_pyr,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pyr,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pyr,ierr)
  call SNESSetFunction(snes_pyr,ret_pyr_vec,F_Pyr_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pyr,jac_pyr,jac_pyr,J_Pyr_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pyr,ierr)
  call SNESGetKSP(snes_pyr,ksp_pyr,ierr)
  call KSPSetType(ksp_pyr,KSPCG,ierr)
  call KSPGetPC(ksp_pyr,pc_pyr,ierr)
  call PCSetType(pc_pyr,PCGAMG,ierr)
  call SNESSolve(snes_pyr,rhs_pyr_vec,pyr_vec,ierr)


  call VecGetArrayF90(pyr_vec,point_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpyr_met(x,y,z+1) =  point_pyr_vec(linindex)-pyr(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)

  call VecDestroy(pyr_vec,ierr)
  call VecDestroy(rhs_pyr_vec,ierr)
  call MatDestroy(jac_pyr,ierr)
  call SNESDestroy(snes_pyr,ierr)
















!--------------------------------------------------------------------!
!-------------------------PYR-PHT------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pyr_vec,ierr)
  call VecSetSizes(ret_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pyr_vec,ierr)
  call VecSetUp(ret_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pyr(x,y,z+1)/dt) - (2*M_pyr_pht*w_pyr*(pht(x,y,z+1)*pht(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pyr_pht*sigma_pyr_pht*(pht(x,y,1)+pht(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pyr(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pyr_pht*sigma_pyr_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)))
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

  call VecAssemblyBegin(pyr_vec,ierr)
  call VecAssemblyEnd(pyr_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pyr_vec,ierr)
  call VecSetSizes(rhs_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pyr_vec,ierr)
  call VecSetUp(rhs_pyr_vec,ierr)
  call VecSetValues(rhs_pyr_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pyr_vec,ierr)
  call VecAssemblyEnd(rhs_pyr_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pyr,ierr)
  call MatSetSizes(jac_pyr,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pyr,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pyr,ierr)
  call SNESSetFunction(snes_pyr,ret_pyr_vec,F_Pyr_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pyr,jac_pyr,jac_pyr,J_Pyr_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pyr,ierr)
  call SNESGetKSP(snes_pyr,ksp_pyr,ierr)
  call KSPSetType(ksp_pyr,KSPCG,ierr)
  call KSPGetPC(ksp_pyr,pc_pyr,ierr)
  call PCSetType(pc_pyr,PCGAMG,ierr)
  call SNESSolve(snes_pyr,rhs_pyr_vec,pyr_vec,ierr)


  call VecGetArrayF90(pyr_vec,point_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpyr_pht(x,y,z+1) =  point_pyr_vec(linindex)-pyr(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)

  call VecDestroy(pyr_vec,ierr)
  call VecDestroy(rhs_pyr_vec,ierr)
  call MatDestroy(jac_pyr,ierr)
  call SNESDestroy(snes_pyr,ierr)





!--------------------------------------------------------------------!
!-------------------------PYR-ENV------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_pyr_vec,ierr)
  call VecSetSizes(ret_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_pyr_vec,ierr)
  call VecSetUp(ret_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (pyr(x,y,z+1)/dt) - (2*M_pyr_env*w_pyr*(env(x,y,z+1)*env(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_pyr_env*sigma_pyr_env*(env(x,y,1)+env(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*pyr(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_pyr_env*sigma_pyr_env*(env(x,y,psz+1)+env(x,y,psz+2)))
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

  call VecAssemblyBegin(pyr_vec,ierr)
  call VecAssemblyEnd(pyr_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_pyr_vec,ierr)
  call VecSetSizes(rhs_pyr_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_pyr_vec,ierr)
  call VecSetUp(rhs_pyr_vec,ierr)
  call VecSetValues(rhs_pyr_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_pyr_vec,ierr)
  call VecAssemblyEnd(rhs_pyr_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_pyr,ierr)
  call MatSetSizes(jac_pyr,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_pyr,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pyr,ierr)
  call SNESSetFunction(snes_pyr,ret_pyr_vec,F_Pyr_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pyr,jac_pyr,jac_pyr,J_Pyr_Env,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pyr,ierr)
  call SNESGetKSP(snes_pyr,ksp_pyr,ierr)
  call KSPSetType(ksp_pyr,KSPCG,ierr)
  call KSPGetPC(ksp_pyr,pc_pyr,ierr)
  call PCSetType(pc_pyr,PCGAMG,ierr)
  call SNESSolve(snes_pyr,rhs_pyr_vec,pyr_vec,ierr)


  call VecGetArrayF90(pyr_vec,point_pyr_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           dpyr_env(x,y,z+1) =  point_pyr_vec(linindex)-pyr(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)

  call VecDestroy(pyr_vec,ierr)
  call VecDestroy(rhs_pyr_vec,ierr)
  call MatDestroy(jac_pyr,ierr)
  call SNESDestroy(snes_pyr,ierr)





















!--------------------------------------------------------------------!
!-------------------------ENV-MET------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_env_vec,ierr)
  call VecSetSizes(ret_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_env_vec,ierr)
  call VecSetUp(ret_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (env(x,y,z+1)/dt) - (2*M_env_met*w_env*(met(x,y,z+1)*met(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_env_met*sigma_env_met*(met(x,y,1)+met(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*env(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_env_met*sigma_env_met*(met(x,y,psz+1)+met(x,y,psz+2)))
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

  call VecAssemblyBegin(env_vec,ierr)
  call VecAssemblyEnd(env_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_env_vec,ierr)
  call VecSetSizes(rhs_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_env_vec,ierr)
  call VecSetUp(rhs_env_vec,ierr)
  call VecSetValues(rhs_env_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_env_vec,ierr)
  call VecAssemblyEnd(rhs_env_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_env,ierr)
  call MatSetSizes(jac_env,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_env,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_env,ierr)
  call SNESSetFunction(snes_env,ret_env_vec,F_Env_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_env,jac_env,jac_env,J_Env_Met,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_env,ierr)
  call SNESGetKSP(snes_env,ksp_env,ierr)
  call KSPSetType(ksp_env,KSPCG,ierr)
  call KSPGetPC(ksp_env,pc_env,ierr)
  call PCSetType(pc_env,PCGAMG,ierr)
  call SNESSolve(snes_env,rhs_env_vec,env_vec,ierr)


  call VecGetArrayF90(env_vec,point_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           denv_met(x,y,z+1) =  point_env_vec(linindex)-env(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(env_vec,point_env_vec,ierr)

  call VecDestroy(env_vec,ierr)
  call VecDestroy(rhs_env_vec,ierr)
  call MatDestroy(jac_env,ierr)
  call SNESDestroy(snes_env,ierr)











!--------------------------------------------------------------------!
!-------------------------ENV-PHT------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_env_vec,ierr)
  call VecSetSizes(ret_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_env_vec,ierr)
  call VecSetUp(ret_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (env(x,y,z+1)/dt) - (2*M_env_pht*w_env*(pht(x,y,z+1)*pht(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_env_pht*sigma_env_pht*(pht(x,y,1)+pht(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*env(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_env_pht*sigma_env_pht*(pht(x,y,psz+1)+pht(x,y,psz+2)))
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

  call VecAssemblyBegin(env_vec,ierr)
  call VecAssemblyEnd(env_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_env_vec,ierr)
  call VecSetSizes(rhs_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_env_vec,ierr)
  call VecSetUp(rhs_env_vec,ierr)
  call VecSetValues(rhs_env_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_env_vec,ierr)
  call VecAssemblyEnd(rhs_env_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_env,ierr)
  call MatSetSizes(jac_env,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_env,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_env,ierr)
  call SNESSetFunction(snes_env,ret_env_vec,F_Env_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_env,jac_env,jac_env,J_Env_Pht,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_env,ierr)
  call SNESGetKSP(snes_env,ksp_env,ierr)
  call KSPSetType(ksp_env,KSPCG,ierr)
  call KSPGetPC(ksp_env,pc_env,ierr)
  call PCSetType(pc_env,PCGAMG,ierr)
  call SNESSolve(snes_env,rhs_env_vec,env_vec,ierr)


  call VecGetArrayF90(env_vec,point_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           denv_pht(x,y,z+1) =  point_env_vec(linindex)-env(x,y,z+1)

        end do
     end do
  end do

  call VecRestoreArrayF90(env_vec,point_env_vec,ierr)

  call VecDestroy(env_vec,ierr)
  call VecDestroy(rhs_env_vec,ierr)
  call MatDestroy(jac_env,ierr)
  call SNESDestroy(snes_env,ierr)























!--------------------------------------------------------------------!
!-------------------------ENV-PYR------------------------------------!
!--------------------------------------------------------------------!


  call VecCreate(PETSC_COMM_SELF,ret_env_vec,ierr)
  call VecSetSizes(ret_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_env_vec,ierr)
  call VecSetUp(ret_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           include 'define_f_w.user'

           B(linindex) =  (env(x,y,z+1)/dt) - (2*M_env_pyr*w_env*(pyr(x,y,z+1)*pyr(x,y,z+1)))

           if (z.eq.1) then

              D_local = 0.5d0*(M_env_pyr*sigma_env_pyr*(pyr(x,y,1)+pyr(x,y,2)))
              B(linindex) = B(linindex) + ((D_local/(dpf*dpf))*env(x,y,z+1-1))

           elseif (z.eq.psz) then

              D_local = 0.5d0*(M_env_pyr*sigma_env_pyr*(pyr(x,y,psz+1)+pyr(x,y,psz+2)))
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

  call VecAssemblyBegin(env_vec,ierr)
  call VecAssemblyEnd(env_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_env_vec,ierr)
  call VecSetSizes(rhs_env_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_env_vec,ierr)
  call VecSetUp(rhs_env_vec,ierr)
  call VecSetValues(rhs_env_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_env_vec,ierr)
  call VecAssemblyEnd(rhs_env_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac_env,ierr)
  call MatSetSizes(jac_env,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac_env,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_env,ierr)
  call SNESSetFunction(snes_env,ret_env_vec,F_Env_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_env,jac_env,jac_env,J_Env_Pyr,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_env,ierr)
  call SNESGetKSP(snes_env,ksp_env,ierr)
  call KSPSetType(ksp_env,KSPCG,ierr)
  call KSPGetPC(ksp_env,pc_env,ierr)
  call PCSetType(pc_env,PCGAMG,ierr)
  call SNESSolve(snes_env,rhs_env_vec,env_vec,ierr)


  call VecGetArrayF90(env_vec,point_env_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           denv_pyr(x,y,z+1) =  point_env_vec(linindex)-env(x,y,z+1)
        end do
     end do
  end do

  call VecRestoreArrayF90(env_vec,point_env_vec,ierr)

  call VecDestroy(env_vec,ierr)
  call VecDestroy(rhs_env_vec,ierr)
  call MatDestroy(jac_env,ierr)
  call SNESDestroy(snes_env,ierr)






!####################################################################################################################################




  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           dmet_pht(x,y,z+1) = min(max((dmet_pht(x,y,z+1) + dpht_met(x,y,z+1)),-0.0005d0),0.0005d0); dpht_met(x,y,z+1) = 0.0d0 - dmet_pht(x,y,z+1)
           dmet_pyr(x,y,z+1) = min(max((dmet_pyr(x,y,z+1) + dpyr_met(x,y,z+1)),-0.0005d0),0.0005d0); dpyr_met(x,y,z+1) = 0.0d0 - dmet_pyr(x,y,z+1)
           dmet_env(x,y,z+1) = min(max((dmet_env(x,y,z+1) + denv_met(x,y,z+1)),-0.0005d0),0.0005d0); denv_met(x,y,z+1) = 0.0d0 - dmet_env(x,y,z+1)
           dpht_pyr(x,y,z+1) = min(max((dpht_pyr(x,y,z+1) + dpyr_pht(x,y,z+1)),-0.0005d0),0.0005d0); dpyr_pht(x,y,z+1) = 0.0d0 - dpht_pyr(x,y,z+1)
           dpht_env(x,y,z+1) = min(max((dpht_env(x,y,z+1) + denv_pht(x,y,z+1)),-0.0005d0),0.0005d0); denv_pht(x,y,z+1) = 0.0d0 - dpht_env(x,y,z+1)
           dpyr_env(x,y,z+1) = min(max((dpyr_env(x,y,z+1) + denv_pyr(x,y,z+1)),-0.0005d0),0.0005d0); denv_pyr(x,y,z+1) = 0.0d0 - dpyr_env(x,y,z+1)

        end do
     end do
  end do





  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           if (rank.eq.4) then
              write(6,*) dmet_pht(x,y,z+1),dpht_met(x,y,z+1)!,dmet_pyr(x,y,z+1),dpyr_met(x,y,z+1),dmet_env(x,y,z+1),denv_met(x,y,z+1)
           end if
        
           newmet(x,y,z+1) = max(min(met(x,y,z+1) + (dmet_pht(x,y,z+1)-dpht_met(x,y,z+1)) + (dmet_pyr(x,y,z+1)-dpyr_met(x,y,z+1)) + (dmet_env(x,y,z+1)-denv_met(x,y,z+1)),1.0d0),0.0d0)
           newpht(x,y,z+1) = max(min(pht(x,y,z+1) + (dpht_met(x,y,z+1)-dmet_pht(x,y,z+1)) + (dpht_pyr(x,y,z+1)-dpyr_pht(x,y,z+1)) + (dpht_env(x,y,z+1)-denv_pht(x,y,z+1)),1.0d0),0.0d0)
           newpyr(x,y,z+1) = max(min(pyr(x,y,z+1) + (dpyr_met(x,y,z+1)-dmet_pyr(x,y,z+1)) + (dpyr_pht(x,y,z+1)-dpht_pyr(x,y,z+1)) + (dpyr_env(x,y,z+1)-denv_pyr(x,y,z+1)),1.0d0),0.0d0)
           newenv(x,y,z+1) = max(min(env(x,y,z+1) + (denv_met(x,y,z+1)-dmet_env(x,y,z+1)) + (denv_pht(x,y,z+1)-dpht_env(x,y,z+1)) + (denv_pyr(x,y,z+1)-dpyr_env(x,y,z+1)),1.0d0),0.0d0)

        end do
     end do
  end do

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




 ! do x = 1,psx
 !     do y = 1,psy
 !        do z = 2,psz+1
 !           sum_fields = newmet(x,y,z)+newpht(x,y,z)+newenv(x,y,z)+newpyr(x,y,z)
 !           newmet(x,y,z) = newmet(x,y,z)/sum_fields
 !           newpht(x,y,z) = newpht(x,y,z)/sum_fields
 !           newenv(x,y,z) = newenv(x,y,z)/sum_fields
 !           newpyr(x,y,z) = newpyr(x,y,z)/sum_fields
 !        end do
 !     end do
 !  end do


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


end subroutine pfsolve
