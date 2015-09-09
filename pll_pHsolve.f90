subroutine para_pHsolve(iter,ksp_pH,simstate)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  PetscErrorCode ierr
  KSP ksp_pH
  KSPConvergedReason pH_converged_reason
  Vec state, solved_pH_vector, state_solved
  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  type(context) simstate
  external computeRHS_pH, computeMatrix_pH, computeInitialGuess_pH

  call KSPSetDM(ksp_pH,simstate%lattval,ierr)
  call KSPSetComputeRHS(ksp_pH,computeRHS_pH,simstate,ierr)
  call KSPSetComputeOperators(ksp_pH,computeMatrix_pH,simstate,ierr)
  call KSPSetComputeInitialGuess(ksp_pH,computeInitialGuess_pH,simstate,ierr)

  call KSPSetFromOptions(ksp_pH,ierr)
  call KSPSolve(ksp_pH,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  call KSPGetSolution(ksp_pH,state_solved,ierr)
  call KSPGetConvergedReason(ksp_pH,pH_converged_reason,ierr)

  if (pH_converged_reason .gt. 0) then
     call DMGetGlobalVector(simstate%lattval,state,ierr)
     call VecCreate(MPI_COMM_WORLD,solved_pH_vector,ierr)
     call VecSetSizes(solved_pH_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(solved_pH_vector,ierr)
     call VecStrideGather(state_solved,npH,solved_pH_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(solved_pH_vector,npH,state,INSERT_VALUES,ierr)
     call DMRestoreGlobalVector(simstate%lattval,state,ierr)
     call VecDestroy(solved_pH_vector,ierr)
  else
     write(6,*) 'pH field evolution did not converge. Reason: ', pH_converged_reason
  end if

end subroutine para_pHsolve














subroutine computeInitialGuess_pH(ksp_pH,b,simstate,ierr)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_pH
  PetscErrorCode ierr
  Vec state, b
  type(context) simstate

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecDuplicate(state,b,ierr)
  call VecCopy(state,b,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)
  return
end subroutine computeInitialGuess_pH
















subroutine computeRHS_pH(ksp_pH,b,simstate,ierr)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_pH
  PetscErrorCode ierr
  Vec state, onlypH, b
  integer :: i,j,k
  type(context) simstate

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecCreate(MPI_COMM_WORLD,onlypH,ierr)
  call VecSetSizes(onlypH,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(onlypH,ierr)
  call VecStrideGather(state,npH,onlypH,INSERT_VALUES,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call VecScale(onlypH,(1.0d0/dt),ierr)
  call VecStrideScatter(onlypH,npH,b,INSERT_VALUES,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call VecDestroy(onlypH,ierr)
  return
end subroutine computeRHS_pH






subroutine ComputeMatrix_pH(ksp_pH,matoper,matprecond,simstate,ierr)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_pH
  PetscErrorCode ierr
  Vec state,statelocal
  Mat matoper, matprecond
  PetscInt     i,j,k
  PetscScalar  v(7)
  MatStencil   row(4,1),col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  real*8 :: D
  integer :: nocols
  real*8 :: add_to_v_ij
  type(context) simstate

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call DMCreateLocalVector(simstate%lattval,statelocal,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,state,INSERT_VALUES,statelocal,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,state,INSERT_VALUES,statelocal,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  do k=simstate%startz,simstate%startz+simstate%widthz-1
     do j=simstate%starty,simstate%starty+simstate%widthy-1
        do i=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = i
           row(MatStencil_j,1) = j
           row(MatStencil_k,1) = k
           row(MatStencil_c,1) = npH


           if (k.eq.psz_g-1) then

              v(1) = 1.0d0
              col(MatStencil_i,1) = i
              col(MatStencil_j,1) = j
              col(MatStencil_k,1) = k
              col(MatStencil_c,1) = npH

              call MatSetValuesStencil(matprecond,1,row,1,col,v,INSERT_VALUES,ierr)

           else

              D = D_H_env*statepointer(nenv,i,j,k)

              nocols = 0
              add_to_v_ij = 0.0d0

              if (k.gt.0) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - (0.5d0*D)
                 col(MatStencil_i,nocols) = i
                 col(MatStencil_j,nocols) = j
                 col(MatStencil_k,nocols) = k-1
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              if (j.gt.0) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - 0.5d0*D
                 col(MatStencil_i,nocols) = i
                 col(MatStencil_j,nocols) = j-1
                 col(MatStencil_k,nocols) = k
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              if (i.gt.0) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - 0.5d0*D
                 col(MatStencil_i,nocols) = i-1
                 col(MatStencil_j,nocols) = j
                 col(MatStencil_k,nocols) = k
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              if (i.lt.psx_g-1) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - 0.5d0*D
                 col(MatStencil_i,nocols) = i+1
                 col(MatStencil_j,nocols) = j
                 col(MatStencil_k,nocols) = k
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              if (j.lt.psy_g-1) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - 0.5d0*D
                 col(MatStencil_i,nocols) = i
                 col(MatStencil_j,nocols) = j+1
                 col(MatStencil_k,nocols) = k
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              if (k.lt.psz_g-1) then
                 nocols = nocols + 1
                 v(nocols) = 0.0d0 - 0.5d0*D
                 col(MatStencil_i,nocols) = i
                 col(MatStencil_j,nocols) = j
                 col(MatStencil_k,nocols) = k+1
                 col(MatStencil_c,nocols) = npH
                 v(nocols) = v(nocols)/(dpf*dpf)
                 add_to_v_ij = add_to_v_ij + v(nocols)
              end if

              nocols = nocols + 1
              v(nocols) = (1.0d0/dt) - add_to_v_ij
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = npH

              call MatSetValuesStencil(matprecond,1,row,nocols,col,v,INSERT_VALUES,ierr)

           end if

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  call MatAssemblyBegin(matprecond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FINAL_ASSEMBLY,ierr)

  call VecDestroy(statelocal,ierr)
  return
end subroutine ComputeMatrix_pH

