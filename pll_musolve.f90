subroutine para_musolve(iter,ksp_mu,simstate)
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
  KSP ksp_mu
  KSPConvergedReason mu_converged_reason
  Vec state, solved_mu_vector, state_solved
  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  type(context) simstate
  external computeRHS_mu, computeMatrix_mu, computeInitialGuess_mu

  call KSPSetDM(ksp_mu,simstate%lattval,ierr)
  call KSPSetComputeRHS(ksp_mu,computeRHS_mu,simstate,ierr)
  call KSPSetComputeOperators(ksp_mu,computeMatrix_mu,simstate,ierr)
  call KSPSetComputeInitialGuess(ksp_mu,computeInitialGuess_mu,simstate,ierr)

  call KSPSetFromOptions(ksp_mu,ierr)
  call KSPSolve(ksp_mu,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  call KSPGetSolution(ksp_mu,state_solved,ierr)
  call KSPGetConvergedReason(ksp_mu,mu_converged_reason,ierr)

  if (mu_converged_reason .gt. 0) then
     call DMGetGlobalVector(simstate%lattval,state,ierr)
     call VecCreate(MPI_COMM_WORLD,solved_mu_vector,ierr)
     call VecSetSizes(solved_mu_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(solved_mu_vector,ierr)
     call VecStrideGather(state_solved,nmus,solved_mu_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(solved_mu_vector,nmus,state,INSERT_VALUES,ierr)
     call DMRestoreGlobalVector(simstate%lattval,state,ierr)
     call VecDestroy(solved_mu_vector,ierr)
  else
     write(6,*) 'Chemical potential field evolution did not converge. Reason: ', mu_converged_reason
  end if

end subroutine para_musolve








subroutine computeInitialGuess_mu(ksp_mu,b,simstate,ierr)
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

  KSP ksp_mu
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
end subroutine computeInitialGuess_mu









subroutine computeRHS_mu(ksp_mu,b,simstate,ierr)
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

  KSP ksp_mu
  PetscErrorCode ierr
  Vec state, b, onlymu
  type(context) simstate

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecCreate(MPI_COMM_WORLD,onlymu,ierr)
  call VecSetSizes(onlymu,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(onlymu,ierr)
  call VecStrideGather(state,nmus,onlymu,INSERT_VALUES,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call VecScale(onlymu,(1.0d0/dt),ierr)
  call VecStrideScatter(onlymu,nmus,b,INSERT_VALUES,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  call VecDestroy(onlymu,ierr)
  return
end subroutine computeRHS_mu






subroutine ComputeMatrix_mu(ksp_mu,matoper,matprecond,simstate,ierr)
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

  KSP ksp_mu
  PetscErrorCode ierr
  Vec state,statelocal
  Mat matoper, matprecond
  PetscInt     i,j,k
  PetscScalar  v(7)
  MatStencil   row(4,1),col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  real*8 :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  integer :: nocols
  real*8 :: add_to_v_ij
  type(context) simstate


  D_inter_met = D_S_met
  D_inter_mkw = max(D_Fe_mkw,D_S_mkw)
  D_inter_pht = max(D_Fe_pht,D_S_pht)
  D_inter_pyr = max(D_Fe_pyr,D_S_pyr)
  D_inter_env = D_S_env


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
           row(MatStencil_c,1) = nmus

           nocols = 0
           add_to_v_ij = 0.0d0

           if (k.gt.0) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i,j,k-1)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i,j,k-1)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i,j,k-1)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i,j,k-1)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i,j,k-1))

              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k-1
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (j.gt.0) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i,j-1,k)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i,j-1,k)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i,j-1,k)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i,j-1,k)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i,j-1,k))

              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j-1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (i.gt.0) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i-1,j,k)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i-1,j,k)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i-1,j,k)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i-1,j,k)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i-1,j,k))

              col(MatStencil_i,nocols) = i-1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (i.lt.psx_g-1) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i+1,j,k)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i+1,j,k)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i+1,j,k)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i+1,j,k)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i+1,j,k))

              col(MatStencil_i,nocols) = i+1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (j.lt.psy_g-1) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i,j+1,k)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i,j+1,k)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i,j+1,k)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i,j+1,k)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i,j+1,k))

              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j+1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (k.lt.psz_g-1) then
              nocols = nocols + 1

              v(nocols) = 0.0d0 - 0.5d0*D_inter_met*(statepointer(nmet,i,j,k)+statepointer(nmet,i,j,k+1)) - &
                   & 0.5d0*D_inter_mkw*(statepointer(nmkw,i,j,k)+statepointer(nmkw,i,j,k+1)) - &
                   & 0.5d0*D_inter_pht*(statepointer(npht,i,j,k)+statepointer(npht,i,j,k+1)) - &
                   & 0.5d0*D_inter_pyr*(statepointer(npyr,i,j,k)+statepointer(npyr,i,j,k+1)) - &
                   & 0.5d0*D_inter_env*(statepointer(nenv,i,j,k)+statepointer(nenv,i,j,k+1))

              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k+1
              col(MatStencil_c,nocols) = nmus
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

              nocols = nocols + 1
              v(nocols) = (1.0d0/dt) - add_to_v_ij
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nmus

              call MatSetValuesStencil(matprecond,1,row,nocols,col,v,INSERT_VALUES,ierr)

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  call MatAssemblyBegin(matprecond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FINAL_ASSEMBLY,ierr)

  call VecDestroy(statelocal,ierr)
  return
end subroutine ComputeMatrix_mu

