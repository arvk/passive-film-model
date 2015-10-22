subroutine solve_pH(iter,ksp_pH,simstate)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  !! **Set up and solve the linear equations for time-evolution of [H+] in the environment**
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  PetscErrorCode :: ierr
  KSP :: ksp_pH  !! Linear pH field solver
  KSPConvergedReason :: pH_converged_reason
  Vec :: state, solved_pH_vector, state_solved !! Vectors to store function values and solutions
  PetscInt, intent(in) :: iter  !! Current iteration number
  PetscInt :: x, y, z           !! Coordinates inside the simulation system
  type(context) simstate       !! Field variables stored in PETSc vectors and DMDA objects
  external computeRHS_pH, computeMatrix_pH, computeInitialGuess_pH

  call KSPSetDM(ksp_pH,simstate%lattval,ierr)
  call KSPSetComputeRHS(ksp_pH,computeRHS_pH,simstate,ierr)
  call KSPSetComputeOperators(ksp_pH,computeMatrix_pH,simstate,ierr)
  call KSPSetComputeInitialGuess(ksp_pH,computeInitialGuess_pH,simstate,ierr)
  call KSPSetOptionsPrefix(ksp_pH,'pH_',ierr)
  call KSPSetFromOptions(ksp_pH,ierr)
  call KSPSolve(ksp_pH,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  call KSPGetSolution(ksp_pH,state_solved,ierr)
  call KSPGetConvergedReason(ksp_pH,pH_converged_reason,ierr)

  if (pH_converged_reason .gt. 0) then
     call VecCreate(MPI_COMM_WORLD,solved_pH_vector,ierr)
     call VecSetSizes(solved_pH_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(solved_pH_vector,ierr)
     call VecStrideGather(state_solved,npH,solved_pH_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(solved_pH_vector,npH,simstate%slice,INSERT_VALUES,ierr)
     call VecDestroy(solved_pH_vector,ierr)
  else
     write(6,*) 'pH field evolution did not converge. Reason: ', pH_converged_reason
  end if

end subroutine solve_pH














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

  KSP :: ksp_pH
  PetscErrorCode :: ierr
  Vec :: b
  type(context) simstate

  call VecCopy(simstate%slice,b,ierr)

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

  KSP :: ksp_pH
  PetscErrorCode :: ierr
  Vec :: b
  PetscInt :: i, j, k, field, fesphase
  type(context) simstate
  PetscScalar, pointer :: statepointer(:,:,:,:), bpointer(:,:,:,:)
  PetscScalar :: local_Efield, sulfid_bias_Efield

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,b,bpointer,ierr)

  do k=simstate%startz,simstate%startz+simstate%widthz-1
     do j=simstate%starty,simstate%starty+simstate%widthy-1
        do i=simstate%startx,simstate%startx+simstate%widthx-1

           do field = 0, (nfields-1)
              bpointer(field,i,j,k) = 0.0d0
           end do

           bpointer(npH,i,j,k) = statepointer(npH,i,j,k) * (1.0d0/dt)

           if ((statepointer(nenv,i,j,min(k+3,simstate%startz+simstate%widthz-1))-statepointer(nenv,i,j,max(k,simstate%startz))).gt.0.1d0) then

              sulfid_bias_Efield = 1.0d0
              if (include_electro) then
                 local_Efield = (statepointer(npot,i,j,min(k+1,simstate%startz+simstate%widthz-1)) - statepointer(npot,i,j,k))/dpf
                 sulfid_bias_Efield = exp((local_Efield*0.143583*1E-10*96500)/(R*T))
              end if

              do fesphase = nmet,npyr
                 statepointer(npH,i,j,min(k+3,simstate%startz+simstate%widthz-1)) = statepointer(npH,i,j,min(k+3,simstate%startz+simstate%widthz-1)) + (sulfid_bias_Efield*(8.0d0/9.0d0)*((rhoS(max(min(fesphase,npyr),nmkw)))*sulf_rate(fesphase)/(dpf))*statepointer(fesphase,i,j,k))
              end do
           end if

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,b,bpointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)
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

  KSP :: ksp_pH
  PetscErrorCode :: ierr
  Vec :: state, statelocal
  Mat :: matoper, matprecond
  PetscInt :: i, j, k, x, y, z
  PetscScalar :: v(7)
  MatStencil :: row(4,1), col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscScalar :: D
  PetscInt :: nocols
  PetscScalar :: add_to_v_ij
  type(context) simstate
  PetscScalar :: zeromatentry(7)
  PetscInt :: matfield


  call DMCreateLocalVector(simstate%lattval,statelocal,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,simstate%slice,INSERT_VALUES,statelocal,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,simstate%slice,INSERT_VALUES,statelocal,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,statelocal,statepointer,ierr)








  zeromatentry(1) = 0.0d0
  zeromatentry(2) = 0.0d0
  zeromatentry(3) = 0.0d0
  zeromatentry(4) = 0.0d0
  zeromatentry(5) = 0.0d0
  zeromatentry(6) = 0.0d0
  zeromatentry(7) = 0.0d0

  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1
           do matfield = 0,nfields-1   ! Row

              row(MatStencil_i,1) = x
              row(MatStencil_j,1) = y
              row(MatStencil_k,1) = z
              row(MatStencil_c,1) = matfield

              nocols = 0

              if (z.ne.0) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x
                 col(MatStencil_j,nocols) = y
                 col(MatStencil_k,nocols) = z-1
                 col(MatStencil_c,nocols) = matfield
              end if

              if (y.ne.0) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x
                 col(MatStencil_j,nocols) = y-1
                 col(MatStencil_k,nocols) = z
                 col(MatStencil_c,nocols) = matfield
              end if

              if (x.ne.0) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x-1
                 col(MatStencil_j,nocols) = y
                 col(MatStencil_k,nocols) = z
                 col(MatStencil_c,nocols) = matfield
              end if

              if (x.ne.psx_g-1) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x+1
                 col(MatStencil_j,nocols) = y
                 col(MatStencil_k,nocols) = z
                 col(MatStencil_c,nocols) = matfield
              end if

              if (y.ne.psy_g-1) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x
                 col(MatStencil_j,nocols) = y+1
                 col(MatStencil_k,nocols) = z
                 col(MatStencil_c,nocols) = matfield
              end if

              if (z.ne.psz_g-1) then    ! z-1
                 nocols = nocols + 1
                 col(MatStencil_i,nocols) = x
                 col(MatStencil_j,nocols) = y
                 col(MatStencil_k,nocols) = z+1
                 col(MatStencil_c,nocols) = matfield
              end if

              call MatSetValuesStencil(matprecond,1,row,nocols,col,zeromatentry,INSERT_VALUES,ierr)

           end do
        end do
     end do
  end do

  call MatAssemblyBegin(matprecond,MAT_FLUSH_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FLUSH_ASSEMBLY,ierr)






  do k=simstate%startz,simstate%startz+simstate%widthz-1
     do j=simstate%starty,simstate%starty+simstate%widthy-1
        do i=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = i
           row(MatStencil_j,1) = j
           row(MatStencil_k,1) = k
           row(MatStencil_c,1) = npH

           nocols = 0
           add_to_v_ij = 0.0d0


           if ((k .ne. (psz_g-1)) .and. (statepointer(nenv,i,j,k) .gt. 0.1d0)) then

              D = D_H_env*6E-2

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

           end if

           nocols = nocols + 1
           v(nocols) = (1.0d0/dt) - add_to_v_ij
           col(MatStencil_i,nocols) = i
           col(MatStencil_j,nocols) = j
           col(MatStencil_k,nocols) = k
           col(MatStencil_c,nocols) = npH

           call MatSetValuesStencil(matprecond,1,row,nocols,col,v,INSERT_VALUES,ierr)


        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  call MatAssemblyBegin(matprecond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FINAL_ASSEMBLY,ierr)

  call VecDestroy(statelocal,ierr)
  return
end subroutine ComputeMatrix_pH

