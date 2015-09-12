subroutine para_angsolve(iter,ksp_ang,simstate)
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
  KSP ksp_ang
  KSPConvergedReason ang_converged_reason
  Vec solved_ang_vector,state,state_solved
  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  type(context) simstate
  external computeRHS_ang, computeMatrix_ang, computeInitialGuess_ang

  call KSPSetDM(ksp_ang,simstate%lattval,ierr)
  call KSPSetComputeRHS(ksp_ang,computeRHS_ang,simstate,ierr)
  call KSPSetComputeOperators(ksp_ang,computeMatrix_ang,simstate,ierr)
  call KSPSetComputeInitialGuess(ksp_ang,computeInitialGuess_ang,simstate,ierr)

  call KSPSetFromOptions(ksp_ang,ierr)
  call KSPSolve(ksp_ang,PETSC_NULL_OBJECT,PETSC_NULL_OBJECT,ierr)
  call KSPGetSolution(ksp_ang,state_solved,ierr)
  call KSPGetConvergedReason(ksp_ang,ang_converged_reason,ierr)

  if (ang_converged_reason .gt. 0) then
     call VecCreate(MPI_COMM_WORLD,solved_ang_vector,ierr)
     call VecSetSizes(solved_ang_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(solved_ang_vector,ierr)
     call VecStrideGather(state_solved,nang,solved_ang_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(solved_ang_vector,nang,simstate%slice,INSERT_VALUES,ierr)
     call VecDestroy(solved_ang_vector,ierr)
  else
     write(6,*) 'Orientation field evolution did not converge. Reason: ', ang_converged_reason
  end if

end subroutine para_angsolve









subroutine computeInitialGuess_ang(ksp_ang,b,simstate,ierr)
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

  KSP ksp_ang
  PetscErrorCode ierr
  Vec b
  type(context) simstate

  call VecCopy(simstate%slice,b,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)
  return
end subroutine computeInitialGuess_ang










subroutine computeRHS_ang(ksp_ang,b,simstate,ierr)
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

  KSP ksp_ang
  PetscErrorCode ierr
  Vec state, b, onlyang
  type(context) simstate

  call VecCreate(MPI_COMM_WORLD,onlyang,ierr)
  call VecSetSizes(onlyang,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(onlyang,ierr)
  call VecStrideGather(simstate%slice,nang,onlyang,INSERT_VALUES,ierr)
  call VecScale(onlyang,(1.0d0/dt),ierr)
  call VecStrideScatter(onlyang,nang,b,INSERT_VALUES,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

  ! orc = odiff(opyr(x,y,z+1),opyr(wrap(x-1,psx),y,z+1))-(opyr(x,y,z+1)-opyr(wrap(x-1,psx),y,z+1))
  ! orc1 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(wrap(x-1,psx),y,z))/(dpf*dpf)

  ! orc = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(x,y,z+1))-(opyr(wrap(x+1,psx),y,z+1)-opyr(x,y,z+1))
  ! orc2 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z))/(dpf*dpf)

  ! orc = odiff(opyr(x,y,z+1),opyr(x,wrap(y-1,psy),z+1))-(opyr(x,y,z+1)-opyr(x,wrap(y-1,psy),z+1))
  ! orc3 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,wrap(y-1,psy),z))/(dpf*dpf)

  ! orc = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,y,z+1))-(opyr(x,wrap(y+1,psy),z+1)-opyr(x,y,z+1))
  ! orc4 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z))/(dpf*dpf)

  ! orc = odiff(opyr(x,y,z+1),opyr(x,y,z+1-1))-(opyr(x,y,z+1)-opyr(x,y,z+1-1))
  ! orc5 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf)

  ! orc = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1))-(opyr(x,y,z+1+1)-opyr(x,y,z+1))
  ! orc6 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,z))/(dpf*dpf)

  ! orc = (orc1+orc3+orc5)-(orc2+orc4+orc6)

  call VecDestroy(onlyang,ierr)
 return
end subroutine computeRHS_ang






subroutine ComputeMatrix_ang(ksp_ang,matoper,matprecond,simstate,ierr)
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

  KSP ksp_ang
  PetscErrorCode ierr
  Vec state,statelocal
  Mat matoper, matprecond
  PetscInt     i,j,k
  PetscScalar  v(7)
  MatStencil   row(4,1),col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  real*8 :: D_opyr_max = 1E1   ! Truncation for the orientation field
  real*8 :: M_opyr_max = 1.0E-26
  real*8 :: M_opyr_min = 1.0E-29
  real*8 :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  integer :: nocols, nox, noy, noz
  real*8 :: add_to_v_ij
  real*8, parameter :: infinitesimal = 1E-15  ! A hard-coded 'small' number
  type(context) simstate

  call DMGetLocalVector(simstate%lattval,statelocal,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,simstate%slice,INSERT_VALUES,statelocal,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,simstate%slice,INSERT_VALUES,statelocal,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  do k=simstate%startz,simstate%startz+simstate%widthz-1
     do j=simstate%starty,simstate%starty+simstate%widthy-1
        do i=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = i
           row(MatStencil_j,1) = j
           row(MatStencil_k,1) = k
           row(MatStencil_c,1) = nang

           nocols = 0
           nox = 0
           noy = 0
           noz = 0
           add_to_v_ij = 0.0d0


           if (k.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i,j,k-1))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k-1
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (j.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i,j-1,k))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j-1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (i.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i-1,j,k))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i-1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (i.lt.psx_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i+1,j,k))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i+1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (j.lt.psy_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i,j+1,k))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j+1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (k.lt.psz_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - (0.5d0*D_opyr_max*(min(M_opyr_min/(((statepointer(nang,i,j,k+1))**4)+infinitesimal),M_opyr_max) + min(M_opyr_min/(((statepointer(nang,i,j,k))**4)+infinitesimal),M_opyr_max))/(dpf*dpf))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k+1
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

              nocols = nocols + 1
              v(nocols) = (1.0d0/dt)
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = V(nocols) - add_to_v_ij

              call MatSetValuesStencil(matprecond,1,row,nocols,col,v,ADD_VALUES,ierr)

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,statelocal,statepointer,ierr)

  call DMRestoreLocalVector(simstate%lattval,statelocal,ierr)

  call MatAssemblyBegin(matprecond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FINAL_ASSEMBLY,ierr)

  return
end subroutine ComputeMatrix_ang










double precision function angdiff(or1,or2)
  implicit none
  real*8, intent(in) :: or1,or2
  real*8 :: Pi = 3.14159265d0

  if (or1>or2) then

     if ((or1-or2).lt.(or2+(Pi/2.0d0)-or1)) then
        angdiff = -(or1-or2)
     else
        angdiff = (or2+(Pi/2.0d0)-or1)
     end if

  else

     if ((or2-or1).lt.(or1+(Pi/2.0d0)-or2)) then
        angdiff = or2-or1
     else
        angdiff = -(or1+(Pi/2.0d0)-or2)
     end if

  end if

  return
end function angdiff


