subroutine solve_ang(iter,ksp_ang,simstate)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  !! **Set up and solve the linear equations for evolution of pyrite crystal shape**
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
  KSP :: ksp_ang !! Linear pyrite crystal shape solver
  KSPConvergedReason :: ang_converged_reason
  Vec :: solved_ang_vector,state,state_solved !! Vectors to store function values and solutions
  PetscInt, intent(in) :: iter  !! Current iteration number
  PetscInt :: x, y, z           !! Coordinates inside the simulation system
  type(context) simstate       !! Field variables stored in PETSc vectors and DMDA objects
  external computeRHS_ang, computeMatrix_ang, computeInitialGuess_ang

  call KSPSetDM(ksp_ang,simstate%lattval,ierr)
  call KSPSetComputeRHS(ksp_ang,computeRHS_ang,simstate,ierr)
  call KSPSetComputeOperators(ksp_ang,computeMatrix_ang,simstate,ierr)
  call KSPSetComputeInitialGuess(ksp_ang,computeInitialGuess_ang,simstate,ierr)
  call KSPSetOptionsPrefix(ksp_ang,'ang_',ierr)
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

end subroutine solve_ang









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

  KSP :: ksp_ang
  PetscErrorCode :: ierr
  Vec :: b
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

  KSP :: ksp_ang
  PetscErrorCode :: ierr
  Vec :: state, b, onlyang
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

  KSP :: ksp_ang
  PetscErrorCode :: ierr
  Vec :: state,statelocal
  Mat :: matoper, matprecond
  PetscInt :: i, j, k, x, y, z
  PetscScalar :: v(7)
  MatStencil :: row(4,1),col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscScalar :: M_opyr_max = 1.0E3
  PetscScalar :: M_opyr_min = 1.0E-2
  PetscScalar :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  PetscInt :: nocols
  PetscScalar :: add_to_v_ij
  PetscScalar, parameter :: infinitesimal = 1E-15  ! A hard-coded 'small' number
  type(context) simstate
  PetscScalar :: zeromatentry(7)
  PetscInt :: matfield
  PetscScalar :: Mob_ang_field
  PetscScalar :: Diff_ang_field
  PetscScalar :: Del_ang_field

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
           row(MatStencil_c,1) = nang

           nocols = 0
           add_to_v_ij = 0.0d0


           Mob_ang_field = M_opyr_min + statepointer(npyr,i,j,k)*(M_opyr_max-M_opyr_min)

           Del_ang_field = abs(statepointer(nang,i,j,min(k+1,psz_g-1))-statepointer(nang,i,j,max(k-1,0))) + &
                & abs(statepointer(nang,i,min(j+1,psy_g-1),k)-statepointer(nang,i,max(j-1,0),k)) + &
                & abs(statepointer(nang,min(i+1,psx_g-1),j,k)-statepointer(nang,max(i-1,0),j,k))
           Del_ang_field = (Del_ang_field/dpf)+1.0d0

           if (k.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i,j,k-1)*statepointer(npyr,i,j,k-1)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k-1
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (j.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i,j-1,k)*statepointer(npyr,i,j-1,k)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j-1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (i.gt.0) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i-1,j,k)*statepointer(npyr,i-1,j,k)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
              col(MatStencil_i,nocols) = i-1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (i.lt.psx_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i+1,j,k)*statepointer(npyr,i+1,j,k)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
              col(MatStencil_i,nocols) = i+1
              col(MatStencil_j,nocols) = j
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if


           if (j.lt.psy_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i,j+1,k)*statepointer(npyr,i,j+1,k)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
              col(MatStencil_i,nocols) = i
              col(MatStencil_j,nocols) = j+1
              col(MatStencil_k,nocols) = k
              col(MatStencil_c,nocols) = nang
              v(nocols) = v(nocols)/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if (k.lt.psz_g-1) then
              nocols = nocols + 1
              v(nocols) = 0.0d0 - ((0.5d0*gb_S*statepointer(npyr,i,j,k+1)*statepointer(npyr,i,j,k+1)*Mob_ang_field/Del_ang_field)+(0.5d0*gb_S*statepointer(npyr,i,j,k)*statepointer(npyr,i,j,k)*Mob_ang_field/Del_ang_field))
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

  call MatAssemblyBegin(matprecond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(matprecond,MAT_FINAL_ASSEMBLY,ierr)

  call VecDestroy(statelocal,ierr)
  return
end subroutine ComputeMatrix_ang










double precision function angdiff(or1,or2)
  implicit none
  PetscScalar, intent(in) :: or1,or2
  PetscScalar :: Pi = 3.14159265d0

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


