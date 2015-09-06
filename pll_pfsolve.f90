subroutine para_pfsolve(iter,snes_pf)
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
  SNES snes_pf
  SNESConvergedReason pf_converged_reason
  DM da
  Vec vec_feval, single_phase_vector
  Vec state,state_unknown,rhs_vec
  Mat mat_jacob
  PetscInt ctx
  external FormFunction_pf, FormJacobian_pf

  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  integer :: fesphase

  call SNESGetDM(snes_pf,da,ierr)

  call DMGetGlobalVector(da,state,ierr)
  call VecDuplicate(state,vec_feval,ierr)
  call VecDuplicate(state,state_unknown,ierr)
  call VecCopy(state,state_unknown,ierr)
  call DMRestoreGlobalVector(da,state,ierr)

  call DMCreateMatrix(da,mat_jacob,ierr)
  call MatSeqAIJSetPreallocation(mat_jacob,7*5,PETSC_NULL_INTEGER,ierr)  !!! ADD MatMPIAIJSetPreallocation

  call SNESSetFunction(snes_pf,vec_feval,FormFunction_pf,ctx,ierr)
  call SNESSetJacobian(snes_pf,mat_jacob,mat_jacob,FormJacobian_pf,ctx,ierr)
  call FormRHS_pf(state_unknown,rhs_vec)


  call SNESSolve(snes_pf,rhs_vec,state_unknown,ierr)




  call VecCreate(MPI_COMM_WORLD,single_phase_vector,ierr)
  call VecSetSizes(single_phase_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(single_phase_vector,ierr)

  call DMGetGlobalVector(da,state,ierr)
  do fesphase = 0,4     !!!!! FIX THIS
     call VecStrideGather(state_unknown,fesphase,single_phase_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(single_phase_vector,fesphase,state,INSERT_VALUES,ierr)
  end do
  call DMRestoreGlobalVector(da,state,ierr)


!  call VecAXPY(state_unknown,-1.0d0,state,ierr)
!  call VecView(state_unknown,PETSC_NULL_OBJECT,ierr)


end subroutine para_pfsolve









subroutine FormRHS_pf(input_state,rhs_vec)
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

  SNES snes_pf
  PetscInt ctx
  PetscErrorCode ierr
  DM da
  Vec input_state, rhs_vec, state_local
  Vec single_phase_vector
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8, parameter :: mysigma = 2.0d0
  real*8 :: Mobility, hill, S
  real*8 :: w(0:nfields), delo(0:nfields)

  call VecCreate(MPI_COMM_WORLD,single_phase_vector,ierr)
  call VecSetSizes(single_phase_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(single_phase_vector,ierr)


  call VecDuplicate(input_state,rhs_vec,ierr)

  do fesphase = 0,4     !!!!! FIX THIS
     call VecStrideGather(input_state,fesphase,single_phase_vector,INSERT_VALUES,ierr)
     call VecScale(single_phase_vector,(1.0d0/dt),ierr)
     call VecStrideScatter(single_phase_vector,fesphase,rhs_vec,INSERT_VALUES,ierr)
  end do

  call VecAssemblyBegin(rhs_vec,ierr)
  call VecAssemblyEnd(rhs_vec,ierr)
  return
end subroutine FormRHS_pf





























subroutine FormFunction_pf(snes_pf,input_state,function_value,ctx,ierr)
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

  SNES snes_pf
  PetscInt ctx
  PetscErrorCode ierr
  DM da
  Vec input_state, function_value, state_local
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8, parameter :: mysigma = 2.0d0
  real*8 :: Mobility, hill, S
  real*8 :: w(0:nfields), delo(0:nfields)


  w = 0.0d0
  delo = 0.0d0

  call SNESGetDM(snes_pf,da,ierr)
  call DMGetLocalVector(da,state_local,ierr)
  call DMGlobalToLocalBegin(da,input_state,INSERT_VALUES,state_local,ierr)
  call DMGlobalToLocalEnd(da,input_state,INSERT_VALUES,state_local,ierr)

  call DMDAGetCorners(da,startx,starty,startz,widthx,widthy,widthz,ierr)

  call VecDuplicate(input_state,function_value,ierr)
  call VecSetUp(function_value,ierr)

  call DMDAVecGetArrayF90(da,state_local,statepointer,ierr)
  call DMDAVecGetArrayF90(da,function_value,functionpointer,ierr)


  do z=startz,startz+widthz-1
     do y=starty,starty+widthy-1
        do x=startx,startx+widthx-1

           do fesphase = nmet,nenv   ! Row
              do fesphase2 = nmet,nenv  ! Column
               if (fesphase .ne. fesphase2) then
                  sum6myD = 0.0d0

                  if (z.ne.0) then    ! z-1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z-1)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x,y,z-1)-statepointer(fesphase2,x,y,z-1))
                  end if

                  if (y.ne.0) then    ! y-1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y-1,z)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x,y-1,z)-statepointer(fesphase2,x,y-1,z))
                  end if

                  if (x.ne.0) then    ! x-1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x-1,y,z)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x-1,y,z)-statepointer(fesphase2,x-1,y,z))
                  end if

                  if (x.ne.psx_g-1) then    ! x+1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x+1,y,z)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x+1,y,z)-statepointer(fesphase2,x+1,y,z))
                  end if

                  if (y.ne.psy_g-1) then    ! y+1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y+1,z)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x,y+1,z)-statepointer(fesphase2,x,y+1,z))
                  end if

                  if (z.ne.psz_g-1) then    ! z+1
                     myD = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z+1)+statepointer(fesphase,x,y,z))
                     sum6myD = sum6myD - myD
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + myD*(statepointer(fesphase,x,y,z+1)-statepointer(fesphase2,x,y,z+1))
                  end if

                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + sum6myD*(statepointer(fesphase,x,y,z)-statepointer(fesphase2,x,y,z))

                     !! LINEAR
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase2,x,y,z)*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (6.0d0*Mob_pf(fesphase,fesphase2)*hill*(w(fesphase)-w(fesphase2))*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (S*Mob_pf(fesphase,fesphase2)*hill*(delo(fesphase)-delo(fesphase2))*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)

                     !! QUADRATIC
                     functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) - (2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)*statepointer(fesphase,x,y,z)

              end if
              end do
           end do

        end do
     end do
  end do


  call DMDAVecRestoreArrayF90(da,state_local,statepointer,ierr)
  call DMDAVecRestoreArrayF90(da,function_value,functionpointer,ierr)

  call VecAssemblyBegin(function_value,ierr)
  call VecAssemblyEnd(function_value,ierr)
  return
end subroutine FormFunction_pf



subroutine FormJacobian_pf(snes_pf,input_state,pf_jacob,pf_precond,ctx,ierr)
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

  SNES snes_pf
  PetscInt ctx
  PetscErrorCode ierr
  DM da
  Vec input_state, function_value, state_local
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  real*8, allocatable :: D_pf(:,:,:,:,:)
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8, parameter :: mysigma = 2.0d0
  real*8 :: Mobility, hill, S
  real*8 :: w(0:nfields), delo(0:nfields)
  PetscScalar  v((9*2*(nfields-1))+1)
  MatStencil   row(4,1),col(4,(9*2*(nfields-1))+1)
  integer :: nocols
  Mat pf_jacob, pf_precond

  w = 0.0d0
  delo = 0.0d0


  call SNESGetDM(snes_pf,da,ierr)
  call DMGetLocalVector(da,state_local,ierr)
  call DMGlobalToLocalBegin(da,input_state,INSERT_VALUES,state_local,ierr)
  call DMGlobalToLocalEnd(da,input_state,INSERT_VALUES,state_local,ierr)
  call DMDAGetCorners(da,startx,starty,startz,widthx,widthy,widthz,ierr)
  call DMDAVecGetArrayF90(da,state_local,statepointer,ierr)


  do z=startz,startz+widthz-1
     do y=starty,starty+widthy-1
        do x=startx,startx+widthx-1

           do fesphase = nmet,nenv   ! Row

           row(MatStencil_i,1) = x
           row(MatStencil_j,1) = y
           row(MatStencil_k,1) = z
           row(MatStencil_c,1) = fesphase

           nocols = 0

              do fesphase2 = nmet,nenv  ! Column

               if (fesphase .ne. fesphase2) then

                  sum6myD = 0.0d0

                  if (z.ne.0) then    ! z-1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z-1)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z-1
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z-1)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z-1
                     col(MatStencil_c,nocols) = fesphase2
                  end if

                  if (y.ne.0) then    ! y-1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y-1,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y-1
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y-1,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y-1
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2
                  end if

                  if (x.ne.0) then    ! x-1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x-1,y,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x-1
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x-1,y,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x-1
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2
                  end if

                  if (x.ne.psx_g-1) then    ! x+1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x+1,y,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x+1
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x+1,y,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x+1
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2
                  end if

                  if (y.ne.psy_g-1) then    ! y+1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y+1,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y+1
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y+1,z)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y+1
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2
                  end if

                  if (z.ne.psz_g-1) then    ! z+1
                     nocols = nocols + 1
                     v(nocols) = 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z+1)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z+1
                     col(MatStencil_c,nocols) = fesphase

                     sum6myD = sum6myD - v(nocols)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - 0.5d0*Mob_pf(fesphase,fesphase2)*mysigma*(statepointer(fesphase,x,y,z+1)+statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z+1
                     col(MatStencil_c,nocols) = fesphase2
                  end if


                     nocols = nocols + 1
                     v(nocols) = sum6myD
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase


                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - sum6myD
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2


                     ! !! LINEAR
                     ! functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (2.0d0*Mobility*hill*statepointer(fesphase2,x,y,z)*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)
                     ! functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (6.0d0*Mobility*hill*(w(fesphase)-w(fesphase2))*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)
                     ! functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) + (S*Mobility*hill*(delo(fesphase)-delo(fesphase2))*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)

                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - (2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase,x,y,z)*statepointer(fesphase,x,y,z)) + &
                          & (6.0d0*Mob_pf(fesphase,fesphase2)*(w(fesphase)-w(fesphase2))*statepointer(fesphase,x,y,z)) + &
                          & (S*Mob_pf(fesphase,fesphase2)*(delo(fesphase)-delo(fesphase2))*statepointer(fesphase,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2


                     nocols = nocols + 1
                     v(nocols) = (2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase2,x,y,z)*statepointer(fesphase2,x,y,z)) + &
                          & (6.0d0*Mob_pf(fesphase,fesphase2)*(w(fesphase)-w(fesphase2))*statepointer(fesphase2,x,y,z)) + &
                          & (S*Mob_pf(fesphase,fesphase2)*(delo(fesphase)-delo(fesphase2))*statepointer(fesphase2,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     ! !! QUADRATIC
                     ! functionpointer(fesphase,x,y,z) = functionpointer(fesphase,x,y,z) - (2.0d0*Mobility*hill*statepointer(fesphase2,x,y,z))*statepointer(fesphase,x,y,z)*statepointer(fesphase,x,y,z)

                     nocols = nocols + 1
                     v(nocols) = (2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase,x,y,z))*2.0d0*statepointer(fesphase2,x,y,z)
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase2


                     nocols = nocols + 1
                     v(nocols) = 0.0d0 - ((2.0d0*Mob_pf(fesphase,fesphase2)*hill*statepointer(fesphase,x,y,z))*2.0d0*statepointer(fesphase2,x,y,z))
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase


              end if
              end do


                     nocols = nocols + 1
                     v(nocols) = 1.0d0/dt
                     col(MatStencil_i,nocols) = x
                     col(MatStencil_j,nocols) = y
                     col(MatStencil_k,nocols) = z
                     col(MatStencil_c,nocols) = fesphase

                     call MatSetValuesStencil(pf_precond,1,row,nocols,col,v,ADD_VALUES,ierr)

           end do

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(da,state_local,statepointer,ierr)

  call MatAssemblyBegin(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  return
end subroutine FormJacobian_pf

