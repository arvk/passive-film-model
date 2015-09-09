subroutine para_potsolve(iter,snes_pot,simstate)
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
  SNES snes_pot
  SNESConvergedReason pot_converged_reason
  Vec vec_feval, elpot_vector
  Vec state,state_unknown,rhs_vec
  Mat mat_jacob
  type(context) simstate
  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  integer :: fesphase
  real*8 :: el_charg = 1.60217657E-19*6.022E23
  external FormFunction_pot, FormJacobian_pot

  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecDuplicate(state,vec_feval,ierr)
  call VecDuplicate(state,state_unknown,ierr)
  call VecDuplicate(state,rhs_vec,ierr)
  call VecCopy(state,state_unknown,ierr)

  call DMCreateMatrix(simstate%lattval,mat_jacob,ierr)

  call DMCreateGlobalVector(simstate%lattval,vec_feval,ierr)
  call VecSet(vec_feval,0.0d0,ierr)

  call SNESSetJacobian(snes_pot,mat_jacob,mat_jacob,FormJacobian_pot,simstate,ierr)
  call SNESSetFunction(snes_pot,vec_feval,FormFunction_pot,simstate,ierr)
  call FormRHS_pot(simstate,state,rhs_vec)

  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call SNESSetDM(snes_pot,simstate%lattval,ierr)
  call SNESSetFromOptions(snes_pot,ierr)

  call SNESSolve(snes_pot,rhs_vec,state_unknown,ierr)
  call SNESGetConvergedReason(snes_pot,pot_converged_reason,ierr)

  if (pot_converged_reason .ge. 0) then
  call DMGetGlobalVector(simstate%lattval,state,ierr)
  call VecCreate(MPI_COMM_WORLD,elpot_vector,ierr)
  call VecSetSizes(elpot_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
  call VecSetUp(elpot_vector,ierr)
  call VecStrideGather(state_unknown,npot,elpot_vector,INSERT_VALUES,ierr)
  call VecStrideScatter(elpot_vector,npot,state,INSERT_VALUES,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)
  else
     write(6,*) 'Potential field did not converge. Reason: ', pot_converged_reason
  end if

  call MatDestroy(mat_jacob,ierr)
  call VecDestroy(vec_feval,ierr)
  call VecDestroy(elpot_vector,ierr)
end subroutine para_potsolve









subroutine FormRHS_pot(simstate,input_state,rhs_vec)
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

  SNES snes_pot
  PetscErrorCode ierr
  Vec input_state, rhs_vec, state_local
  Vec single_phase_vector
  PetscScalar, pointer :: statepointer(:,:,:,:), rhspointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8 :: Mobility
  real*8 :: w(0:nfields), delo(0:nfields)
  type(context) simstate
  integer :: myi
  MatNullSpace nullspace

  call DMDAVecGetArrayF90(simstate%lattval,input_state,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,rhs_vec,rhspointer,ierr)

  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1

           do myi = 0,(nfields-1)
              rhspointer(myi,x,y,z) = 0.0d0
           end do

           if (statepointer(nenv,x,y,z).gt.0.97) then
              rhspointer(npot,x,y,z) = 0.0d0
           else
              rhspointer(npot,x,y,z) = metal_potential
           end if

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,rhs_vec,rhspointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,input_state,statepointer,ierr)

  call VecAssemblyBegin(rhs_vec,ierr)
  call VecAssemblyEnd(rhs_vec,ierr)

  return
end subroutine FormRHS_pot





























subroutine FormFunction_pot(snes_pot,input_state,function_value,simstate,ierr)
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

  SNES snes_pot
  PetscErrorCode ierr
  Vec input_state, function_value, state_local, state
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:), staticpointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8 :: Mobility
  real*8 :: w(0:nfields), delo(0:nfields)
  real*8 :: grady, gradz
  real*8 :: exponent
  real*8 :: c0 = 1.0d0 !! Moles/m^3
  real*8 :: el_charg = 1.60217657E-19*6.022E23
  type(context) simstate
  real*8, parameter :: maxconc = 1E2
  real*8 :: min
  MatNullSpace nullspace
  integer :: myi

  w = 0.0d0
  delo = 0.0d0

  call DMGetLocalVector(simstate%lattval,state_local,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)

  call DMGetGlobalVector(simstate%lattval,state,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,state,staticpointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,state_local,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,function_value,functionpointer,ierr)


  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1


           do myi = 0,(nfields-1)
              functionpointer(myi,x,y,z) = 0.0d0
           end do


           if (statepointer(nenv,x,y,z).gt.0.97) then

              if (z.eq.psz_g-1) then    ! z+1

                 functionpointer(npot,x,y,z) = statepointer(npot,x,y,z)

              else

                 if (z.ne.0) then    ! z-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y,z-1)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y,z-1)-statepointer(npot,x,y,z))
                 end if

                 if (y.ne.0) then    ! y-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y-1,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y-1,z)-statepointer(npot,x,y,z))
                 end if

                 if (x.ne.0) then    ! x-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x-1,y,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x-1,y,z)-statepointer(npot,x,y,z))
                 end if

                 if (x.ne.psx_g-1) then    ! x+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x+1,y,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x+1,y,z)-statepointer(npot,x,y,z))
                 end if

                 if (y.ne.psy_g-1) then    ! y+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y+1,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y+1,z)-statepointer(npot,x,y,z))
                 end if

                 if (z.ne.psz_g-1) then    ! z+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y,z+1)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y,z+1)-statepointer(npot,x,y,z))
                 end if

                 exponent = (statepointer(npot,x,y,z)*96485)/(R*T)
                 functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + c0*el_charg*(exp(0.0d0-exponent)-exp(0.0d0+exponent))
              end if

           else

              functionpointer(npot,x,y,z) = statepointer(npot,x,y,z)

           end if




        end do
     end do
  end do

  call DMRestoreGlobalVector(simstate%lattval,state,ierr)

  call DMDAVecRestoreArrayF90(simstate%lattval,state_local,statepointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,function_value,functionpointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,state,staticpointer,ierr)

  call DMRestoreLocalVector(simstate%lattval,state_local,ierr)

  call VecAssemblyBegin(function_value,ierr)
  call VecAssemblyEnd(function_value,ierr)

  return
end subroutine FormFunction_pot




















subroutine FormJacobian_pot(snes_pot,input_state,pf_jacob,pf_precond,simstate,ierr)
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

  SNES snes_pot
  PetscErrorCode ierr
  Vec input_state, function_value, state_local, state
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:), staticpointer(:,:,:,:)
  integer :: startx,starty,startz,widthx,widthy,widthz
  real*8, allocatable :: D_pot(:,:,:,:,:)
  integer :: fesphase,fesphase2
  integer :: x,y,z
  real*8 :: myD, sum6myD
  real*8, parameter :: myM = 1.0d0
  real*8 :: Mobility
  real*8 :: w(0:nfields), delo(0:nfields)
  PetscScalar  v((6*2)+1)
  MatStencil   row(4,1),col(4,(6*2)+1)
  integer :: nocols
  Mat pf_jacob, pf_precond
  real*8 :: grady, gradz
  real*8 :: exponent
  real*8 :: c0 = 1.0d0 !! Moles/m^3
  real*8 :: el_charg = 1.60217657E-19*6.022E23
  type(context) simstate
  real*8, parameter :: maxconc = 1E2
  real*8 :: min
  MatNullSpace nullspace


!  call DMCreateMatrix(simstate%lattval,pf_jacob,ierr)

  call DMGetLocalVector(simstate%lattval,state_local,ierr)
  call DMGetGlobalVector(simstate%lattval,state,ierr)

  call DMGlobalToLocalBegin(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,state,staticpointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,state_local,statepointer,ierr)



  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = x
           row(MatStencil_j,1) = y
           row(MatStencil_k,1) = z
           row(MatStencil_c,1) = npot

           v = 0

           nocols = 0
           if (z.ne.0) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x
              col(MatStencil_j,nocols) = y
              col(MatStencil_k,nocols) = z-1
              col(MatStencil_c,nocols) = npot
           end if

           if (y.ne.0) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x
              col(MatStencil_j,nocols) = y-1
              col(MatStencil_k,nocols) = z
              col(MatStencil_c,nocols) = npot
           end if

           if (x.ne.0) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x-1
              col(MatStencil_j,nocols) = y
              col(MatStencil_k,nocols) = z
              col(MatStencil_c,nocols) = npot
           end if

           if (x.ne.psx_g-1) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x+1
              col(MatStencil_j,nocols) = y
              col(MatStencil_k,nocols) = z
              col(MatStencil_c,nocols) = npot
           end if

           if (y.ne.psy_g-1) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x
              col(MatStencil_j,nocols) = y+1
              col(MatStencil_k,nocols) = z
              col(MatStencil_c,nocols) = npot
           end if

           if (z.ne.psz_g-1) then    ! z-1
              nocols = nocols + 1
              col(MatStencil_i,nocols) = x
              col(MatStencil_j,nocols) = y
              col(MatStencil_k,nocols) = z+1
              col(MatStencil_c,nocols) = npot
           end if

           nocols = nocols + 1
           col(MatStencil_i,nocols) = x
           col(MatStencil_j,nocols) = y
           col(MatStencil_k,nocols) = z
           col(MatStencil_c,nocols) = npot

           call MatSetValuesStencil(pf_precond,1,row,nocols,col,v,INSERT_VALUES,ierr)

        end do
     end do
  end do

  call MatAssemblyBegin(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FINAL_ASSEMBLY,ierr)


  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = x
           row(MatStencil_j,1) = y
           row(MatStencil_k,1) = z
           row(MatStencil_c,1) = npot

           if (statepointer(nenv,x,y,z).gt.0.97) then

              if (z.eq.psz_g-1) then    ! z-1

                 v(1) = 1.0d0
                 col(MatStencil_i,1) = x
                 col(MatStencil_j,1) = y
                 col(MatStencil_k,1) = z
                 col(MatStencil_c,1) = npot

                 call MatSetValuesStencil(pf_precond,1,row,1,col,v,ADD_VALUES,ierr)

              else

                 nocols = 0

                 if (z.ne.0) then    ! z-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y,z-1)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z-1
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 if (y.ne.0) then    ! y-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y-1,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y-1
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 if (x.ne.0) then    ! x-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x-1,y,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x-1
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 if (x.ne.psx_g-1) then    ! x+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x+1,y,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x+1
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 if (y.ne.psy_g-1) then    ! y+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y+1,z)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y+1
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 if (z.ne.psz_g-1) then    ! z+1
                    myD = 0.0d0
                   do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(statepointer(fesphase,x,y,z+1)+statepointer(fesphase,x,y,z))/(dpf*dpf)
                    end do

                    nocols = nocols + 1
                    v(nocols) = myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z+1
                    col(MatStencil_c,nocols) = npot

                    nocols = nocols + 1
                    v(nocols) = 0.0d0 - myD
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = npot
                 end if

                 nocols = nocols + 1
                 exponent = (statepointer(npot,x,y,z)*96485)/(R*T)
                 v(nocols) = 0.0d0 - (c0*el_charg*(exp(0.0d0-exponent)+exp(0.0d0+exponent))*(96485/(R*T)))
                 col(MatStencil_i,nocols) = x
                 col(MatStencil_j,nocols) = y
                 col(MatStencil_k,nocols) = z
                 col(MatStencil_c,nocols) = npot

                 call MatSetValuesStencil(pf_precond,1,row,nocols,col,v,ADD_VALUES,ierr)

              end if

           else

              v(1) = 1.0d0
              col(MatStencil_i,1) = x
              col(MatStencil_j,1) = y
              col(MatStencil_k,1) = z
              col(MatStencil_c,1) = npot

              call MatSetValuesStencil(pf_precond,1,row,1,col,v,ADD_VALUES,ierr)

           end if


        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,state_local,statepointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,state,staticpointer,ierr)


  call DMRestoreLocalVector(simstate%lattval,state_local,ierr)
  call DMRestoreGlobalVector(simstate%lattval,state,ierr)


  call MatAssemblyBegin(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FINAL_ASSEMBLY,ierr)

  if (pf_precond .ne. pf_jacob) then
     call MatAssemblyBegin(pf_jacob,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(pf_jacob,MAT_FINAL_ASSEMBLY,ierr)
  end if

  return
end subroutine FormJacobian_pot
