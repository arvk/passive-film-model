subroutine solve_pot(iter,snes_pot,simstate)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  !! **Set up and solve the Poisson-Boltzmann equations in the environment region on the simulation cell**
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  PetscErrorCode :: ierr
  SNES :: snes_pot !! Non-linear electrical potential solver
  SNESConvergedReason :: pot_converged_reason
  Vec :: vec_feval, elpot_vector, solution_vec,rhs_vec !! Vectors to store function values and solutions
  Mat :: mat_jacob !! Jacobian for electrical potential evolution
  type(context) simstate !! Field variables stored in PETSc vectors and DMDA objects
  PetscInt, intent(in) :: iter  !! Current iteration number
  PetscInt :: fesphase !! Index corresponding to FeS phase
  external FormFunction_pot, FormJacobian_pot


  call DMCreateGlobalVector(simstate%lattval,solution_vec,ierr)
  call VecCopy(simstate%slice,solution_vec,ierr)

  call DMCreateGlobalVector(simstate%lattval,vec_feval,ierr)
  call VecSet(vec_feval,0.0d0,ierr)

  call DMCreateGlobalVector(simstate%lattval,rhs_vec,ierr)
  call FormRHS_pot(rhs_vec,simstate)

  call DMCreateMatrix(simstate%lattval,mat_jacob,ierr)
  call DMSetMatType(simstate%lattval,MATAIJ,ierr)


  call SNESSetJacobian(snes_pot,mat_jacob,mat_jacob,FormJacobian_pot,simstate,ierr)
  call SNESSetFunction(snes_pot,vec_feval,FormFunction_pot,simstate,ierr)
  call SNESSetDM(snes_pot,simstate%lattval,ierr)
  call SNESSetOptionsPrefix(snes_pot,'pot_',ierr)
  call SNESSetFromOptions(snes_pot,ierr)
  call SNESSolve(snes_pot,rhs_vec,solution_vec,ierr)
  call SNESGetConvergedReason(snes_pot,pot_converged_reason,ierr)


  if (pot_converged_reason .ge. 0) then
     call VecCreate(MPI_COMM_WORLD,elpot_vector,ierr)
     call VecSetSizes(elpot_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(elpot_vector,ierr)
     call VecStrideGather(solution_vec,npot,elpot_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(elpot_vector,npot,simstate%slice,INSERT_VALUES,ierr)
     call VecDestroy(elpot_vector,ierr)
  else
     write(6,*) 'Potential field did not converge. Reason: ', pot_converged_reason
  end if


  call MatDestroy(mat_jacob,ierr)
  call VecDestroy(rhs_vec,ierr)
  call VecDestroy(vec_feval,ierr)
  call VecDestroy(solution_vec,ierr)

end subroutine solve_pot









subroutine FormRHS_pot(rhs_vec,simstate)
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

  PetscErrorCode :: ierr
  Vec :: rhs_vec
  PetscScalar, pointer :: statepointer(:,:,:,:), rhspointer(:,:,:,:)
  PetscInt :: startx,starty,startz,widthx,widthy,widthz
  PetscInt :: x,y,z
  type(context) simstate
  PetscInt :: field

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,rhs_vec,rhspointer,ierr)

  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1

           do field = 0,(nfields-1)
              rhspointer(field,x,y,z) = 0.0d0
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
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

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

  SNES :: snes_pot
  PetscErrorCode :: ierr
  Vec :: input_state, function_value, state_local, static_local
  PetscScalar, pointer :: statepointer(:,:,:,:), functionpointer(:,:,:,:), staticpointer(:,:,:,:)
  PetscInt :: fesphase
  PetscInt :: x,y,z
  PetscScalar :: myD
  PetscScalar :: exponent
  PetscScalar :: c0 !! Average density of counter charge (moles/m^3)
  PetscScalar :: el_charg = 1.60217657E-19*6.022E23
  type(context) simstate
  PetscScalar, parameter :: maxconc = 10**(0-1)*1E3
  PetscInt :: myi

  c0 = 10**(0.0d0-pH_in)*1E3

  call DMGetLocalVector(simstate%lattval,state_local,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,input_state,INSERT_VALUES,state_local,ierr)

  call DMGetLocalVector(simstate%lattval,static_local,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,simstate%slice,INSERT_VALUES,static_local,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,simstate%slice,INSERT_VALUES,static_local,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,state_local,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,function_value,functionpointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,static_local,staticpointer,ierr)


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
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y,z-1)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y,z-1)-statepointer(npot,x,y,z))
                 end if

                 if (y.ne.0) then    ! y-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y-1,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y-1,z)-statepointer(npot,x,y,z))
                 end if

                 if (x.ne.0) then    ! x-1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x-1,y,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x-1,y,z)-statepointer(npot,x,y,z))
                 end if

                 if (x.ne.psx_g-1) then    ! x+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x+1,y,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x+1,y,z)-statepointer(npot,x,y,z))
                 end if

                 if (y.ne.psy_g-1) then    ! y+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y+1,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
                    end do
                    functionpointer(npot,x,y,z) = functionpointer(npot,x,y,z) + myD*(statepointer(npot,x,y+1,z)-statepointer(npot,x,y,z))
                 end if

                 if (z.ne.psz_g-1) then    ! z+1
                    myD = 0.0d0
                    do fesphase = nmet,nenv
                       myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y,z+1)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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


  call DMDAVecRestoreArrayF90(simstate%lattval,state_local,statepointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,function_value,functionpointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,static_local,staticpointer,ierr)

  call DMRestoreLocalVector(simstate%lattval,state_local,ierr)
  call DMRestoreLocalVector(simstate%lattval,static_local,ierr)

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

  SNES :: snes_pot
  PetscErrorCode :: ierr
  Vec :: input_state, function_value, state_local, state, static_local
  PetscScalar, pointer :: statepointer(:,:,:,:), staticpointer(:,:,:,:)
  PetscInt :: fesphase
  PetscInt :: x,y,z
  PetscScalar :: myD
  PetscScalar  v((6*2)+1)
  MatStencil   row(4,1),col(4,(6*2)+1)
  PetscInt :: nocols
  Mat pf_jacob, pf_precond
  PetscScalar :: exponent
  PetscScalar :: c0 !! Average density of counter charge (moles/m^3)
  PetscScalar :: el_charg = 1.60217657E-19*6.022E23
  PetscScalar, parameter :: maxconc = 10**(0-1)*1E2
  PetscScalar zeromatentry(7)
  PetscInt :: matfield, matfield2
  type(context) simstate

  c0 = 10**(0.0d0-pH_in)*1E3

  call DMGetLocalVector(simstate%lattval,static_local,ierr)
  call DMGlobalToLocalBegin(simstate%lattval,simstate%slice,INSERT_VALUES,static_local,ierr)
  call DMGlobalToLocalEnd(simstate%lattval,simstate%slice,INSERT_VALUES,static_local,ierr)

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

              do matfield2 = 0,nfields-1   ! Column

                 nocols = 0

                 if (z.ne.0) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z-1
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 if (y.ne.0) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y-1
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 if (x.ne.0) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x-1
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 if (x.ne.psx_g-1) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x+1
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 if (y.ne.psy_g-1) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y+1
                    col(MatStencil_k,nocols) = z
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 if (z.ne.psz_g-1) then    ! z-1
                    nocols = nocols + 1
                    col(MatStencil_i,nocols) = x
                    col(MatStencil_j,nocols) = y
                    col(MatStencil_k,nocols) = z+1
                    col(MatStencil_c,nocols) = matfield2
                 end if

                 call MatSetValuesStencil(pf_precond,1,row,nocols,col,zeromatentry,INSERT_VALUES,ierr)

              end do


           end do
        end do
     end do
  end do

  call MatAssemblyBegin(pf_precond,MAT_FLUSH_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FLUSH_ASSEMBLY,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,input_state,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,static_local,staticpointer,ierr)

  do z=simstate%startz,simstate%startz+simstate%widthz-1
     do y=simstate%starty,simstate%starty+simstate%widthy-1
        do x=simstate%startx,simstate%startx+simstate%widthx-1

           row(MatStencil_i,1) = x
           row(MatStencil_j,1) = y
           row(MatStencil_k,1) = z
           row(MatStencil_c,1) = npot

           if ((staticpointer(nenv,x,y,z).gt.0.97) .and. (z.ne.(psz_g-1)))  then

              nocols = 0

              if (z.ne.0) then    ! z-1
                 myD = 0.0d0
                 do fesphase = nmet,nenv
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y,z-1)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y-1,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x-1,y,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x+1,y,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y+1,z)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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
                    myD = myD + 0.5d0*permittivity(fesphase)*epsilon0*(staticpointer(fesphase,x,y,z+1)+staticpointer(fesphase,x,y,z))/(dpf*dpf)
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

  call MatAssemblyBegin(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FINAL_ASSEMBLY,ierr)

  if (pf_precond .ne. pf_jacob) then
     call MatAssemblyBegin(pf_jacob,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(pf_jacob,MAT_FINAL_ASSEMBLY,ierr)
  end if

  call mpi_barrier(MPI_COMM_WORLD,ierr) ! Barrier before beginning time loop
  call DMRestoreLocalVector(simstate%lattval,static_local,ierr)

  return
end subroutine FormJacobian_pot
