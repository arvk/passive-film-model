subroutine solve_mu(iter,ksp_mu,simstate)
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  !! **Set up and solve the linear equations for time-evolution of sulfur chemical potential in different \(FeS\) phases**
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  PetscErrorCode ierr
  KSP ksp_mu !! Linear chemical potential solver
  KSPConvergedReason mu_converged_reason
  Vec state, solved_mu_vector, state_solved !! Vectors to store function values and solutions
  PetscInt, intent(in) :: iter  !! Current iteration number
  PetscInt :: x, y, z           !! Coordinates inside the simulation system
  type(context) simstate       !! Field variables stored in PETSc vectors and DMDA objects
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
     call VecCreate(MPI_COMM_WORLD,solved_mu_vector,ierr)
     call VecSetSizes(solved_mu_vector,PETSC_DECIDE,psx_g*psy_g*psz_g,ierr)
     call VecSetUp(solved_mu_vector,ierr)
     call VecStrideGather(state_solved,nmus,solved_mu_vector,INSERT_VALUES,ierr)
     call VecStrideScatter(solved_mu_vector,nmus,simstate%slice,INSERT_VALUES,ierr)
     call VecDestroy(solved_mu_vector,ierr)
  else
     write(6,*) 'Chemical potential field evolution did not converge. Reason: ', mu_converged_reason
  end if

end subroutine solve_mu








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
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_mu
  PetscErrorCode ierr
  Vec b
  type(context) simstate

  call VecCopy(simstate%slice,b,ierr)

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
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_mu
  PetscErrorCode ierr
  Vec state, exstate, b
  PetscScalar, pointer :: statepointer(:,:,:,:), exstatepointer(:,:,:,:), bpointer(:,:,:,:)
  type(context) simstate
  PetscInt :: i, j, k, fesphase, field
  PetscScalar :: Chi, S_source_sink

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,simstate%exslice,exstatepointer,ierr)
  call DMDAVecGetArrayF90(simstate%lattval,b,bpointer,ierr)

  if (include_dissolve) then
     sulf_rate(nmet) = sulf_rate_liq(nmet)
     sulf_rate(nmkw) = sulf_rate_liq(nmkw)
     sulf_rate(npht) = sulf_rate_liq(npht)
     sulf_rate(npyr) = sulf_rate_liq(npyr)
  else
     sulf_rate(nmet) = sulf_rate_gas(nmet)
     sulf_rate(nmkw) = sulf_rate_gas(nmkw)
     sulf_rate(npht) = sulf_rate_gas(npht)
     sulf_rate(npyr) = sulf_rate_gas(npyr)
  end if

  do k=simstate%startz,simstate%startz+simstate%widthz-1
     do j=simstate%starty,simstate%starty+simstate%widthy-1
        do i=simstate%startx,simstate%startx+simstate%widthx-1

           do field = 0, (nfields-1)
              bpointer(field,i,j,k) = 0.0d0
           end do

           bpointer(nmus,i,j,k) = statepointer(nmus,i,j,k) * (1.0d0/dt)

           ! Calculate derivative of sulfur concentration with chemical potential
           Chi = 0.0d0
           do fesphase = 0,nphases-1
              Chi = Chi + drho_dmu(fesphase)*statepointer(fesphase,i,j,k)
           end do

           ! Calculate sulfur source and sinks due to phase transformations
           S_source_sink = 0.0d0
           do fesphase = 0,nphases-1
              S_source_sink = S_source_sink + rhoS(fesphase)*(statepointer(fesphase,i,j,k)-exstatepointer(fesphase,i,j,k))
           end do

           bpointer(nmus,i,j,k) = bpointer(nmus,i,j,k) - (S_source_sink/Chi)

           ! INCLUDE SULFIDATION
           if ((statepointer(nenv,i,j,min(k+3,simstate%startz+simstate%widthz-1))-statepointer(nenv,i,j,max(k,simstate%startz))).gt.0.1d0) then
              do fesphase = nmet,npyr
                 bpointer(nmus,i,j,k) = bpointer(nmus,i,j,k) + ((rhoS(max(min(fesphase,npyr),nmkw))-rhoS(max(min(fesphase-1,npyr),nmet)))*sulf_rate(fesphase)/(dpf*drho_dmu(fesphase)))*statepointer(fesphase,i,j,k)
              end do
              bpointer(nmus,i,j,k) = min(bpointer(nmus,i,j,k),avg_mu_env/dt)
           end if

           if ((k.eq.psz_g-1).or.statepointer(nenv,i,j,k).gt.0.97d0) then
              bpointer(nmus,i,j,k) = avg_mu_env * (1.0d0/dt)
           end if

          if (k.eq.0) then
             bpointer(nmus,i,j,k) = (mus_met_mkw_eqb - (R*T*0.5d0)) * (1.0d0/dt)
          end if

        end do
     end do
  end do

  call DMDAVecRestoreArrayF90(simstate%lattval,b,bpointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%exslice,exstatepointer,ierr)
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

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
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>

  KSP ksp_mu
  PetscErrorCode ierr
  Vec statelocal
  Mat matoper, matprecond
  PetscInt     i,j,k,x,y,z
  PetscScalar  v(7)
  MatStencil   row(4,1),col(4,7)
  PetscScalar, pointer :: statepointer(:,:,:,:)
  PetscScalar :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  PetscInt :: nocols
  PetscScalar :: add_to_v_ij
  type(context) simstate
  PetscScalar zeromatentry(7)
  PetscInt :: matfield


  D_inter_met = D_S_met
  D_inter_mkw = max(D_Fe_mkw,D_S_mkw)
  D_inter_pht = max(D_Fe_pht,D_S_pht)
  D_inter_pyr = max(D_Fe_pyr,D_S_pyr)
  D_inter_env = D_S_env

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
           row(MatStencil_c,1) = nmus

           nocols = 0
           add_to_v_ij = 0.0d0

           if ((k.ne.psz_g-1).and.(k.ne.0).and.(statepointer(nenv,i,j,k).lt.0.97d0)) then

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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

           if ((statepointer(nenv,i,j,min(k+2,simstate%startz+simstate%widthz-1))-statepointer(nenv,i,j,max(k,simstate%startz))).lt.0.1d0) then

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
              v(nocols) = (v(nocols)*statepointer(nvoi,i,j,k))/(dpf*dpf)
              add_to_v_ij = add_to_v_ij + v(nocols)
           end if

        end if

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

