subroutine elpotsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  external electroFunction, electroJacobian

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>




  integer :: x, y, z   ! Loop variables
  real*8 :: correct_rounding     ! Error correction variable used while updating fields
  integer :: status(MPI_STATUS_SIZE)
  integer, intent(in) :: iter

  real*8, dimension(psx,psy,psz+2) :: epsilonr

  real*8, dimension(psx,psy,psz+2) :: newelpot

  ! A/B/JA matrices for implicit solver
  real*8, dimension(psx*psy*psz) :: B
  real*8, dimension(psx*psy*psz) :: approxsol
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8, dimension(psx*psy*psz) :: vecread
  real*8, dimension(psx*psy*psz) :: scratch1,scratch2,scratch3

  integer :: linindex, contindex
  integer :: iterations, solver_info

  real*8 :: epsilon0, epsilon_met, epsilon_mkw, epsilon_pht, epsilon_pyr, epsilon_env
  

  PetscErrorCode ierr
  Vec elpot_vec,rhs_vec,ret_vec
  PetscScalar, pointer :: point_elpot_vec(:)
  Mat jac
  SNES snes_elpot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call swap_electro()

  epsilon_met = 500.0d0
  epsilon_mkw = 500.0d0
  epsilon_pht = 2.0d0
  epsilon_pyr = 11.0d0
  epsilon_env = 80.0d0

  epsilon0 = 8.854187817E-12 !! Define vacuum permittivity

  do z = 1,psz+2
     do y = 1,psy
        do x = 1,psx

           epsilonr(x,y,z) = epsilon_met*met(x,y,z) + epsilon_mkw*mkw(x,y,z) + epsilon_pht*pht(x,y,z) + epsilon_pyr*pyr(x,y,z) + epsilon_env*env(x,y,z) 
           epsilonr(x,y,z) = epsilonr(x,y,z)*epsilon0*(1.0d0-voids(x,y,z))

        end do
     end do
  end do




     


  call VecCreate(PETSC_COMM_SELF,ret_vec,ierr)
  call VecSetSizes(ret_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(ret_vec,ierr)
  call VecSetUp(ret_vec,ierr)






  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  0.0d0

           if (z.eq.1) then
              B(linindex) = B(linindex) - ((0.5d0*(epsilonr(x,y,z+1)+epsilonr(x,y,z+1-1))/(dpf*dpf))*elpot(x,y,z+1-1))
           elseif (z.eq.psz) then
              B(linindex) = B(linindex) - ((0.5d0*(epsilonr(x,y,z+1+1)+epsilonr(x,y,z+1))/(dpf*dpf))*elpot(x,y,z+1+1))
           end if

           approxsol(linindex) = elpot(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do



  call VecCreate(PETSC_COMM_SELF,elpot_vec,ierr)
  call VecSetSizes(elpot_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(elpot_vec,ierr)
  call VecSetUp(elpot_vec,ierr)
  call VecSetValues(elpot_vec,psx*psy*psz,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecAssemblyBegin(elpot_vec,ierr)
  call VecAssemblyEnd(elpot_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_vec,ierr)
  call VecAssemblyEnd(rhs_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac,ierr)
  call MatSetSizes(jac,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz,psx*psy*psz,ierr)
  call MatSetUp(jac,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_elpot,ierr)
  call SNESSetFunction(snes_elpot,ret_vec,electroFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_elpot,jac,jac,electroJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_elpot,ierr)
  call SNESSolve(snes_elpot,rhs_vec,elpot_vec,ierr)


  call VecGetArrayF90(elpot_vec,point_elpot_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           newelpot(x,y,z+1) = point_elpot_vec(linindex)

        end do
     end do
  end do
  call VecRestoreArrayF90(elpot_vec,point_elpot_vec,ierr)


  !! Apply boundary conditions to phase field(s) update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           newelpot(x,y,2) = elpot(x,y,2)
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           newelpot(x,y,psz+1) = elpot(x,y,psz+1)
        end do
     end do
  end if


  !! Apply boundary conditions to the environment
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (met(x,y,z).gt.0.25) then
              newelpot(x,y,z) = -0.50d0
           end if
        end do
     end do
  end do




  !###########################################################
  !##################--UPDATE ALL FIELDS--####################
  !###########################################################

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           elpot(x,y,z) = newelpot(x,y,z)
        end do
     end do
  end do


  call VecDestroy(elpot_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call VecDestroy(ret_vec,ierr)
  call MatDestroy(jac,ierr)
  call SNESDestroy(snes_elpot,ierr)





end subroutine elpotsolve

