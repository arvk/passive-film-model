subroutine pfsolve(iter)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none
  external pfFunction, pfJacobian

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

  !!---------------PF evolution-----------------!!

  !! Bulk free energy 
  real*8 :: f_met, f_mkw, f_pht, f_pyr, f_env
  real*8 :: w_met, w_mkw, w_pht, w_pyr, w_env

  real*8 :: sigma_pyr_met, sigma_pyr_mkw, sigma_pyr_pht, sigma_pyr_env
  real*8 :: sigma_met_pyr, sigma_mkw_pyr, sigma_pht_pyr, sigma_env_pyr


  integer, parameter :: imet = 1
  integer, parameter :: imkw = 2
  integer, parameter :: ipht = 3
  integer, parameter :: ipyr = 4
  integer, parameter :: ienv = 5

  integer, parameter :: no_fields = 5

  real*8, dimension(psx,psy,psz+2) :: D_met, D_mkw, D_pht, D_pyr, D_env


  ! A/B/JA matrices for implicit solver
  real*8, dimension(psx*psy*psz*no_fields) :: B
  real*8, dimension(psx*psy*psz*no_fields) :: approxsol
  integer, dimension(psx*psy*psz*no_fields) :: vector_locator
  real*8, dimension(psx*psy*psz*no_fields) :: vecread
  real*8, dimension(psx*psy*psz*no_fields) :: scratch1,scratch2,scratch3

  integer :: linindex, contindex
  integer :: iterations, solver_info


  PetscErrorCode ierr
  Vec pf_vec,rhs_vec,ret_vec
  PetscScalar, pointer :: point_pf_vec(:)
  Mat jac
  SNES snes_pf











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if (mod(iter,swap_freq_pf).eq.1) then
        call swap_pf()
     end if

     sigma_pyr_met = sigma_pyr_met_0
     sigma_pyr_mkw = sigma_pyr_mkw_0
     sigma_pyr_pht = sigma_pyr_pht_0
     sigma_pyr_env = sigma_pyr_env_0

     sigma_met_pyr = sigma_pyr_met
     sigma_mkw_pyr = sigma_pyr_mkw
     sigma_pht_pyr = sigma_pyr_pht
     sigma_env_pyr = sigma_pyr_env


     do z = 1,psz+2
        do y = 1,psy
           do x = 1,psx
              D_met(x,y,z) = (M_met_mkw*sigma_met_mkw*mkw(x,y,z))+(M_met_pht*sigma_met_pht*pht(x,y,z))+(M_met_pyr*sigma_met_pyr*pyr(x,y,z))+(M_met_env*sigma_met_env*env(x,y,z))
              D_mkw(x,y,z) = (M_mkw_met*sigma_mkw_met*met(x,y,z))+(M_mkw_pht*sigma_mkw_pht*pht(x,y,z))+(M_mkw_pyr*sigma_mkw_pyr*pyr(x,y,z))+(M_mkw_env*sigma_mkw_env*env(x,y,z))
              D_pht(x,y,z) = (M_pht_met*sigma_pht_met*met(x,y,z))+(M_pht_mkw*sigma_pht_mkw*mkw(x,y,z))+(M_pht_pyr*sigma_pht_pyr*pyr(x,y,z))+(M_pht_env*sigma_pht_env*env(x,y,z))
              D_pyr(x,y,z) = (M_pyr_met*sigma_pyr_met*met(x,y,z))+(M_pyr_mkw*sigma_pyr_mkw*mkw(x,y,z))+(M_pyr_pht*sigma_pyr_pht*pht(x,y,z))+(M_pyr_env*sigma_pyr_env*env(x,y,z))
              D_env(x,y,z) = (M_env_met*sigma_env_met*met(x,y,z))+(M_env_mkw*sigma_env_mkw*mkw(x,y,z))+(M_env_pht*sigma_env_pht*pht(x,y,z))+(M_env_pyr*sigma_env_pyr*pyr(x,y,z))
           end do
        end do
     end do
     


  call VecCreate(PETSC_COMM_SELF,ret_vec,ierr)
  call VecSetSizes(ret_vec,PETSC_DECIDE,psx*psy*psz*no_fields,ierr)
  call VecSetFromOptions(ret_vec,ierr)
  call VecSetUp(ret_vec,ierr)






  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex+((imet-1)*psx*psy*psz)) =  (met(x,y,z+1)/dt)
           B(linindex+((imkw-1)*psx*psy*psz)) =  (mkw(x,y,z+1)/dt)
           B(linindex+((ipht-1)*psx*psy*psz)) =  (pht(x,y,z+1)/dt)
           B(linindex+((ipyr-1)*psx*psy*psz)) =  (pyr(x,y,z+1)/dt)
           B(linindex+((ienv-1)*psx*psy*psz)) =  (env(x,y,z+1)/dt)

           if (z.eq.1) then
              B(linindex+((imet-1)*psx*psy*psz)) = B(linindex+((imet-1)*psx*psy*psz)) + ((0.5d0*(D_met(x,y,z+1)+D_met(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*psz)) = B(linindex+((imkw-1)*psx*psy*psz)) + ((0.5d0*(D_mkw(x,y,z+1)+D_mkw(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*psz)) = B(linindex+((ipht-1)*psx*psy*psz)) + ((0.5d0*(D_pht(x,y,z+1)+D_pht(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*psz)) = B(linindex+((ipyr-1)*psx*psy*psz)) + ((0.5d0*(D_pyr(x,y,z+1)+D_pyr(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*psz)) = B(linindex+((ienv-1)*psx*psy*psz)) + ((0.5d0*(D_env(x,y,z+1)+D_env(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))
           elseif (z.eq.psz) then
              B(linindex+((imet-1)*psx*psy*psz)) = B(linindex+((imet-1)*psx*psy*psz)) + ((0.5d0*(D_met(x,y,z+1+1)+D_met(x,y,z+1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*psz)) = B(linindex+((imkw-1)*psx*psy*psz)) + ((0.5d0*(D_mkw(x,y,z+1+1)+D_mkw(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*psz)) = B(linindex+((ipht-1)*psx*psy*psz)) + ((0.5d0*(D_pht(x,y,z+1+1)+D_pht(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*psz)) = B(linindex+((ipyr-1)*psx*psy*psz)) + ((0.5d0*(D_pyr(x,y,z+1+1)+D_pyr(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*psz)) = B(linindex+((ienv-1)*psx*psy*psz)) + ((0.5d0*(D_env(x,y,z+1+1)+D_env(x,y,z+1))/(dpf*dpf))*env(x,y,z+1-1))
           end if

           approxsol(linindex+((imet-1)*psx*psy*psz)) = met(x,y,z+1)
           approxsol(linindex+((imkw-1)*psx*psy*psz)) = mkw(x,y,z+1)
           approxsol(linindex+((ipht-1)*psx*psy*psz)) = pht(x,y,z+1)
           approxsol(linindex+((ipyr-1)*psx*psy*psz)) = pyr(x,y,z+1)
           approxsol(linindex+((ienv-1)*psx*psy*psz)) = env(x,y,z+1)

           vector_locator(linindex+((imet-1)*psx*psy*psz)) = linindex+((imet-1)*psx*psy*psz)-1
           vector_locator(linindex+((imkw-1)*psx*psy*psz)) = linindex+((imkw-1)*psx*psy*psz)-1
           vector_locator(linindex+((ipht-1)*psx*psy*psz)) = linindex+((ipht-1)*psx*psy*psz)-1
           vector_locator(linindex+((ipyr-1)*psx*psy*psz)) = linindex+((ipyr-1)*psx*psy*psz)-1
           vector_locator(linindex+((ienv-1)*psx*psy*psz)) = linindex+((ienv-1)*psx*psy*psz)-1

        end do
     end do
  end do



  call VecCreate(PETSC_COMM_SELF,pf_vec,ierr)
  call VecSetSizes(pf_vec,PETSC_DECIDE,psx*psy*psz*no_fields,ierr)
  call VecSetFromOptions(pf_vec,ierr)
  call VecSetUp(pf_vec,ierr)
  call VecSetValues(pf_vec,psx*psy*psz*no_fields,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecAssemblyBegin(pf_vec,ierr)
  call VecAssemblyEnd(pf_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*psz*no_fields,ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*psz*no_fields,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_vec,ierr)
  call VecAssemblyEnd(rhs_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac,ierr)
  call MatSetSizes(jac,PETSC_DECIDE,PETSC_DECIDE,psx*psy*psz*no_fields,psx*psy*psz*no_fields,ierr)
  call MatSetUp(jac,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pf,ierr)
  call SNESSetFunction(snes_pf,ret_vec,pfFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pf,jac,jac,pfJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pf,ierr)
  call SNESSolve(snes_pf,rhs_vec,pf_vec,ierr)

  call VecDestroy(pf_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call VecDestroy(ret_vec,ierr)
  call MatDestroy(jac,ierr)
  call SNESDestroy(snes_pf,ierr)


end subroutine pfsolve

