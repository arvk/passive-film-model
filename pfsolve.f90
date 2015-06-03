subroutine pfsolve(iter)
  use commondata
  use fields
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

  real*8 :: sumfields

  integer, parameter :: imet = 1
  integer, parameter :: imkw = 2
  integer, parameter :: ipht = 3
  integer, parameter :: ipyr = 4
  integer, parameter :: ienv = 5

  integer, parameter :: no_fields = 5

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_met, D_met_mkw, D_met_pht, D_met_pyr, D_met_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_mkw_met, D_mkw, D_mkw_pht, D_mkw_pyr, D_mkw_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_pht_met, D_pht_mkw, D_pht, D_pht_pyr, D_pht_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_pyr_met, D_pyr_mkw, D_pyr_pht, D_pyr, D_pyr_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_env_met, D_env_mkw, D_env_pht, D_env_pyr, D_env

  ! A/B/JA matrices for implicit solver
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: B
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: approxsol
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: vector_locator
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: vecread
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: scratch1,scratch2,scratch3

  integer :: linindex, contindex
  integer :: iterations, solver_info


  PetscErrorCode ierr
  Vec pf_vec,rhs_vec,ret_vec
  PetscScalar, pointer :: point_pf_vec(:)
  Mat jac
  SNES snes_pf
  SNESConvergedReason pf_converged_reason










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call calc_grad_pf()

  do z = 1,psz+(2*ghost_width)
     do y = 1,psy
        do x = 1,psx

           sigma_pyr_met = sigma_pyr_met_0*(1+0.85d0*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_mkw = sigma_pyr_mkw_0*(1+0.85d0*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.85d0*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.85d0*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_mkw_pyr = sigma_pyr_mkw
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

           D_met_mkw(x,y,z) = 0.0d0 - (M_met_mkw*sigma_met_mkw*mkw(x,y,z))
           D_met_pht(x,y,z) = 0.0d0 - (M_met_pht*sigma_met_pht*pht(x,y,z))
           D_met_pyr(x,y,z) = 0.0d0 - (M_met_pyr*sigma_met_pyr*pyr(x,y,z))
           D_met_env(x,y,z) = 0.0d0 - (M_met_env*sigma_met_env*env(x,y,z))
           D_met(x,y,z) = 0.0d0 - (D_met_mkw(x,y,z)+D_met_pht(x,y,z)+D_met_pyr(x,y,z)+D_met_env(x,y,z))

           D_mkw_met(x,y,z) = 0.0d0 - (M_mkw_met*sigma_mkw_met*met(x,y,z))
           D_mkw_pht(x,y,z) = 0.0d0 - (M_mkw_pht*sigma_mkw_pht*pht(x,y,z))
           D_mkw_pyr(x,y,z) = 0.0d0 - (M_mkw_pyr*sigma_mkw_pyr*pyr(x,y,z))
           D_mkw_env(x,y,z) = 0.0d0 - (M_mkw_env*sigma_mkw_env*env(x,y,z))
           D_mkw(x,y,z) = 0.0d0 - (D_mkw_met(x,y,z)+D_mkw_pht(x,y,z)+D_mkw_pyr(x,y,z)+D_mkw_env(x,y,z))

           D_pht_met(x,y,z) = 0.0d0 - (M_pht_met*sigma_pht_met*met(x,y,z))
           D_pht_mkw(x,y,z) = 0.0d0 - (M_pht_mkw*sigma_pht_mkw*mkw(x,y,z))
           D_pht_pyr(x,y,z) = 0.0d0 - (M_pht_pyr*sigma_pht_pyr*pyr(x,y,z))
           D_pht_env(x,y,z) = 0.0d0 - (M_pht_env*sigma_pht_env*env(x,y,z))
           D_pht(x,y,z) = 0.0d0 - (D_pht_met(x,y,z)+D_pht_mkw(x,y,z)+D_pht_pyr(x,y,z)+D_pht_env(x,y,z))

           D_pyr_met(x,y,z) = 0.0d0 - (M_pyr_met*sigma_pyr_met*met(x,y,z))
           D_pyr_mkw(x,y,z) = 0.0d0 - (M_pyr_mkw*sigma_pyr_mkw*mkw(x,y,z))
           D_pyr_pht(x,y,z) = 0.0d0 - (M_pyr_pht*sigma_pyr_pht*pht(x,y,z))
           D_pyr_env(x,y,z) = 0.0d0 - (M_pyr_env*sigma_pyr_env*env(x,y,z))
           D_pyr(x,y,z) = 0.0d0 - (D_pyr_met(x,y,z)+D_pyr_mkw(x,y,z)+D_pyr_pht(x,y,z)+D_pyr_env(x,y,z))

           D_env_met(x,y,z) = 0.0d0 - (M_env_met*sigma_env_met*met(x,y,z))
           D_env_mkw(x,y,z) = 0.0d0 - (M_env_mkw*sigma_env_mkw*mkw(x,y,z))
           D_env_pht(x,y,z) = 0.0d0 - (M_env_pht*sigma_env_pht*pht(x,y,z))
           D_env_pyr(x,y,z) = 0.0d0 - (M_env_pyr*sigma_env_pyr*pyr(x,y,z))
           D_env(x,y,z) = 0.0d0 - (D_env_met(x,y,z)+D_env_mkw(x,y,z)+D_env_pht(x,y,z)+D_env_pyr(x,y,z))

        end do
     end do
  end do







  call VecCreate(PETSC_COMM_SELF,ret_vec,ierr)
  call VecSetSizes(ret_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2)*no_fields,ierr)
  call VecSetFromOptions(ret_vec,ierr)
  call VecSetUp(ret_vec,ierr)






  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) =  (met(x,y,z+1)/dt)
           B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) =  (mkw(x,y,z+1)/dt)
           B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) =  (pht(x,y,z+1)/dt)
           B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) =  (pyr(x,y,z+1)/dt)
           B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) =  (env(x,y,z+1)/dt)

           if (z.eq.1) then

              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met(x,y,z+1)+D_met(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_met(x,y,z+1)+D_mkw_met(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_met(x,y,z+1)+D_pht_met(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_met(x,y,z+1)+D_pyr_met(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_met(x,y,z+1)+D_env_met(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))

              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_mkw(x,y,z+1)+D_met_mkw(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw(x,y,z+1)+D_mkw(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_mkw(x,y,z+1)+D_pht_mkw(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_mkw(x,y,z+1)+D_pyr_mkw(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_mkw(x,y,z+1)+D_env_mkw(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))

              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_pht(x,y,z+1)+D_met_pht(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_pht(x,y,z+1)+D_mkw_pht(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht(x,y,z+1)+D_pht(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_pht(x,y,z+1)+D_pyr_pht(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_pht(x,y,z+1)+D_env_pht(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))

              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_pyr(x,y,z+1)+D_met_pyr(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_pyr(x,y,z+1)+D_mkw_pyr(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_pyr(x,y,z+1)+D_pht_pyr(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr(x,y,z+1)+D_pyr(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_pyr(x,y,z+1)+D_env_pyr(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))

              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_env(x,y,z+1)+D_met_env(x,y,z+1-1))/(dpf*dpf))*met(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_env(x,y,z+1)+D_mkw_env(x,y,z+1-1))/(dpf*dpf))*mkw(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_env(x,y,z+1)+D_pht_env(x,y,z+1-1))/(dpf*dpf))*pht(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_env(x,y,z+1)+D_pyr_env(x,y,z+1-1))/(dpf*dpf))*pyr(x,y,z+1-1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env(x,y,z+1)+D_env(x,y,z+1-1))/(dpf*dpf))*env(x,y,z+1-1))

           elseif (z.eq.(psz+(2*ghost_width)-2)) then

              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met(x,y,z+1+1)+D_met(x,y,z+1))/(dpf*dpf))*met(x,y,z+1+1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_met(x,y,z+1+1)+D_mkw_met(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1+1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_met(x,y,z+1+1)+D_pht_met(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1+1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_met(x,y,z+1+1)+D_pyr_met(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1+1))
              B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_met(x,y,z+1+1)+D_env_met(x,y,z+1))/(dpf*dpf))*env(x,y,z+1+1))

              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_mkw(x,y,z+1+1)+D_met_mkw(x,y,z+1))/(dpf*dpf))*met(x,y,z+1+1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw(x,y,z+1+1)+D_mkw(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1+1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_mkw(x,y,z+1+1)+D_pht_mkw(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1+1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_mkw(x,y,z+1+1)+D_pyr_mkw(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1+1))
              B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_mkw(x,y,z+1+1)+D_env_mkw(x,y,z+1))/(dpf*dpf))*env(x,y,z+1+1))

              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_pht(x,y,z+1+1)+D_met_pht(x,y,z+1))/(dpf*dpf))*met(x,y,z+1+1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_pht(x,y,z+1+1)+D_mkw_pht(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1+1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht(x,y,z+1+1)+D_pht(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1+1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_pht(x,y,z+1+1)+D_pyr_pht(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1+1))
              B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_pht(x,y,z+1+1)+D_env_pht(x,y,z+1))/(dpf*dpf))*env(x,y,z+1+1))

              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_pyr(x,y,z+1+1)+D_met_pyr(x,y,z+1))/(dpf*dpf))*met(x,y,z+1+1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_pyr(x,y,z+1+1)+D_mkw_pyr(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1+1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_pyr(x,y,z+1+1)+D_pht_pyr(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1+1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr(x,y,z+1+1)+D_pyr(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1+1))
              B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env_pyr(x,y,z+1+1)+D_env_pyr(x,y,z+1))/(dpf*dpf))*env(x,y,z+1+1))

              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_met_env(x,y,z+1+1)+D_met_env(x,y,z+1))/(dpf*dpf))*met(x,y,z+1+1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_mkw_env(x,y,z+1+1)+D_mkw_env(x,y,z+1))/(dpf*dpf))*mkw(x,y,z+1+1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pht_env(x,y,z+1+1)+D_pht_env(x,y,z+1))/(dpf*dpf))*pht(x,y,z+1+1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_pyr_env(x,y,z+1+1)+D_pyr_env(x,y,z+1))/(dpf*dpf))*pyr(x,y,z+1+1))
              B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = B(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) + ((0.5d0*(D_env(x,y,z+1+1)+D_env(x,y,z+1))/(dpf*dpf))*env(x,y,z+1+1))

           end if

           approxsol(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = met(x,y,z+1)
           approxsol(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = mkw(x,y,z+1)
           approxsol(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = pht(x,y,z+1)
           approxsol(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = pyr(x,y,z+1)
           approxsol(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = env(x,y,z+1)

           vector_locator(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))) = linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2))-1
           vector_locator(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))) = linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))-1
           vector_locator(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))) = linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))-1
           vector_locator(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))) = linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))-1
           vector_locator(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))) = linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))-1

        end do
     end do
  end do



  call VecCreate(PETSC_COMM_SELF,pf_vec,ierr)
  call VecSetSizes(pf_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2)*no_fields,ierr)
  call VecSetFromOptions(pf_vec,ierr)
  call VecSetUp(pf_vec,ierr)
  call VecSetValues(pf_vec,psx*psy*(psz+(2*ghost_width)-2)*no_fields,vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecAssemblyBegin(pf_vec,ierr)
  call VecAssemblyEnd(pf_vec,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2)*no_fields,ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*(psz+(2*ghost_width)-2)*no_fields,vector_locator,B,INSERT_VALUES,ierr)

  call VecAssemblyBegin(rhs_vec,ierr)
  call VecAssemblyEnd(rhs_vec,ierr)

  call MatCreate(PETSC_COMM_SELF,jac,ierr)
  call MatSetSizes(jac,PETSC_DECIDE,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2)*no_fields,psx*psy*(psz+(2*ghost_width)-2)*no_fields,ierr)
  call MatSetUp(jac,ierr)

  call SNESCreate(PETSC_COMM_SELF,snes_pf,ierr)
  call SNESSetFunction(snes_pf,ret_vec,pfFunction,PETSC_NULL_OBJECT,ierr)
  call SNESSetJacobian(snes_pf,jac,jac,pfJacobian,PETSC_NULL_OBJECT,ierr)
  call SNESSetFromOptions(snes_pf,ierr)
  call SNESSolve(snes_pf,rhs_vec,pf_vec,ierr)

  call SNESGetConvergedReason(snes_pf,pf_converged_reason,ierr)

  if (pf_converged_reason.gt.0) then

     call VecGetArrayF90(pf_vec,point_pf_vec,ierr)
     do z = 1,psz+(2*ghost_width)-2
        do y = 1,psy
           do x = 1,psx

              linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

              newmet(x,y,z+1) = point_pf_vec(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2)))
              newmkw(x,y,z+1) = point_pf_vec(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2)))
              newpht(x,y,z+1) = point_pf_vec(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2)))
              newpyr(x,y,z+1) = point_pf_vec(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2)))
              newenv(x,y,z+1) = env(x,y,z+1) !point_pf_vec(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))),1.0d0),0.0d0)

           end do
        end do
     end do
     call VecRestoreArrayF90(pf_vec,point_pf_vec,ierr)

  else ! if pf_converged_reason < 0

     newmet = met
     newmkw = mkw
     newpht = pht
     newpyr = pyr
     newenv = env

  end if


  !! Apply boundary conditions to phase field(s) update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           newmet(x,y,ghost_width+1) = met(x,y,ghost_width+1)
           newmkw(x,y,ghost_width+1) = mkw(x,y,ghost_width+1)
           newpht(x,y,ghost_width+1) = pht(x,y,ghost_width+1)
           newpyr(x,y,ghost_width+1) = pyr(x,y,ghost_width+1)
           newenv(x,y,ghost_width+1) = env(x,y,ghost_width+1)
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           newmet(x,y,psz+ghost_width) = met(x,y,psz+ghost_width)
           newmkw(x,y,psz+ghost_width) = mkw(x,y,psz+ghost_width)
           newpht(x,y,psz+ghost_width) = pht(x,y,psz+ghost_width)
           newpyr(x,y,psz+ghost_width) = pyr(x,y,psz+ghost_width)
           newenv(x,y,psz+ghost_width) = env(x,y,psz+ghost_width)
        end do
     end do
  end if


  !! Apply boundary conditions to the environment
  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           if (env(x,y,z).gt.0.99) then
              newmet(x,y,z) = 0.0d0
              newmkw(x,y,z) = 0.0d0
              newpht(x,y,z) = 0.0d0
              newpyr(x,y,z) = 0.0d0
              newenv(x,y,z) = 1.0d0
           end if
        end do
     end do
  end do







!!! Update phase fields


  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           dmet_dt(x,y,z) = ((newmet(x,y,z) - met(x,y,z))/dt)*(1.0d0-voids(x,y,z))
           dmkw_dt(x,y,z) = ((newmkw(x,y,z) - mkw(x,y,z))/dt)*(1.0d0-voids(x,y,z))
           dpht_dt(x,y,z) = ((newpht(x,y,z) - pht(x,y,z))/dt)*(1.0d0-voids(x,y,z))
           dpyr_dt(x,y,z) = ((newpyr(x,y,z) - pyr(x,y,z))/dt)*(1.0d0-voids(x,y,z))
           denv_dt(x,y,z) = ((newenv(x,y,z) - env(x,y,z))/dt)*(1.0d0-voids(x,y,z))
        end do
        do z = 1,ghost_width
           dmet_dt(x,y,z) = 0.0d0 ; dmkw_dt(x,y,z) = 0.0d0 ; dpht_dt(x,y,z) = 0.0d0 ; dpyr_dt(x,y,z) = 0.0d0 ; denv_dt(x,y,z) = 0.0d0
           dmet_dt(x,y,psz+ghost_width+z) = 0.0d0 ; dmkw_dt(x,y,psz+ghost_width+z) = 0.0d0 ; dpht_dt(x,y,psz+ghost_width+z) = 0.0d0 ; dpyr_dt(x,y,psz+ghost_width+z) = 0.0d0 ; denv_dt(x,y,psz+ghost_width+z) = 0.0d0
        end do
     end do
  end do

  !###########################################################
  !##################--UPDATE ALL FIELDS--####################
  !###########################################################

  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           met(x,y,z) = newmet(x,y,z)
           mkw(x,y,z) = newmkw(x,y,z)
           pht(x,y,z) = newpht(x,y,z)
           pyr(x,y,z) = newpyr(x,y,z)
           env(x,y,z) = newenv(x,y,z)
        end do
     end do
  end do


  call VecDestroy(pf_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call VecDestroy(ret_vec,ierr)
  call MatDestroy(jac,ierr)
  call SNESDestroy(snes_pf,ierr)



end subroutine pfsolve

