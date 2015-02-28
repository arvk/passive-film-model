subroutine pfFunction(snes,pf_vec,ret_vec,dummy,ierr)
  use commondata
  use fields
  use laplacians
  use gradients
  use thermo_constants
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>

  SNES snes
  Vec pf_vec,ret_vec
  Mat lhs_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pf_vec(:)

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8, dimension(psx*psy*psz) :: lnr,sqr

  real*8, dimension(:), allocatable :: A
  integer, dimension(:), allocatable :: JA, IA

  integer :: x, y, z   ! Loop variables
  integer :: linindex, contindex
  integer :: wrap

  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap


  !! Bulk free energy 
  real*8 :: f_pht, f_env, f_met, f_pyr
  real*8 :: w_pht, w_env, w_met, w_pyr


  real*8 :: sigma_pyr_met, sigma_pyr_mkw, sigma_pyr_pht, sigma_pyr_env
  real*8 :: sigma_met_pyr, sigma_mkw_pyr, sigma_pht_pyr, sigma_env_pyr

  real*8, dimension(psx,psy,psz+2) :: D_met, D_met_mkw, D_met_pht, D_met_pyr, D_met_env
  real*8, dimension(psx,psy,psz+2) :: D_mkw_met, D_mkw, D_mkw_pht, D_mkw_pyr, D_mkw_env
  real*8, dimension(psx,psy,psz+2) :: D_pht_met, D_pht_mkw, D_pht, D_pht_pyr, D_pht_env
  real*8, dimension(psx,psy,psz+2) :: D_pyr_met, D_pyr_mkw, D_pyr_pht, D_pyr, D_pyr_env
  real*8, dimension(psx,psy,psz+2) :: D_env_met, D_env_mkw, D_env_pht, D_env_pyr, D_env

  real*8, dimension(psx,psy,psz+2) :: loc_met, loc_mkw, loc_pht, loc_pyr, loc_env


  integer, parameter :: imet = 1
  integer, parameter :: imkw = 2
  integer, parameter :: ipht = 3
  integer, parameter :: ipyr = 4
  integer, parameter :: ienv = 5

  integer, parameter :: no_fields = 5






     sigma_pyr_met = sigma_pyr_met_0
     sigma_pyr_mkw = sigma_pyr_mkw_0
     sigma_pyr_pht = sigma_pyr_pht_0
     sigma_pyr_env = sigma_pyr_env_0

     sigma_met_pyr = sigma_pyr_met
     sigma_mkw_pyr = sigma_pyr_mkw
     sigma_pht_pyr = sigma_pyr_pht
     sigma_env_pyr = sigma_pyr_env





  call VecGetArrayF90(pf_vec,point_pf_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           loc_met(x,y,z+1) = point_pf_vec(linindex+((imet-1)*psx*psy*psz))
           loc_mkw(x,y,z+1) = point_pf_vec(linindex+((imkw-1)*psx*psy*psz))
           loc_pht(x,y,z+1) = point_pf_vec(linindex+((ipht-1)*psx*psy*psz))
           loc_pyr(x,y,z+1) = point_pf_vec(linindex+((ipyr-1)*psx*psy*psz))
           loc_env(x,y,z+1) = point_pf_vec(linindex+((ienv-1)*psx*psy*psz))

        end do
     end do
  end do
  call VecRestoreArrayF90(pf_vec,point_pf_vec,ierr)


  do y = 1,psy
     do x = 1,psx
        loc_met(x,y,1) = loc_met(x,y,1) ; loc_met(x,y,psz+2) = loc_met(x,y,psz+1) 
        loc_mkw(x,y,1) = loc_mkw(x,y,1) ; loc_mkw(x,y,psz+2) = loc_mkw(x,y,psz+1) 
        loc_pht(x,y,1) = loc_pht(x,y,1) ; loc_pht(x,y,psz+2) = loc_pht(x,y,psz+1) 
        loc_pyr(x,y,1) = loc_pyr(x,y,1) ; loc_pyr(x,y,psz+2) = loc_pyr(x,y,psz+1) 
        loc_env(x,y,1) = loc_env(x,y,1) ; loc_env(x,y,psz+2) = loc_env(x,y,psz+1) 
     end do
  end do


  do z = 1,psz+2
     do y = 1,psy
        do x = 1,psx

           D_met_mkw(x,y,z) = (M_met_mkw*sigma_met_mkw*loc_mkw(x,y,z))
           D_met_pht(x,y,z) = (M_met_pht*sigma_met_pht*loc_pht(x,y,z))
           D_met_pyr(x,y,z) = (M_met_pyr*sigma_met_pyr*loc_pyr(x,y,z))
           D_met_env(x,y,z) = (M_met_env*sigma_met_env*loc_env(x,y,z))
           D_met(x,y,z) = D_met_mkw(x,y,z)+D_met_pht(x,y,z)+D_met_pyr(x,y,z)+D_met_env(x,y,z)

           D_mkw_met(x,y,z) = (M_mkw_met*sigma_mkw_met*loc_met(x,y,z))
           D_mkw_pht(x,y,z) = (M_mkw_pht*sigma_mkw_pht*loc_pht(x,y,z))
           D_mkw_pyr(x,y,z) = (M_mkw_pyr*sigma_mkw_pyr*loc_pyr(x,y,z))
           D_mkw_env(x,y,z) = (M_mkw_env*sigma_mkw_env*loc_env(x,y,z))
           D_mkw(x,y,z) = D_mkw_met(x,y,z)+D_mkw_pht(x,y,z)+D_mkw_pyr(x,y,z)+D_mkw_env(x,y,z)

           D_pht_met(x,y,z) = (M_pht_met*sigma_pht_met*loc_met(x,y,z))
           D_pht_mkw(x,y,z) = (M_pht_mkw*sigma_pht_mkw*loc_mkw(x,y,z))
           D_pht_pyr(x,y,z) = (M_pht_pyr*sigma_pht_pyr*loc_pyr(x,y,z))
           D_pht_env(x,y,z) = (M_pht_env*sigma_pht_env*loc_env(x,y,z))
           D_pht(x,y,z) = D_pht_met(x,y,z)+D_pht_mkw(x,y,z)+D_pht_pyr(x,y,z)+D_pht_env(x,y,z)

           D_pyr_met(x,y,z) = (M_pyr_met*sigma_pyr_met*loc_met(x,y,z))
           D_pyr_mkw(x,y,z) = (M_pyr_mkw*sigma_pyr_mkw*loc_mkw(x,y,z))
           D_pyr_pht(x,y,z) = (M_pyr_pht*sigma_pyr_pht*loc_pht(x,y,z))
           D_pyr_env(x,y,z) = (M_pyr_env*sigma_pyr_env*loc_env(x,y,z))
           D_pyr(x,y,z) = D_pyr_met(x,y,z)+D_pyr_mkw(x,y,z)+D_pyr_pht(x,y,z)+D_pyr_env(x,y,z)

           D_env_met(x,y,z) = (M_env_met*sigma_env_met*loc_met(x,y,z))
           D_env_mkw(x,y,z) = (M_env_mkw*sigma_env_mkw*loc_mkw(x,y,z))
           D_env_pht(x,y,z) = (M_env_pht*sigma_env_pht*loc_pht(x,y,z))
           D_env_pyr(x,y,z) = (M_env_pyr*sigma_env_pyr*loc_pyr(x,y,z))
           D_env(x,y,z) = D_env_met(x,y,z)+D_env_mkw(x,y,z)+D_env_pht(x,y,z)+D_env_pyr(x,y,z)

        end do
     end do
  end do


  allocate(A(((7*psx*psy*psz)-(2*psx*psy))*no_fields))
  allocate(JA(((7*psx*psy*psz)-(2*psx*psy))*no_fields))
  allocate(IA((no_fields*psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0


!!!! FOR MET
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!!!! MET-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)-1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y-1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x-1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_met(x,y,(z+1)+1)+D_met(x,y,(z+1)-1)+&
                &D_met(x,wrap(y+1,psy),z+1)+D_met(x,wrap(y-1,psy),z+1)+&
                &D_met(wrap(x+1,psx),y,z+1)+D_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x+1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y+1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)+1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if



!!! MET-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,y,(z+1)-1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,wrap(y-1,psy),z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(wrap(x-1,psx),y,z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_mkw_met(x,y,(z+1)+1)+D_mkw_met(x,y,(z+1)-1)+&
                &D_mkw_met(x,wrap(y+1,psy),z+1)+D_mkw_met(x,wrap(y-1,psy),z+1)+&
                &D_mkw_met(wrap(x+1,psx),y,z+1)+D_mkw_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(wrap(x+1,psx),y,z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,wrap(y+1,psy),z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,y,(z+1)+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if



!!! MET-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,y,(z+1)-1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,wrap(y-1,psy),z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(wrap(x-1,psx),y,z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pht_met(x,y,(z+1)+1)+D_pht_met(x,y,(z+1)-1)+&
                &D_pht_met(x,wrap(y+1,psy),z+1)+D_pht_met(x,wrap(y-1,psy),z+1)+&
                &D_pht_met(wrap(x+1,psx),y,z+1)+D_pht_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(wrap(x+1,psx),y,z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,wrap(y+1,psy),z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,y,(z+1)+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if






!!! MET-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,y,(z+1)-1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,wrap(y-1,psy),z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(wrap(x-1,psx),y,z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pyr_met(x,y,(z+1)+1)+D_pyr_met(x,y,(z+1)-1)+&
                &D_pyr_met(x,wrap(y+1,psy),z+1)+D_pyr_met(x,wrap(y-1,psy),z+1)+&
                &D_pyr_met(wrap(x+1,psx),y,z+1)+D_pyr_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(wrap(x+1,psx),y,z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,wrap(y+1,psy),z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,y,(z+1)+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if






!!! MET-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,y,(z+1)-1)+D_env_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,wrap(y-1,psy),z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(wrap(x-1,psx),y,z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_env_met(x,y,(z+1)+1)+D_env_met(x,y,(z+1)-1)+&
                &D_env_met(x,wrap(y+1,psy),z+1)+D_env_met(x,wrap(y-1,psy),z+1)+&
                &D_env_met(wrap(x+1,psx),y,z+1)+D_env_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(wrap(x+1,psx),y,z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,wrap(y+1,psy),z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,y,(z+1)+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

        end do
     end do
  end do
















!!!! FOR MKW
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! MKW-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,y,(z+1)-1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,wrap(y-1,psy),z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(wrap(x-1,psx),y,z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_met_mkw(x,y,(z+1)+1)+D_met_mkw(x,y,(z+1)-1)+&
                &D_met_mkw(x,wrap(y+1,psy),z+1)+D_met_mkw(x,wrap(y-1,psy),z+1)+&
                &D_met_mkw(wrap(x+1,psx),y,z+1)+D_met_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(wrap(x+1,psx),y,z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,wrap(y+1,psy),z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,y,(z+1)+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if



!!!!! MKW-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,y,(z+1)-1)+D_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,wrap(y-1,psy),z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(wrap(x-1,psx),y,z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_mkw(x,y,(z+1)+1)+D_mkw(x,y,(z+1)-1)+&
                &D_mkw(x,wrap(y+1,psy),z+1)+D_mkw(x,wrap(y-1,psy),z+1)+&
                &D_mkw(wrap(x+1,psx),y,z+1)+D_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(wrap(x+1,psx),y,z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,wrap(y+1,psy),z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,y,(z+1)+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if



!!! MKW-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,y,(z+1)-1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,wrap(y-1,psy),z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(wrap(x-1,psx),y,z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pht_mkw(x,y,(z+1)+1)+D_pht_mkw(x,y,(z+1)-1)+&
                &D_pht_mkw(x,wrap(y+1,psy),z+1)+D_pht_mkw(x,wrap(y-1,psy),z+1)+&
                &D_pht_mkw(wrap(x+1,psx),y,z+1)+D_pht_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(wrap(x+1,psx),y,z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,wrap(y+1,psy),z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,y,(z+1)+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if






!!! MKW-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,y,(z+1)-1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,wrap(y-1,psy),z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(wrap(x-1,psx),y,z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pyr_mkw(x,y,(z+1)+1)+D_pyr_mkw(x,y,(z+1)-1)+&
                &D_pyr_mkw(x,wrap(y+1,psy),z+1)+D_pyr_mkw(x,wrap(y-1,psy),z+1)+&
                &D_pyr_mkw(wrap(x+1,psx),y,z+1)+D_pyr_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(wrap(x+1,psx),y,z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,wrap(y+1,psy),z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,y,(z+1)+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if






!!! MKW-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,y,(z+1)-1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,wrap(y-1,psy),z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(wrap(x-1,psx),y,z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_env_mkw(x,y,(z+1)+1)+D_env_mkw(x,y,(z+1)-1)+&
                &D_env_mkw(x,wrap(y+1,psy),z+1)+D_env_mkw(x,wrap(y-1,psy),z+1)+&
                &D_env_mkw(wrap(x+1,psx),y,z+1)+D_env_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(wrap(x+1,psx),y,z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,wrap(y+1,psy),z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,y,(z+1)+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

        end do
     end do
  end do

















!!!! FOR PHT
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PHT-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,y,(z+1)-1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,wrap(y-1,psy),z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(wrap(x-1,psx),y,z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_met_pht(x,y,(z+1)+1)+D_met_pht(x,y,(z+1)-1)+&
                &D_met_pht(x,wrap(y+1,psy),z+1)+D_met_pht(x,wrap(y-1,psy),z+1)+&
                &D_met_pht(wrap(x+1,psx),y,z+1)+D_met_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(wrap(x+1,psx),y,z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,wrap(y+1,psy),z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,y,(z+1)+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if


!!! PHT-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,y,(z+1)-1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,wrap(y-1,psy),z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(wrap(x-1,psx),y,z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_mkw_pht(x,y,(z+1)+1)+D_mkw_pht(x,y,(z+1)-1)+&
                &D_mkw_pht(x,wrap(y+1,psy),z+1)+D_mkw_pht(x,wrap(y-1,psy),z+1)+&
                &D_mkw_pht(wrap(x+1,psx),y,z+1)+D_mkw_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(wrap(x+1,psx),y,z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,wrap(y+1,psy),z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,y,(z+1)+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if


!!!!! PHT-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)-1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y-1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x-1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pht(x,y,(z+1)+1)+D_pht(x,y,(z+1)-1)+&
                &D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,wrap(y-1,psy),z+1)+&
                &D_pht(wrap(x+1,psx),y,z+1)+D_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x+1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)+1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if






!!! PHT-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,y,(z+1)-1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,wrap(y-1,psy),z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(wrap(x-1,psx),y,z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pyr_pht(x,y,(z+1)+1)+D_pyr_pht(x,y,(z+1)-1)+&
                &D_pyr_pht(x,wrap(y+1,psy),z+1)+D_pyr_pht(x,wrap(y-1,psy),z+1)+&
                &D_pyr_pht(wrap(x+1,psx),y,z+1)+D_pyr_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(wrap(x+1,psx),y,z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,wrap(y+1,psy),z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,y,(z+1)+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if






!!! PHT-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,y,(z+1)-1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,wrap(y-1,psy),z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(wrap(x-1,psx),y,z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_env_pht(x,y,(z+1)+1)+D_env_pht(x,y,(z+1)-1)+&
                &D_env_pht(x,wrap(y+1,psy),z+1)+D_env_pht(x,wrap(y-1,psy),z+1)+&
                &D_env_pht(wrap(x+1,psx),y,z+1)+D_env_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(wrap(x+1,psx),y,z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,wrap(y+1,psy),z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,y,(z+1)+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

        end do
     end do
  end do



























!!!! FOR PYR
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PYR-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,y,(z+1)-1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,wrap(y-1,psy),z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(wrap(x-1,psx),y,z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_met_pyr(x,y,(z+1)+1)+D_met_pyr(x,y,(z+1)-1)+&
                &D_met_pyr(x,wrap(y+1,psy),z+1)+D_met_pyr(x,wrap(y-1,psy),z+1)+&
                &D_met_pyr(wrap(x+1,psx),y,z+1)+D_met_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(wrap(x+1,psx),y,z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,wrap(y+1,psy),z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,y,(z+1)+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if


!!! PYR-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,y,(z+1)-1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,wrap(y-1,psy),z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(wrap(x-1,psx),y,z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_mkw_pyr(x,y,(z+1)+1)+D_mkw_pyr(x,y,(z+1)-1)+&
                &D_mkw_pyr(x,wrap(y+1,psy),z+1)+D_mkw_pyr(x,wrap(y-1,psy),z+1)+&
                &D_mkw_pyr(wrap(x+1,psx),y,z+1)+D_mkw_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(wrap(x+1,psx),y,z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,wrap(y+1,psy),z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,y,(z+1)+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if



!!! PYR-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,y,(z+1)-1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,wrap(y-1,psy),z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(wrap(x-1,psx),y,z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pht_pyr(x,y,(z+1)+1)+D_pht_pyr(x,y,(z+1)-1)+&
                &D_pht_pyr(x,wrap(y+1,psy),z+1)+D_pht_pyr(x,wrap(y-1,psy),z+1)+&
                &D_pht_pyr(wrap(x+1,psx),y,z+1)+D_pht_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(wrap(x+1,psx),y,z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,wrap(y+1,psy),z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,y,(z+1)+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if



!!!!! PYR-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)-1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y-1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x-1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pyr(x,y,(z+1)+1)+D_pyr(x,y,(z+1)-1)+&
                &D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,wrap(y-1,psy),z+1)+&
                &D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if


!!! PYR-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,y,(z+1)-1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,wrap(y-1,psy),z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(wrap(x-1,psx),y,z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_env_pyr(x,y,(z+1)+1)+D_env_pyr(x,y,(z+1)-1)+&
                &D_env_pyr(x,wrap(y+1,psy),z+1)+D_env_pyr(x,wrap(y-1,psy),z+1)+&
                &D_env_pyr(wrap(x+1,psx),y,z+1)+D_env_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(wrap(x+1,psx),y,z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,wrap(y+1,psy),z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,y,(z+1)+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

        end do
     end do
  end do




























!!!! FOR ENV
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! ENV-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,y,(z+1)-1)+D_met_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,wrap(y-1,psy),z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(wrap(x-1,psx),y,z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_met_env(x,y,(z+1)+1)+D_met_env(x,y,(z+1)-1)+&
                &D_met_env(x,wrap(y+1,psy),z+1)+D_met_env(x,wrap(y-1,psy),z+1)+&
                &D_met_env(wrap(x+1,psx),y,z+1)+D_met_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(wrap(x+1,psx),y,z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,wrap(y+1,psy),z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,y,(z+1)+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           end if


!!! ENV-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,y,(z+1)-1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,wrap(y-1,psy),z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(wrap(x-1,psx),y,z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_mkw_env(x,y,(z+1)+1)+D_mkw_env(x,y,(z+1)-1)+&
                &D_mkw_env(x,wrap(y+1,psy),z+1)+D_mkw_env(x,wrap(y-1,psy),z+1)+&
                &D_mkw_env(wrap(x+1,psx),y,z+1)+D_mkw_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(wrap(x+1,psx),y,z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,wrap(y+1,psy),z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,y,(z+1)+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           end if



!!! ENV-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,y,(z+1)-1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,wrap(y-1,psy),z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(wrap(x-1,psx),y,z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pht_env(x,y,(z+1)+1)+D_pht_env(x,y,(z+1)-1)+&
                &D_pht_env(x,wrap(y+1,psy),z+1)+D_pht_env(x,wrap(y-1,psy),z+1)+&
                &D_pht_env(wrap(x+1,psx),y,z+1)+D_pht_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(wrap(x+1,psx),y,z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,wrap(y+1,psy),z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,y,(z+1)+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           end if





!!! ENV-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,y,(z+1)-1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,wrap(y-1,psy),z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(wrap(x-1,psx),y,z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_pyr_env(x,y,(z+1)+1)+D_pyr_env(x,y,(z+1)-1)+&
                &D_pyr_env(x,wrap(y+1,psy),z+1)+D_pyr_env(x,wrap(y-1,psy),z+1)+&
                &D_pyr_env(wrap(x+1,psx),y,z+1)+D_pyr_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(wrap(x+1,psx),y,z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,wrap(y+1,psy),z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,y,(z+1)+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           end if




!!!!! ENV-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)-1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y-1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x-1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = D_env(x,y,(z+1)+1)+D_env(x,y,(z+1)-1)+&
                &D_env(x,wrap(y+1,psy),z+1)+D_env(x,wrap(y-1,psy),z+1)+&
                &D_env(wrap(x+1,psx),y,z+1)+D_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x+1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*psz)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y+1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*psz)

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)+1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           end if

        end do
     end do
  end do
































































  IA((psx*psy*psz*no_fields)+1)= IA(psx*psy*psz*no_fields)+6









end subroutine pfFunction
