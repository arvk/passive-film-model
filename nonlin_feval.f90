subroutine pfFunction(snes,pf_vec,ret_vec,dummy,ierr)
  use commondata
  use fields
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

  integer, parameter :: imet = 1
  integer, parameter :: imkw = 2
  integer, parameter :: ipht = 3
  integer, parameter :: ipyr = 4
  integer, parameter :: ienv = 5

  integer, parameter :: no_fields = 5

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: vector_locator
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)*no_fields) :: const,lnr,sqr

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
  real*8 :: f_met, f_mkw, f_pht, f_pyr, f_env
  real*8 :: w_met, w_mkw, w_pht, w_pyr, w_env

  real*8 :: sigma_pyr_met, sigma_pyr_mkw, sigma_pyr_pht, sigma_pyr_env
  real*8 :: sigma_met_pyr, sigma_mkw_pyr, sigma_pht_pyr, sigma_env_pyr

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_met, D_met_mkw, D_met_pht, D_met_pyr, D_met_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_mkw_met, D_mkw, D_mkw_pht, D_mkw_pyr, D_mkw_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_pht_met, D_pht_mkw, D_pht, D_pht_pyr, D_pht_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_pyr_met, D_pyr_mkw, D_pyr_pht, D_pyr, D_pyr_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D_env_met, D_env_mkw, D_env_pht, D_env_pyr, D_env

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: loc_met, loc_mkw, loc_pht, loc_pyr, loc_env

  real*8 :: S = 0.0d0!2E1

  real*8 :: odiff
  real*8 :: delx,dely,delz
  real*8 :: del_omet,del_omkw,del_opht,del_opyr,del_oenv

  call calc_grad_pf()

  call VecGetArrayF90(pf_vec,point_pf_vec,ierr)
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           loc_met(x,y,z+1) = point_pf_vec(linindex+((imet-1)*psx*psy*(psz+(2*ghost_width)-2)))
           loc_mkw(x,y,z+1) = point_pf_vec(linindex+((imkw-1)*psx*psy*(psz+(2*ghost_width)-2)))
           loc_pht(x,y,z+1) = point_pf_vec(linindex+((ipht-1)*psx*psy*(psz+(2*ghost_width)-2)))
           loc_pyr(x,y,z+1) = point_pf_vec(linindex+((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2)))
           loc_env(x,y,z+1) = point_pf_vec(linindex+((ienv-1)*psx*psy*(psz+(2*ghost_width)-2)))

        end do
     end do
  end do
  call VecRestoreArrayF90(pf_vec,point_pf_vec,ierr)


  do z = 1,ghost_width
     do y = 1,psy
        do x = 1,psx
           loc_met(x,y,z) = loc_met(x,y,1+ghost_width) ; loc_met(x,y,psz+z+ghost_width) = loc_met(x,y,psz+ghost_width) 
           loc_mkw(x,y,z) = loc_mkw(x,y,1+ghost_width) ; loc_mkw(x,y,psz+z+ghost_width) = loc_mkw(x,y,psz+ghost_width) 
           loc_pht(x,y,z) = loc_pht(x,y,1+ghost_width) ; loc_pht(x,y,psz+z+ghost_width) = loc_pht(x,y,psz+ghost_width) 
           loc_pyr(x,y,z) = loc_pyr(x,y,1+ghost_width) ; loc_pyr(x,y,psz+z+ghost_width) = loc_pyr(x,y,psz+ghost_width) 
           loc_env(x,y,z) = loc_env(x,y,1+ghost_width) ; loc_env(x,y,psz+z+ghost_width) = loc_env(x,y,psz+ghost_width) 
        end do
     end do
  end do


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

           D_met_mkw(x,y,z) = 0.0d0 - (M_met_mkw*sigma_met_mkw*abs(loc_mkw(x,y,z)))
           D_met_pht(x,y,z) = 0.0d0 - (M_met_pht*sigma_met_pht*abs(loc_pht(x,y,z)))
           D_met_pyr(x,y,z) = 0.0d0 - (M_met_pyr*sigma_met_pyr*abs(loc_pyr(x,y,z)))
           D_met_env(x,y,z) = 0.0d0 - (M_met_env*sigma_met_env*abs(loc_env(x,y,z)))
           D_met(x,y,z) = 0.0d0 - (D_met_mkw(x,y,z)+D_met_pht(x,y,z)+D_met_pyr(x,y,z)+D_met_env(x,y,z))

           D_mkw_met(x,y,z) = 0.0d0 - (M_mkw_met*sigma_mkw_met*abs(loc_met(x,y,z)))
           D_mkw_pht(x,y,z) = 0.0d0 - (M_mkw_pht*sigma_mkw_pht*abs(loc_pht(x,y,z)))
           D_mkw_pyr(x,y,z) = 0.0d0 - (M_mkw_pyr*sigma_mkw_pyr*abs(loc_pyr(x,y,z)))
           D_mkw_env(x,y,z) = 0.0d0 - (M_mkw_env*sigma_mkw_env*abs(loc_env(x,y,z)))
           D_mkw(x,y,z) = 0.0d0 - (D_mkw_met(x,y,z)+D_mkw_pht(x,y,z)+D_mkw_pyr(x,y,z)+D_mkw_env(x,y,z))

           D_pht_met(x,y,z) = 0.0d0 - (M_pht_met*sigma_pht_met*abs(loc_met(x,y,z)))
           D_pht_mkw(x,y,z) = 0.0d0 - (M_pht_mkw*sigma_pht_mkw*abs(loc_mkw(x,y,z)))
           D_pht_pyr(x,y,z) = 0.0d0 - (M_pht_pyr*sigma_pht_pyr*abs(loc_pyr(x,y,z)))
           D_pht_env(x,y,z) = 0.0d0 - (M_pht_env*sigma_pht_env*abs(loc_env(x,y,z)))
           D_pht(x,y,z) = 0.0d0 - (D_pht_met(x,y,z)+D_pht_mkw(x,y,z)+D_pht_pyr(x,y,z)+D_pht_env(x,y,z))

           D_pyr_met(x,y,z) = 0.0d0 - (M_pyr_met*sigma_pyr_met*abs(loc_met(x,y,z)))
           D_pyr_mkw(x,y,z) = 0.0d0 - (M_pyr_mkw*sigma_pyr_mkw*abs(loc_mkw(x,y,z)))
           D_pyr_pht(x,y,z) = 0.0d0 - (M_pyr_pht*sigma_pyr_pht*abs(loc_pht(x,y,z)))
           D_pyr_env(x,y,z) = 0.0d0 - (M_pyr_env*sigma_pyr_env*abs(loc_env(x,y,z)))
           D_pyr(x,y,z) = 0.0d0 - (D_pyr_met(x,y,z)+D_pyr_mkw(x,y,z)+D_pyr_pht(x,y,z)+D_pyr_env(x,y,z))

           D_env_met(x,y,z) = 0.0d0 - (M_env_met*sigma_env_met*abs(loc_met(x,y,z)))
           D_env_mkw(x,y,z) = 0.0d0 - (M_env_mkw*sigma_env_mkw*abs(loc_mkw(x,y,z)))
           D_env_pht(x,y,z) = 0.0d0 - (M_env_pht*sigma_env_pht*abs(loc_pht(x,y,z)))
           D_env_pyr(x,y,z) = 0.0d0 - (M_env_pyr*sigma_env_pyr*abs(loc_pyr(x,y,z)))
           D_env(x,y,z) = 0.0d0 - (D_env_met(x,y,z)+D_env_mkw(x,y,z)+D_env_pht(x,y,z)+D_env_pyr(x,y,z))

        end do
     end do
  end do


  allocate(A(((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy))*no_fields*no_fields))
  allocate(JA(((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy))*no_fields*no_fields))
  allocate(IA((no_fields*psx*psy*(psz+(2*ghost_width)-2))+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

!!!! FOR MET
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0






           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!!!! MET-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)-1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y-1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x-1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_met(x,y,(z+1)+1)+D_met(x,y,(z+1)-1)+&
                &D_met(x,wrap(y+1,psy),z+1)+D_met(x,wrap(y-1,psy),z+1)+&
                &D_met(wrap(x+1,psx),y,z+1)+D_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x+1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y+1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)+1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!! MET-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,y,(z+1)-1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,wrap(y-1,psy),z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(wrap(x-1,psx),y,z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_mkw_met(x,y,(z+1)+1)+D_mkw_met(x,y,(z+1)-1)+&
                &D_mkw_met(x,wrap(y+1,psy),z+1)+D_mkw_met(x,wrap(y-1,psy),z+1)+&
                &D_mkw_met(wrap(x+1,psx),y,z+1)+D_mkw_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(wrap(x+1,psx),y,z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,wrap(y+1,psy),z+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_met(x,y,(z+1)+1)+D_mkw_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!! MET-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,y,(z+1)-1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,wrap(y-1,psy),z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(wrap(x-1,psx),y,z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pht_met(x,y,(z+1)+1)+D_pht_met(x,y,(z+1)-1)+&
                &D_pht_met(x,wrap(y+1,psy),z+1)+D_pht_met(x,wrap(y-1,psy),z+1)+&
                &D_pht_met(wrap(x+1,psx),y,z+1)+D_pht_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(wrap(x+1,psx),y,z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,wrap(y+1,psy),z+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_met(x,y,(z+1)+1)+D_pht_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! MET-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,y,(z+1)-1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,wrap(y-1,psy),z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(wrap(x-1,psx),y,z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pyr_met(x,y,(z+1)+1)+D_pyr_met(x,y,(z+1)-1)+&
                &D_pyr_met(x,wrap(y+1,psy),z+1)+D_pyr_met(x,wrap(y-1,psy),z+1)+&
                &D_pyr_met(wrap(x+1,psx),y,z+1)+D_pyr_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(wrap(x+1,psx),y,z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,wrap(y+1,psy),z+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_met(x,y,(z+1)+1)+D_pyr_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! MET-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,y,(z+1)-1)+D_env_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,wrap(y-1,psy),z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(wrap(x-1,psx),y,z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_env_met(x,y,(z+1)+1)+D_env_met(x,y,(z+1)-1)+&
                &D_env_met(x,wrap(y+1,psy),z+1)+D_env_met(x,wrap(y-1,psy),z+1)+&
                &D_env_met(wrap(x+1,psx),y,z+1)+D_env_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(wrap(x+1,psx),y,z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,wrap(y+1,psy),z+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_met(x,y,(z+1)+1)+D_env_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

        end do
     end do
  end do








!!!! FOR MKW
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx



           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0





           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! MKW-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,y,(z+1)-1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,wrap(y-1,psy),z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(wrap(x-1,psx),y,z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_met_mkw(x,y,(z+1)+1)+D_met_mkw(x,y,(z+1)-1)+&
                &D_met_mkw(x,wrap(y+1,psy),z+1)+D_met_mkw(x,wrap(y-1,psy),z+1)+&
                &D_met_mkw(wrap(x+1,psx),y,z+1)+D_met_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(wrap(x+1,psx),y,z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,wrap(y+1,psy),z+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_mkw(x,y,(z+1)+1)+D_met_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!!!! MKW-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,y,(z+1)-1)+D_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,wrap(y-1,psy),z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(wrap(x-1,psx),y,z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_mkw(x,y,(z+1)+1)+D_mkw(x,y,(z+1)-1)+&
                &D_mkw(x,wrap(y+1,psy),z+1)+D_mkw(x,wrap(y-1,psy),z+1)+&
                &D_mkw(wrap(x+1,psx),y,z+1)+D_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(wrap(x+1,psx),y,z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,wrap(y+1,psy),z+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw(x,y,(z+1)+1)+D_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!! MKW-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,y,(z+1)-1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,wrap(y-1,psy),z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(wrap(x-1,psx),y,z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pht_mkw(x,y,(z+1)+1)+D_pht_mkw(x,y,(z+1)-1)+&
                &D_pht_mkw(x,wrap(y+1,psy),z+1)+D_pht_mkw(x,wrap(y-1,psy),z+1)+&
                &D_pht_mkw(wrap(x+1,psx),y,z+1)+D_pht_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(wrap(x+1,psx),y,z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,wrap(y+1,psy),z+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_mkw(x,y,(z+1)+1)+D_pht_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! MKW-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,y,(z+1)-1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,wrap(y-1,psy),z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(wrap(x-1,psx),y,z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pyr_mkw(x,y,(z+1)+1)+D_pyr_mkw(x,y,(z+1)-1)+&
                &D_pyr_mkw(x,wrap(y+1,psy),z+1)+D_pyr_mkw(x,wrap(y-1,psy),z+1)+&
                &D_pyr_mkw(wrap(x+1,psx),y,z+1)+D_pyr_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(wrap(x+1,psx),y,z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,wrap(y+1,psy),z+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_mkw(x,y,(z+1)+1)+D_pyr_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! MKW-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,y,(z+1)-1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,wrap(y-1,psy),z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(wrap(x-1,psx),y,z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_env_mkw(x,y,(z+1)+1)+D_env_mkw(x,y,(z+1)-1)+&
                &D_env_mkw(x,wrap(y+1,psy),z+1)+D_env_mkw(x,wrap(y-1,psy),z+1)+&
                &D_env_mkw(wrap(x+1,psx),y,z+1)+D_env_mkw(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_mkw(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(wrap(x+1,psx),y,z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,wrap(y+1,psy),z+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_mkw(x,y,(z+1)+1)+D_env_mkw(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

        end do
     end do
  end do












!!!! FOR PHT
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx




           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0





           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PHT-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,y,(z+1)-1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,wrap(y-1,psy),z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(wrap(x-1,psx),y,z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_met_pht(x,y,(z+1)+1)+D_met_pht(x,y,(z+1)-1)+&
                &D_met_pht(x,wrap(y+1,psy),z+1)+D_met_pht(x,wrap(y-1,psy),z+1)+&
                &D_met_pht(wrap(x+1,psx),y,z+1)+D_met_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(wrap(x+1,psx),y,z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,wrap(y+1,psy),z+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pht(x,y,(z+1)+1)+D_met_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if


!!! PHT-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,y,(z+1)-1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,wrap(y-1,psy),z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(wrap(x-1,psx),y,z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_mkw_pht(x,y,(z+1)+1)+D_mkw_pht(x,y,(z+1)-1)+&
                &D_mkw_pht(x,wrap(y+1,psy),z+1)+D_mkw_pht(x,wrap(y-1,psy),z+1)+&
                &D_mkw_pht(wrap(x+1,psx),y,z+1)+D_mkw_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(wrap(x+1,psx),y,z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,wrap(y+1,psy),z+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pht(x,y,(z+1)+1)+D_mkw_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if


!!!!! PHT-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)-1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y-1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x-1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pht(x,y,(z+1)+1)+D_pht(x,y,(z+1)-1)+&
                &D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,wrap(y-1,psy),z+1)+&
                &D_pht(wrap(x+1,psx),y,z+1)+D_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x+1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)+1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! PHT-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,y,(z+1)-1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,wrap(y-1,psy),z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(wrap(x-1,psx),y,z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pyr_pht(x,y,(z+1)+1)+D_pyr_pht(x,y,(z+1)-1)+&
                &D_pyr_pht(x,wrap(y+1,psy),z+1)+D_pyr_pht(x,wrap(y-1,psy),z+1)+&
                &D_pyr_pht(wrap(x+1,psx),y,z+1)+D_pyr_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(wrap(x+1,psx),y,z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,wrap(y+1,psy),z+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_pht(x,y,(z+1)+1)+D_pyr_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if






!!! PHT-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,y,(z+1)-1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,wrap(y-1,psy),z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(wrap(x-1,psx),y,z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_env_pht(x,y,(z+1)+1)+D_env_pht(x,y,(z+1)-1)+&
                &D_env_pht(x,wrap(y+1,psy),z+1)+D_env_pht(x,wrap(y-1,psy),z+1)+&
                &D_env_pht(wrap(x+1,psx),y,z+1)+D_env_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(wrap(x+1,psx),y,z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,wrap(y+1,psy),z+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pht(x,y,(z+1)+1)+D_env_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

        end do
     end do
  end do

























!!!! FOR PYR
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0



           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PYR-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,y,(z+1)-1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,wrap(y-1,psy),z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(wrap(x-1,psx),y,z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_met_pyr(x,y,(z+1)+1)+D_met_pyr(x,y,(z+1)-1)+&
                &D_met_pyr(x,wrap(y+1,psy),z+1)+D_met_pyr(x,wrap(y-1,psy),z+1)+&
                &D_met_pyr(wrap(x+1,psx),y,z+1)+D_met_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(wrap(x+1,psx),y,z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,wrap(y+1,psy),z+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_pyr(x,y,(z+1)+1)+D_met_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if


!!! PYR-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,y,(z+1)-1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,wrap(y-1,psy),z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(wrap(x-1,psx),y,z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_mkw_pyr(x,y,(z+1)+1)+D_mkw_pyr(x,y,(z+1)-1)+&
                &D_mkw_pyr(x,wrap(y+1,psy),z+1)+D_mkw_pyr(x,wrap(y-1,psy),z+1)+&
                &D_mkw_pyr(wrap(x+1,psx),y,z+1)+D_mkw_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(wrap(x+1,psx),y,z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,wrap(y+1,psy),z+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_pyr(x,y,(z+1)+1)+D_mkw_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!! PYR-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,y,(z+1)-1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,wrap(y-1,psy),z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(wrap(x-1,psx),y,z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pht_pyr(x,y,(z+1)+1)+D_pht_pyr(x,y,(z+1)-1)+&
                &D_pht_pyr(x,wrap(y+1,psy),z+1)+D_pht_pyr(x,wrap(y-1,psy),z+1)+&
                &D_pht_pyr(wrap(x+1,psx),y,z+1)+D_pht_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(wrap(x+1,psx),y,z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,wrap(y+1,psy),z+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_pyr(x,y,(z+1)+1)+D_pht_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!!!! PYR-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)-1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y-1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x-1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pyr(x,y,(z+1)+1)+D_pyr(x,y,(z+1)-1)+&
                &D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,wrap(y-1,psy),z+1)+&
                &D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if


!!! PYR-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,y,(z+1)-1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,wrap(y-1,psy),z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(wrap(x-1,psx),y,z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_env_pyr(x,y,(z+1)+1)+D_env_pyr(x,y,(z+1)-1)+&
                &D_env_pyr(x,wrap(y+1,psy),z+1)+D_env_pyr(x,wrap(y-1,psy),z+1)+&
                &D_env_pyr(wrap(x+1,psx),y,z+1)+D_env_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(wrap(x+1,psx),y,z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,wrap(y+1,psy),z+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env_pyr(x,y,(z+1)+1)+D_env_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

        end do
     end do
  end do
























!!!! FOR ENV
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0




           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! ENV-MET

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,y,(z+1)-1)+D_met_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,wrap(y-1,psy),z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(wrap(x-1,psx),y,z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_met_env(x,y,(z+1)+1)+D_met_env(x,y,(z+1)-1)+&
                &D_met_env(x,wrap(y+1,psy),z+1)+D_met_env(x,wrap(y-1,psy),z+1)+&
                &D_met_env(wrap(x+1,psx),y,z+1)+D_met_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(wrap(x+1,psx),y,z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,wrap(y+1,psy),z+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met_env(x,y,(z+1)+1)+D_met_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if


!!! ENV-MKW

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,y,(z+1)-1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,wrap(y-1,psy),z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(wrap(x-1,psx),y,z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_mkw_env(x,y,(z+1)+1)+D_mkw_env(x,y,(z+1)-1)+&
                &D_mkw_env(x,wrap(y+1,psy),z+1)+D_mkw_env(x,wrap(y-1,psy),z+1)+&
                &D_mkw_env(wrap(x+1,psx),y,z+1)+D_mkw_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_mkw_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(wrap(x+1,psx),y,z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,wrap(y+1,psy),z+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_mkw_env(x,y,(z+1)+1)+D_mkw_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if



!!! ENV-PHT

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,y,(z+1)-1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,wrap(y-1,psy),z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(wrap(x-1,psx),y,z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pht_env(x,y,(z+1)+1)+D_pht_env(x,y,(z+1)-1)+&
                &D_pht_env(x,wrap(y+1,psy),z+1)+D_pht_env(x,wrap(y-1,psy),z+1)+&
                &D_pht_env(wrap(x+1,psx),y,z+1)+D_pht_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(wrap(x+1,psx),y,z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,wrap(y+1,psy),z+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht_env(x,y,(z+1)+1)+D_pht_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if





!!! ENV-PYR

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,y,(z+1)-1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,wrap(y-1,psy),z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(wrap(x-1,psx),y,z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_pyr_env(x,y,(z+1)+1)+D_pyr_env(x,y,(z+1)-1)+&
                &D_pyr_env(x,wrap(y+1,psy),z+1)+D_pyr_env(x,wrap(y-1,psy),z+1)+&
                &D_pyr_env(wrap(x+1,psx),y,z+1)+D_pyr_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(wrap(x+1,psx),y,z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,wrap(y+1,psy),z+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr_env(x,y,(z+1)+1)+D_pyr_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if




!!!!! ENV-ENV

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)-1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y-1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x-1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = D_env(x,y,(z+1)+1)+D_env(x,y,(z+1)-1)+&
                &D_env(x,wrap(y+1,psy),z+1)+D_env(x,wrap(y-1,psy),z+1)+&
                &D_env(wrap(x+1,psx),y,z+1)+D_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x+1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx) + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y+1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)+1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))
           end if

        end do
     end do
  end do

  IA((psx*psy*(psz+(2*ghost_width)-2)*no_fields)+1)= IA(psx*psy*(psz+(2*ghost_width)-2)*no_fields)+(6*no_fields)


  do linindex = 1,psx*psy*(psz+(2*ghost_width)-2)*no_fields
     is_sorted = .FALSE.

     do while (is_sorted .eqv. .FALSE.)

        do rowindex = IA(linindex),IA(linindex+1)-2

           is_sorted = .TRUE.

           JAleft = JA(rowindex)
           JAright = JA(rowindex+1)

           if (JAleft .gt. JAright) then

              JAswap = JA(rowindex)
              JA(rowindex) = JA(rowindex+1)
              JA(rowindex+1) = JAswap

              Aswap = A(rowindex)
              A(rowindex) = A(rowindex+1)
              A(rowindex+1) = Aswap

              is_sorted = .FALSE.
              exit
           end if

        end do


     end do

  end do

  IA = IA - 1
  JA = JA - 1

  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*(psz+(2*ghost_width)-2)*no_fields,psx*psy*(psz+(2*ghost_width)-2)*no_fields,IA,JA,A,lhs_mat,ierr)
  call MatMult(lhs_mat,pf_vec,ret_vec,ierr)

  lnr = 0.0d0
  sqr = 0.0d0

  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0
           w_met = f_met - (mu(x,y,z+1)*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)
           w_mkw = f_mkw - (mu(x,y,z+1)*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)
           w_pht = f_pht - (mu(x,y,z+1))

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)
           w_pyr = f_pyr - (mu(x,y,z+1)*2)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T)))
           w_env = 0.0d0


           del_omet = 0.0d0
           del_omkw = 0.0d0
           del_opht = 0.0d0
           del_oenv = 0.0d0

           delx = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(wrap(x-1,psx),y,z+1))
           dely = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,wrap(y-1,psy),z+1))
           delz = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1-1))
           del_opyr = (sqrt((delx*delx)+(dely*dely)+(delz*delz))/dpf)


           hill_met_mkw = (16.0d0/3.0d0)*double_well_barrier; hill_mkw_met = (16.0d0/3.0d0)*double_well_barrier
           hill_met_pht = (16.0d0/3.0d0)*double_well_barrier; hill_pht_met = (16.0d0/3.0d0)*double_well_barrier
           hill_met_pyr = (16.0d0/3.0d0)*double_well_barrier; hill_pyr_met = (16.0d0/3.0d0)*double_well_barrier
           hill_met_env = (16.0d0/3.0d0)*double_well_barrier; hill_env_met = (16.0d0/3.0d0)*double_well_barrier

           hill_mkw_pht = (16.0d0/3.0d0)*double_well_barrier; hill_pht_mkw = (16.0d0/3.0d0)*double_well_barrier
           hill_mkw_pyr = (16.0d0/3.0d0)*double_well_barrier; hill_pyr_mkw = (16.0d0/3.0d0)*double_well_barrier
           hill_mkw_env = (16.0d0/3.0d0)*double_well_barrier; hill_env_mkw = (16.0d0/3.0d0)*double_well_barrier

           hill_pht_pyr = (16.0d0/3.0d0)*double_well_barrier; hill_pyr_pht = (16.0d0/3.0d0)*double_well_barrier
           hill_pht_env = (16.0d0/3.0d0)*double_well_barrier; hill_env_pht = (16.0d0/3.0d0)*double_well_barrier

           hill_pyr_env = (16.0d0/3.0d0)*double_well_barrier; hill_env_pyr = (16.0d0/3.0d0)*double_well_barrier


!!! FOR MET           
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*(psz+(2*ghost_width)-2))

           const(linindex) = 0.0d0

           lnr(linindex) = (2*M_met_mkw*hill_met_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (S*(del_omet-del_omkw)*loc_mkw(x,y,z+1)*M_met_mkw) + (6*(w_met-w_mkw)*loc_mkw(x,y,z+1)*M_met_mkw) + &
                         & (2*M_met_pht*hill_met_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (S*(del_omet-del_opht)*loc_pht(x,y,z+1)*M_met_pht) + (6*(w_met-w_pht)*loc_pht(x,y,z+1)*M_met_pht) + &
                         & (2*M_met_pyr*hill_met_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (S*(del_omet-del_opyr)*loc_pyr(x,y,z+1)*M_met_pyr) + (6*(w_met-w_pyr)*loc_pyr(x,y,z+1)*M_met_pyr) + &
                         & (2*M_met_env*hill_met_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (S*(del_omet-del_oenv)*loc_env(x,y,z+1)*M_met_env) + (6*(w_met-w_env)*loc_env(x,y,z+1)*M_met_env) 

           lnr(linindex) = lnr(linindex)*loc_met(x,y,z+1)

           sqr(linindex) = 0.0d0 - (M_met_mkw*2*hill_met_mkw*loc_mkw(x,y,z+1)) &
                         &       - (M_met_pht*2*hill_met_pht*loc_pht(x,y,z+1)) &
                         &       - (M_met_pyr*2*hill_met_pyr*loc_pyr(x,y,z+1)) &
                         &       - (M_met_env*2*hill_met_env*loc_env(x,y,z+1)) 

           sqr(linindex) = sqr(linindex)*loc_met(x,y,z+1)*loc_met(x,y,z+1)


!!! FOR MKW           
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*(psz+(2*ghost_width)-2))

           const(linindex) = 0.0d0

           lnr(linindex) = (2*M_mkw_met*hill_mkw_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (S*(del_omkw-del_omet)*loc_met(x,y,z+1)*M_mkw_met) + (6*(w_mkw-w_met)*loc_met(x,y,z+1)*M_mkw_met) + &
                         & (2*M_mkw_pht*hill_mkw_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (S*(del_omkw-del_opht)*loc_pht(x,y,z+1)*M_mkw_pht) + (6*(w_mkw-w_pht)*loc_pht(x,y,z+1)*M_mkw_pht) + &
                         & (2*M_mkw_pyr*hill_mkw_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (S*(del_omkw-del_opyr)*loc_pyr(x,y,z+1)*M_mkw_pyr) + (6*(w_mkw-w_pyr)*loc_pyr(x,y,z+1)*M_mkw_pyr) + &
                         & (2*M_mkw_env*hill_mkw_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (S*(del_omkw-del_oenv)*loc_env(x,y,z+1)*M_mkw_env) + (6*(w_mkw-w_env)*loc_env(x,y,z+1)*M_mkw_env) 

           lnr(linindex) = lnr(linindex)*loc_mkw(x,y,z+1)

           sqr(linindex) = 0.0d0 - (M_mkw_met*2*hill_mkw_met*loc_met(x,y,z+1)) &
                         &       - (M_mkw_pht*2*hill_mkw_pht*loc_pht(x,y,z+1)) &
                         &       - (M_mkw_pyr*2*hill_mkw_pyr*loc_pyr(x,y,z+1)) &
                         &       - (M_mkw_env*2*hill_mkw_env*loc_env(x,y,z+1)) 

           sqr(linindex) = sqr(linindex)*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)


!!! FOR PHT           
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*(psz+(2*ghost_width)-2))

           const(linindex) = 0.0d0

           lnr(linindex) = (2*M_pht_met*hill_pht_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (S*(del_opht-del_omet)*loc_met(x,y,z+1)*M_pht_met) + (6*(w_pht-w_met)*loc_met(x,y,z+1)*M_pht_met) + &
                         & (2*M_pht_mkw*hill_pht_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (S*(del_opht-del_omkw)*loc_mkw(x,y,z+1)*M_pht_mkw) + (6*(w_pht-w_mkw)*loc_mkw(x,y,z+1)*M_pht_mkw) + &
                         & (2*M_pht_pyr*hill_pht_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (S*(del_opht-del_opyr)*loc_pyr(x,y,z+1)*M_pht_pyr) + (6*(w_pht-w_pyr)*loc_pyr(x,y,z+1)*M_pht_pyr) + &
                         & (2*M_pht_env*hill_pht_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (S*(del_opht-del_oenv)*loc_env(x,y,z+1)*M_pht_env) + (6*(w_pht-w_env)*loc_env(x,y,z+1)*M_pht_env) 

           lnr(linindex) = lnr(linindex)*loc_pht(x,y,z+1)

           sqr(linindex) = 0.0d0 - (M_pht_met*2*hill_pht_met*loc_met(x,y,z+1)) &
                         &       - (M_pht_mkw*2*hill_pht_mkw*loc_mkw(x,y,z+1)) &
                         &       - (M_pht_pyr*2*hill_pht_pyr*loc_pyr(x,y,z+1)) &
                         &       - (M_pht_env*2*hill_pht_env*loc_env(x,y,z+1)) 

           sqr(linindex) = sqr(linindex)*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)




!!! FOR PYR           
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*(psz+(2*ghost_width)-2))

           const(linindex) = 0.0d0

           lnr(linindex) = (2*M_pyr_met*hill_pyr_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (S*(del_opyr-del_omet)*loc_met(x,y,z+1)*M_pyr_met) + (6*(w_pyr-w_met)*loc_met(x,y,z+1)*M_pyr_met) + &
                         & (2*M_pyr_mkw*hill_pyr_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (S*(del_opyr-del_omkw)*loc_mkw(x,y,z+1)*M_pyr_mkw) + (6*(w_pyr-w_mkw)*loc_mkw(x,y,z+1)*M_pyr_mkw) + &
                         & (2*M_pyr_pht*hill_pyr_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (S*(del_opyr-del_opht)*loc_pht(x,y,z+1)*M_pyr_pht) + (6*(w_pyr-w_pht)*loc_pht(x,y,z+1)*M_pyr_pht) + &
                         & (2*M_pyr_env*hill_pyr_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (S*(del_opyr-del_oenv)*loc_env(x,y,z+1)*M_pyr_env) + (6*(w_pyr-w_env)*loc_env(x,y,z+1)*M_pyr_env) 

           lnr(linindex) = lnr(linindex)*loc_pyr(x,y,z+1)

           sqr(linindex) = 0.0d0 - (M_pyr_met*2*hill_pyr_met*loc_met(x,y,z+1)) &
                         &       - (M_pyr_mkw*2*hill_pyr_mkw*loc_mkw(x,y,z+1)) &
                         &       - (M_pyr_pht*2*hill_pyr_pht*loc_pht(x,y,z+1)) &
                         &       - (M_pyr_env*2*hill_pyr_env*loc_env(x,y,z+1)) 

           sqr(linindex) = sqr(linindex)*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)



!!! FOR ENV           
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*(psz+(2*ghost_width)-2))

           const(linindex) = 0.0d0

           lnr(linindex) = (2*M_env_met*hill_env_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (S*(del_oenv-del_omet)*loc_met(x,y,z+1)*M_env_met) + (6*(w_env-w_met)*loc_met(x,y,z+1)*M_env_met) + &
                         & (2*M_env_mkw*hill_env_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (S*(del_oenv-del_omkw)*loc_mkw(x,y,z+1)*M_env_mkw) + (6*(w_env-w_mkw)*loc_mkw(x,y,z+1)*M_env_mkw) + &
                         & (2*M_env_pht*hill_env_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (S*(del_oenv-del_opht)*loc_pht(x,y,z+1)*M_env_pht) + (6*(w_env-w_pht)*loc_pht(x,y,z+1)*M_env_pht) + &
                         & (2*M_env_pyr*hill_env_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (S*(del_oenv-del_opyr)*loc_pyr(x,y,z+1)*M_env_pyr) + (6*(w_env-w_pyr)*loc_pyr(x,y,z+1)*M_env_pyr)

           lnr(linindex) = lnr(linindex)*loc_env(x,y,z+1)

           sqr(linindex) = 0.0d0 - (M_env_met*2*hill_env_met*loc_met(x,y,z+1)) &
                         &       - (M_env_mkw*2*hill_env_mkw*loc_mkw(x,y,z+1)) &
                         &       - (M_env_pht*2*hill_env_pht*loc_pht(x,y,z+1)) &
                         &       - (M_env_pyr*2*hill_env_pyr*loc_pyr(x,y,z+1))

           sqr(linindex) = sqr(linindex)*loc_env(x,y,z+1)*loc_env(x,y,z+1)


        end do
     end do
  end do
  

  call VecSetValues(ret_vec,psx*psy*(psz+(2*ghost_width)-2)*no_fields,vector_locator,const,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*(psz+(2*ghost_width)-2)*no_fields,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*(psz+(2*ghost_width)-2)*no_fields,vector_locator,sqr,ADD_VALUES,ierr)

  call VecAssemblyBegin(ret_vec,ierr)
  call VecAssemblyEnd(ret_vec,ierr)

  call MatDestroy(lhs_mat,ierr)

  return

end subroutine pfFunction
