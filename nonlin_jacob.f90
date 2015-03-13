subroutine pfJacobian(snes,pf_vec,pf_jacob,pf_precond,dummy,ierr)
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
  Mat lhs_mat, pf_jacob, pf_precond
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pf_vec(:)

  PetscInt rowval,colval
  PetscScalar val

  integer, parameter :: imet = 1
  integer, parameter :: imkw = 2
  integer, parameter :: ipht = 3
  integer, parameter :: ipyr = 4
  integer, parameter :: ienv = 5

  integer, parameter :: no_fields = 5

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*psz*no_fields) :: vector_locator
  real*8 :: lnr,sqr

  integer :: io,jo


  real*8, dimension(:), allocatable :: A
  integer, dimension(:), allocatable :: JA, IA

  integer :: x, y, z   ! Loop variables
  integer :: linindex, contindex
  integer :: wrap

  !! Bulk free energy 
  real*8 :: f_met, f_mkw, f_pht, f_pyr, f_env
  real*8 :: w_met, w_mkw, w_pht, w_pyr, w_env

  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  real*8 :: sigma_pyr_met, sigma_pyr_mkw, sigma_pyr_pht, sigma_pyr_env
  real*8 :: sigma_met_pyr, sigma_mkw_pyr, sigma_pht_pyr, sigma_env_pyr

  real*8, dimension(psx,psy,psz+2) :: D_met, D_met_mkw, D_met_pht, D_met_pyr, D_met_env
  real*8, dimension(psx,psy,psz+2) :: D_mkw_met, D_mkw, D_mkw_pht, D_mkw_pyr, D_mkw_env
  real*8, dimension(psx,psy,psz+2) :: D_pht_met, D_pht_mkw, D_pht, D_pht_pyr, D_pht_env
  real*8, dimension(psx,psy,psz+2) :: D_pyr_met, D_pyr_mkw, D_pyr_pht, D_pyr, D_pyr_env
  real*8, dimension(psx,psy,psz+2) :: D_env_met, D_env_mkw, D_env_pht, D_env_pyr, D_env

  real*8, dimension(psx,psy,psz+2) :: loc_met, loc_mkw, loc_pht, loc_pyr, loc_env


  call calc_grad_pf()

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
        loc_met(x,y,1) = loc_met(x,y,2) ; loc_met(x,y,psz+2) = loc_met(x,y,psz+1) 
        loc_mkw(x,y,1) = loc_mkw(x,y,2) ; loc_mkw(x,y,psz+2) = loc_mkw(x,y,psz+1) 
        loc_pht(x,y,1) = loc_pht(x,y,2) ; loc_pht(x,y,psz+2) = loc_pht(x,y,psz+1) 
        loc_pyr(x,y,1) = loc_pyr(x,y,2) ; loc_pyr(x,y,psz+2) = loc_pyr(x,y,psz+1) 
        loc_env(x,y,1) = loc_env(x,y,2) ; loc_env(x,y,psz+2) = loc_env(x,y,psz+1) 
     end do
  end do


  do z = 1,psz+2
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


  allocate(A(((7*psx*psy*psz)-(2*psx*psy))*no_fields*no_fields))
  allocate(JA(((7*psx*psy*psz)-(2*psx*psy))*no_fields*no_fields))
  allocate(IA((no_fields*psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0


!!!! FOR MET
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx



           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z+1)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z+1)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z+1)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z+1)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T))*(13303/T))
           w_env = 0.0d0


           hill_met_mkw = 3.0d0*abs(w_met-w_mkw); hill_mkw_met = 3.0d0*abs(w_mkw-w_met)
           hill_met_pht = 3.0d0*abs(w_met-w_pht); hill_pht_met = 3.0d0*abs(w_pht-w_met)
           hill_met_pyr = 3.0d0*abs(w_met-w_pyr); hill_pyr_met = 3.0d0*abs(w_pyr-w_met)
           hill_met_env = 3.0d0*abs(w_met-w_env); hill_env_met = 3.0d0*abs(w_env-w_met)

           hill_mkw_pht = 3.0d0*abs(w_mkw-w_pht); hill_pht_mkw = 3.0d0*abs(w_pht-w_mkw)
           hill_mkw_pyr = 3.0d0*abs(w_mkw-w_pyr); hill_pyr_mkw = 3.0d0*abs(w_pyr-w_mkw)
           hill_mkw_env = 3.0d0*abs(w_mkw-w_env); hill_env_mkw = 3.0d0*abs(w_env-w_mkw)

           hill_pht_pyr = 3.0d0*abs(w_pht-w_pyr); hill_pyr_pht = 3.0d0*abs(w_pyr-w_pht)
           hill_pht_env = 3.0d0*abs(w_pht-w_env); hill_env_pht = 3.0d0*abs(w_env-w_pht)

           hill_pyr_env = 3.0d0*abs(w_pyr-w_env); hill_env_pyr = 3.0d0*abs(w_env-w_pyr)



           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imet-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!!!! MET-MET
           lnr = (2*M_met_mkw*hill_met_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (4*(w_met-w_mkw)*loc_mkw(x,y,z+1)*M_met_mkw) + &
                & (2*M_met_pht*hill_met_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (4*(w_met-w_pht)*loc_pht(x,y,z+1)*M_met_pht) + &
                & (2*M_met_pyr*hill_met_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (4*(w_met-w_pyr)*loc_pyr(x,y,z+1)*M_met_pyr) + &
                & (2*M_met_env*hill_met_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (4*(w_met-w_env)*loc_env(x,y,z+1)*M_met_env) 

           sqr = (2*M_met_mkw*(w_met-w_mkw)) - (M_met_mkw*2*hill_met_mkw*loc_mkw(x,y,z+1)) + &
                & (2*M_met_pht*(w_met-w_pht)) - (M_met_pht*2*hill_met_pht*loc_pht(x,y,z+1)) + &
                & (2*M_met_pyr*(w_met-w_pyr)) - (M_met_pyr*2*hill_met_pyr*loc_pyr(x,y,z+1)) + &
                & (2*M_met_env*(w_met-w_env)) - (M_met_env*2*hill_met_env*loc_env(x,y,z+1)) 

           sqr = sqr*2.0d0*loc_met(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr           
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
           lnr = (4*(w_met-w_mkw)*loc_met(x,y,z+1)*M_met_mkw) - (2*M_met_mkw*hill_met_mkw*loc_met(x,y,z+1)*loc_met(x,y,z+1))
           sqr = (2*M_met_mkw*(w_met-w_mkw)) + (M_met_mkw*2*hill_met_mkw*loc_met(x,y,z+1))
           sqr = sqr*2.0d0*loc_mkw(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_met-w_pht)*loc_met(x,y,z+1)*M_met_pht) - (2*M_met_pht*hill_met_pht*loc_met(x,y,z+1)*loc_met(x,y,z+1))
           sqr = (2*M_met_pht*(w_met-w_pht)) + (M_met_pht*2*hill_met_pht*loc_met(x,y,z+1))
           sqr = sqr*2.0d0*loc_pht(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_met-w_pyr)*loc_met(x,y,z+1)*M_met_pyr) - (2*M_met_pyr*hill_met_pyr*loc_met(x,y,z+1)*loc_met(x,y,z+1))
           sqr = (2*M_met_pyr*(w_met-w_pyr)) + (M_met_pyr*2*hill_met_pyr*loc_met(x,y,z+1))
           sqr = sqr*2.0d0*loc_pyr(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_met-w_env)*loc_met(x,y,z+1)*M_met_env) - (2*M_met_env*hill_met_env*loc_met(x,y,z+1)*loc_met(x,y,z+1))
           sqr = (2*M_met_env*(w_met-w_env)) + (M_met_env*2*hill_met_env*loc_met(x,y,z+1))
           sqr = sqr*2.0d0*loc_env(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z+1)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z+1)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z+1)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z+1)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T))*(13303/T))
           w_env = 0.0d0



           hill_met_mkw = 3.0d0*abs(w_met-w_mkw); hill_mkw_met = 3.0d0*abs(w_mkw-w_met)
           hill_met_pht = 3.0d0*abs(w_met-w_pht); hill_pht_met = 3.0d0*abs(w_pht-w_met)
           hill_met_pyr = 3.0d0*abs(w_met-w_pyr); hill_pyr_met = 3.0d0*abs(w_pyr-w_met)
           hill_met_env = 3.0d0*abs(w_met-w_env); hill_env_met = 3.0d0*abs(w_env-w_met)

           hill_mkw_pht = 3.0d0*abs(w_mkw-w_pht); hill_pht_mkw = 3.0d0*abs(w_pht-w_mkw)
           hill_mkw_pyr = 3.0d0*abs(w_mkw-w_pyr); hill_pyr_mkw = 3.0d0*abs(w_pyr-w_mkw)
           hill_mkw_env = 3.0d0*abs(w_mkw-w_env); hill_env_mkw = 3.0d0*abs(w_env-w_mkw)

           hill_pht_pyr = 3.0d0*abs(w_pht-w_pyr); hill_pyr_pht = 3.0d0*abs(w_pyr-w_pht)
           hill_pht_env = 3.0d0*abs(w_pht-w_env); hill_env_pht = 3.0d0*abs(w_env-w_pht)

           hill_pyr_env = 3.0d0*abs(w_pyr-w_env); hill_env_pyr = 3.0d0*abs(w_env-w_pyr)




           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((imkw-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! MKW-MET
           lnr = (4*(w_mkw-w_met)*loc_mkw(x,y,z+1)*M_mkw_met) - (2*M_mkw_met*hill_mkw_met*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1))
           sqr = (2*M_mkw_met*(w_mkw-w_met)) + (M_mkw_met*2*hill_mkw_met*loc_mkw(x,y,z+1))
           sqr = sqr*2.0d0*loc_met(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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


           lnr = (2*M_mkw_met*hill_mkw_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (4*(w_mkw-w_met)*loc_met(x,y,z+1)*M_mkw_met) + &
                & (2*M_mkw_pht*hill_mkw_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (4*(w_mkw-w_pht)*loc_pht(x,y,z+1)*M_mkw_pht) + &
                & (2*M_mkw_pyr*hill_mkw_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (4*(w_mkw-w_pyr)*loc_pyr(x,y,z+1)*M_mkw_pyr) + &
                & (2*M_mkw_env*hill_mkw_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (4*(w_mkw-w_env)*loc_env(x,y,z+1)*M_mkw_env) 

           sqr = (2*M_mkw_met*(w_mkw-w_met)) - (M_mkw_met*2*hill_mkw_met*loc_met(x,y,z+1)) + &
                & (2*M_mkw_pht*(w_mkw-w_pht)) - (M_mkw_pht*2*hill_mkw_pht*loc_pht(x,y,z+1)) + &
                & (2*M_mkw_pyr*(w_mkw-w_pyr)) - (M_mkw_pyr*2*hill_mkw_pyr*loc_pyr(x,y,z+1)) + &
                & (2*M_mkw_env*(w_mkw-w_env)) - (M_mkw_env*2*hill_mkw_env*loc_env(x,y,z+1)) 

           sqr = sqr*2.0d0*loc_mkw(x,y,z+1)


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
           A(contindex) = A(contindex) + lnr + sqr 
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
           lnr = (4*(w_mkw-w_pht)*loc_mkw(x,y,z+1)*M_mkw_pht) - (2*M_mkw_pht*hill_mkw_pht*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1))
           sqr = (2*M_mkw_pht*(w_mkw-w_pht)) + (M_mkw_pht*2*hill_mkw_pht*loc_mkw(x,y,z+1))
           sqr = sqr*2.0d0*loc_pht(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_mkw-w_pyr)*loc_mkw(x,y,z+1)*M_mkw_pyr) - (2*M_mkw_pyr*hill_mkw_pyr*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1))
           sqr = (2*M_mkw_pyr*(w_mkw-w_pyr)) + (M_mkw_pyr*2*hill_mkw_pyr*loc_mkw(x,y,z+1))
           sqr = sqr*2.0d0*loc_pyr(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_mkw-w_env)*loc_mkw(x,y,z+1)*M_mkw_env) - (2*M_mkw_env*hill_mkw_env*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1))
           sqr = (2*M_mkw_env*(w_mkw-w_env)) + (M_mkw_env*2*hill_mkw_env*loc_mkw(x,y,z+1))
           sqr = sqr*2.0d0*loc_env(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z+1)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z+1)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z+1)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z+1)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T))*(13303/T))
           w_env = 0.0d0


           hill_met_mkw = 3.0d0*abs(w_met-w_mkw); hill_mkw_met = 3.0d0*abs(w_mkw-w_met)
           hill_met_pht = 3.0d0*abs(w_met-w_pht); hill_pht_met = 3.0d0*abs(w_pht-w_met)
           hill_met_pyr = 3.0d0*abs(w_met-w_pyr); hill_pyr_met = 3.0d0*abs(w_pyr-w_met)
           hill_met_env = 3.0d0*abs(w_met-w_env); hill_env_met = 3.0d0*abs(w_env-w_met)

           hill_mkw_pht = 3.0d0*abs(w_mkw-w_pht); hill_pht_mkw = 3.0d0*abs(w_pht-w_mkw)
           hill_mkw_pyr = 3.0d0*abs(w_mkw-w_pyr); hill_pyr_mkw = 3.0d0*abs(w_pyr-w_mkw)
           hill_mkw_env = 3.0d0*abs(w_mkw-w_env); hill_env_mkw = 3.0d0*abs(w_env-w_mkw)

           hill_pht_pyr = 3.0d0*abs(w_pht-w_pyr); hill_pyr_pht = 3.0d0*abs(w_pyr-w_pht)
           hill_pht_env = 3.0d0*abs(w_pht-w_env); hill_env_pht = 3.0d0*abs(w_env-w_pht)

           hill_pyr_env = 3.0d0*abs(w_pyr-w_env); hill_env_pyr = 3.0d0*abs(w_env-w_pyr)



           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipht-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PHT-MET
           lnr = (4*(w_pht-w_met)*loc_pht(x,y,z+1)*M_pht_met) - (2*M_pht_met*hill_pht_met*loc_pht(x,y,z+1)*loc_pht(x,y,z+1))
           sqr = (2*M_pht_met*(w_pht-w_met)) + (M_pht_met*2*hill_pht_met*loc_pht(x,y,z+1))
           sqr = sqr*2.0d0*loc_met(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_pht-w_mkw)*loc_pht(x,y,z+1)*M_pht_mkw) - (2*M_pht_mkw*hill_pht_mkw*loc_pht(x,y,z+1)*loc_pht(x,y,z+1))
           sqr = (2*M_pht_mkw*(w_pht-w_mkw)) + (M_pht_mkw*2*hill_pht_mkw*loc_pht(x,y,z+1))
           sqr = sqr*2.0d0*loc_mkw(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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

           lnr = (2*M_pht_met*hill_pht_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (4*(w_pht-w_met)*loc_met(x,y,z+1)*M_pht_met) + &
                & (2*M_pht_mkw*hill_pht_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (4*(w_pht-w_mkw)*loc_mkw(x,y,z+1)*M_pht_mkw) + &
                & (2*M_pht_pyr*hill_pht_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (4*(w_pht-w_pyr)*loc_pyr(x,y,z+1)*M_pht_pyr) + &
                & (2*M_pht_env*hill_pht_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (4*(w_pht-w_env)*loc_env(x,y,z+1)*M_pht_env) 

           sqr = (2*M_pht_met*(w_pht-w_met)) - (M_pht_met*2*hill_pht_met*loc_met(x,y,z+1)) + &
                & (2*M_pht_mkw*(w_pht-w_mkw)) - (M_pht_mkw*2*hill_pht_mkw*loc_mkw(x,y,z+1)) + &
                & (2*M_pht_pyr*(w_pht-w_pyr)) - (M_pht_pyr*2*hill_pht_pyr*loc_pyr(x,y,z+1)) + &
                & (2*M_pht_env*(w_pht-w_env)) - (M_pht_env*2*hill_pht_env*loc_env(x,y,z+1)) 

           sqr = sqr*2.0d0*loc_pht(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr 
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
           lnr = (4*(w_pht-w_pyr)*loc_pht(x,y,z+1)*M_pht_pyr) - (2*M_pht_pyr*hill_pht_pyr*loc_pht(x,y,z+1)*loc_pht(x,y,z+1))
           sqr = (2*M_pht_pyr*(w_pht-w_pyr)) + (M_pht_pyr*2*hill_pht_pyr*loc_pht(x,y,z+1))
           sqr = sqr*2.0d0*loc_pyr(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_pht-w_env)*loc_pht(x,y,z+1)*M_pht_env) - (2*M_pht_env*hill_pht_env*loc_pht(x,y,z+1)*loc_pht(x,y,z+1))
           sqr = (2*M_pht_env*(w_pht-w_env)) + (M_pht_env*2*hill_pht_env*loc_pht(x,y,z+1))
           sqr = sqr*2.0d0*loc_env(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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




           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z+1)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z+1)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z+1)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z+1)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T))*(13303/T))
           w_env = 0.0d0


           hill_met_mkw = 3.0d0*abs(w_met-w_mkw); hill_mkw_met = 3.0d0*abs(w_mkw-w_met)
           hill_met_pht = 3.0d0*abs(w_met-w_pht); hill_pht_met = 3.0d0*abs(w_pht-w_met)
           hill_met_pyr = 3.0d0*abs(w_met-w_pyr); hill_pyr_met = 3.0d0*abs(w_pyr-w_met)
           hill_met_env = 3.0d0*abs(w_met-w_env); hill_env_met = 3.0d0*abs(w_env-w_met)

           hill_mkw_pht = 3.0d0*abs(w_mkw-w_pht); hill_pht_mkw = 3.0d0*abs(w_pht-w_mkw)
           hill_mkw_pyr = 3.0d0*abs(w_mkw-w_pyr); hill_pyr_mkw = 3.0d0*abs(w_pyr-w_mkw)
           hill_mkw_env = 3.0d0*abs(w_mkw-w_env); hill_env_mkw = 3.0d0*abs(w_env-w_mkw)

           hill_pht_pyr = 3.0d0*abs(w_pht-w_pyr); hill_pyr_pht = 3.0d0*abs(w_pyr-w_pht)
           hill_pht_env = 3.0d0*abs(w_pht-w_env); hill_env_pht = 3.0d0*abs(w_env-w_pht)

           hill_pyr_env = 3.0d0*abs(w_pyr-w_env); hill_env_pyr = 3.0d0*abs(w_env-w_pyr)



           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ipyr-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! PYR-MET
           lnr = (4*(w_pyr-w_met)*loc_pyr(x,y,z+1)*M_pyr_met) - (2*M_pyr_met*hill_pyr_met*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1))
           sqr = (2*M_pyr_met*(w_pyr-w_met)) + (M_pyr_met*2*hill_pyr_met*loc_pyr(x,y,z+1))
           sqr = sqr*2.0d0*loc_met(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_pyr-w_mkw)*loc_pyr(x,y,z+1)*M_pyr_mkw) - (2*M_pyr_mkw*hill_pyr_mkw*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1))
           sqr = (2*M_pyr_mkw*(w_pyr-w_mkw)) + (M_pyr_mkw*2*hill_pyr_mkw*loc_pyr(x,y,z+1))
           sqr = sqr*2.0d0*loc_mkw(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_pyr-w_pht)*loc_pyr(x,y,z+1)*M_pyr_pht) - (2*M_pyr_pht*hill_pyr_pht*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1))
           sqr = (2*M_pyr_pht*(w_pyr-w_pht)) + (M_pyr_pht*2*hill_pyr_pht*loc_pyr(x,y,z+1))
           sqr = sqr*2.0d0*loc_pht(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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

           lnr = (2*M_pyr_met*hill_pyr_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (4*(w_pyr-w_met)*loc_met(x,y,z+1)*M_pyr_met) + &
                         & (2*M_pyr_mkw*hill_pyr_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (4*(w_pyr-w_mkw)*loc_mkw(x,y,z+1)*M_pyr_mkw) + &
                         & (2*M_pyr_pht*hill_pyr_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (4*(w_pyr-w_pht)*loc_pht(x,y,z+1)*M_pyr_pht) + &
                         & (2*M_pyr_env*hill_pyr_env*loc_env(x,y,z+1)*loc_env(x,y,z+1)) + (4*(w_pyr-w_env)*loc_env(x,y,z+1)*M_pyr_env) 

           sqr = (2*M_pyr_met*(w_pyr-w_met)) - (M_pyr_met*2*hill_pyr_met*loc_met(x,y,z+1)) + &
                         & (2*M_pyr_mkw*(w_pyr-w_mkw)) - (M_pyr_mkw*2*hill_pyr_mkw*loc_mkw(x,y,z+1)) + &
                         & (2*M_pyr_pht*(w_pyr-w_pht)) - (M_pyr_pht*2*hill_pyr_pht*loc_pht(x,y,z+1)) + &
                         & (2*M_pyr_env*(w_pyr-w_env)) - (M_pyr_env*2*hill_pyr_env*loc_env(x,y,z+1)) 

           sqr = sqr*2.0d0*loc_pyr(x,y,z+1)


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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_pyr-w_env)*loc_pyr(x,y,z+1)*M_pyr_env) - (2*M_pyr_env*hill_pyr_env*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1))
           sqr = (2*M_pyr_env*(w_pyr-w_env)) + (M_pyr_env*2*hill_pyr_env*loc_pyr(x,y,z+1))
           sqr = sqr*2.0d0*loc_env(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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


           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401
           w_met = f_met - (mu(x,y,z+1)*140401*0.0015d0)

           !! Volume per mole of Mkw = 19.130 cc assuming a density of 4.28 g/cc from Lennie et. al. Mineralogical Magazine, December, Vol. 59, pp. 677-683
           !! Number of moles per m^3 = 48683
           f_mkw = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 65060)*48683.0d0
           w_mkw = f_mkw - (mu(x,y,z+1)*48683.0d0*0.80d0)  ! 0.95 because mackinawite is sligtly iron-rich and sulfur deficient

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z+1)*52275)

           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z+1)*mu(x,y,z+1)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z+1)*2*41667)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z+1)*(exp(mu(x,y,z+1)/(R*T))*(13303/T))
           w_env = 0.0d0





           hill_met_mkw = 3.0d0*abs(w_met-w_mkw); hill_mkw_met = 3.0d0*abs(w_mkw-w_met)
           hill_met_pht = 3.0d0*abs(w_met-w_pht); hill_pht_met = 3.0d0*abs(w_pht-w_met)
           hill_met_pyr = 3.0d0*abs(w_met-w_pyr); hill_pyr_met = 3.0d0*abs(w_pyr-w_met)
           hill_met_env = 3.0d0*abs(w_met-w_env); hill_env_met = 3.0d0*abs(w_env-w_met)

           hill_mkw_pht = 3.0d0*abs(w_mkw-w_pht); hill_pht_mkw = 3.0d0*abs(w_pht-w_mkw)
           hill_mkw_pyr = 3.0d0*abs(w_mkw-w_pyr); hill_pyr_mkw = 3.0d0*abs(w_pyr-w_mkw)
           hill_mkw_env = 3.0d0*abs(w_mkw-w_env); hill_env_mkw = 3.0d0*abs(w_env-w_mkw)

           hill_pht_pyr = 3.0d0*abs(w_pht-w_pyr); hill_pyr_pht = 3.0d0*abs(w_pyr-w_pht)
           hill_pht_env = 3.0d0*abs(w_pht-w_env); hill_env_pht = 3.0d0*abs(w_env-w_pht)

           hill_pyr_env = 3.0d0*abs(w_pyr-w_env); hill_env_pyr = 3.0d0*abs(w_env-w_pyr)





           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x + ((ienv-1)*psx*psy*psz)
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

!!! ENV-MET
           lnr = (4*(w_env-w_met)*loc_env(x,y,z+1)*M_env_met) - (2*M_env_met*hill_env_met*loc_env(x,y,z+1)*loc_env(x,y,z+1))
           sqr = (2*M_env_met*(w_env-w_met)) + (M_env_met*2*hill_env_met*loc_env(x,y,z+1))
           sqr = sqr*2.0d0*loc_met(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_env-w_mkw)*loc_env(x,y,z+1)*M_env_mkw) - (2*M_env_mkw*hill_env_mkw*loc_env(x,y,z+1)*loc_env(x,y,z+1))
           sqr = (2*M_env_mkw*(w_env-w_mkw)) + (M_env_mkw*2*hill_env_mkw*loc_env(x,y,z+1))
           sqr = sqr*2.0d0*loc_mkw(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_env-w_pht)*loc_env(x,y,z+1)*M_env_pht) - (2*M_env_pht*hill_env_pht*loc_env(x,y,z+1)*loc_env(x,y,z+1))
           sqr = (2*M_env_pht*(w_env-w_pht)) + (M_env_pht*2*hill_env_pht*loc_env(x,y,z+1))
           sqr = sqr*2.0d0*loc_pht(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           lnr = (4*(w_env-w_pyr)*loc_env(x,y,z+1)*M_env_pyr) - (2*M_env_pyr*hill_env_pyr*loc_env(x,y,z+1)*loc_env(x,y,z+1))
           sqr = (2*M_env_pyr*(w_env-w_pyr)) + (M_env_pyr*2*hill_env_pyr*loc_env(x,y,z+1))
           sqr = sqr*2.0d0*loc_pyr(x,y,z+1)

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
           A(contindex) = A(contindex) + lnr + sqr
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
           
           lnr = (2*M_env_met*hill_env_met*loc_met(x,y,z+1)*loc_met(x,y,z+1)) + (4*(w_env-w_met)*loc_met(x,y,z+1)*M_env_met) + &
                & (2*M_env_mkw*hill_env_mkw*loc_mkw(x,y,z+1)*loc_mkw(x,y,z+1)) + (4*(w_env-w_mkw)*loc_mkw(x,y,z+1)*M_env_mkw) + &
                & (2*M_env_pht*hill_env_pht*loc_pht(x,y,z+1)*loc_pht(x,y,z+1)) + (4*(w_env-w_pht)*loc_pht(x,y,z+1)*M_env_pht) + &
                & (2*M_env_pyr*hill_env_pyr*loc_pyr(x,y,z+1)*loc_pyr(x,y,z+1)) + (4*(w_env-w_pyr)*loc_pyr(x,y,z+1)*M_env_pyr)

           sqr = (2*M_env_met*(w_env-w_met)) - (M_env_met*2*hill_env_met*loc_met(x,y,z+1)) + &
                & (2*M_env_mkw*(w_env-w_mkw)) - (M_env_mkw*2*hill_env_mkw*loc_mkw(x,y,z+1)) + &
                & (2*M_env_pht*(w_env-w_pht)) - (M_env_pht*2*hill_env_pht*loc_pht(x,y,z+1)) + &
                & (2*M_env_pyr*(w_env-w_pyr)) - (M_env_pyr*2*hill_env_pyr*loc_pyr(x,y,z+1))

           sqr = sqr*2.0d0*loc_env(x,y,z+1)


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
           A(contindex) = A(contindex) + lnr + sqr
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

  IA((psx*psy*psz*no_fields)+1)= IA(psx*psy*psz*no_fields)+(6*no_fields)


  do linindex = 1,psx*psy*psz*no_fields
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

  do io = 1,psx*psy*psz*no_fields! Nrows

     rowval = io-1

     do jo = 1,IA(io+1)-IA(io)

        colval = JA(IA(io) + jo)
        val = A(IA(io) + jo)

        call MatSetValue(pf_precond,rowval,colval,val,INSERT_VALUES,ierr)

     end do

  end do

  call MatAssemblyBegin(pf_precond,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(pf_precond,MAT_FINAL_ASSEMBLY,ierr)

  return

end subroutine pfJacobian

