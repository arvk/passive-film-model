subroutine MetFunction(snes,met_vec,ret_vec,dummy,ierr)
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
  Vec met_vec,ret_vec
  Mat lhs_met_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_met_vec(:)

  real*8, dimension(psx,psy,psz+2) :: D_met

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


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)


           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D_met(x,y,z) = (M_met_pht*sigma_met_pht*pht(x,y,z)) + (M_met_pyr*sigma_met_pyr*pyr(x,y,z)) + (M_met_env*sigma_met_env*env(x,y,z))
        end do
     end do
  end do



  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)-1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y-1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x-1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_met(x,y,(z+1)+1)+D_met(x,y,(z+1)-1)+&
                &D_met(x,wrap(y+1,psy),z+1)+D_met(x,wrap(y-1,psy),z+1)+&
                &D_met(wrap(x+1,psx),y,z+1)+D_met(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_met(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(wrap(x+1,psx),y,z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_met(x,wrap(y+1,psy),z+1)+D_met(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_met(x,y,(z+1)+1)+D_met(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6








  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

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


  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_met_mat,ierr)
  call MatMult(lhs_met_mat,met_vec,ret_vec,ierr)


  call VecGetArrayF90(met_vec,point_met_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

       !!! Calculate linear prefactor
           lnr(linindex) = (M_met_pht*sigma_met_pht*del2pht(x,y,z+1)) + (2*M_met_pht*hill_met_pht*pht(x,y,z+1)*pht(x,y,z+1)) + (4*(w_met-w_pht)*pht(x,y,z+1)*M_met_pht) + &
                & (M_met_pyr*sigma_met_pyr*del2pyr(x,y,z+1)) + (2*M_met_pyr*hill_met_pyr*pyr(x,y,z+1)*pyr(x,y,z+1)) + (4*(w_met-w_pyr)*pyr(x,y,z+1)*M_met_pyr) + &
                & (M_met_env*sigma_met_env*del2env(x,y,z+1)) + (2*M_met_env*hill_met_env*env(x,y,z+1)*env(x,y,z+1)) + (4*(w_met-w_env)*env(x,y,z+1)*M_met_env) 

           lnr(linindex) = lnr(linindex)*point_met_vec(linindex)

       !!! Calculate quadratic prefactor
           sqr(linindex) = - ((2*M_met_pht*hill_met_pht*pht(x,y,z+1)) + (2*M_met_pht*(w_met-w_pht))) - &
                & (M_met_pyr*((2*hill_met_pyr*pyr(x,y,z+1)) + (2*M_met_pyr*(w_met-w_pyr)))) - &
                & (M_met_env*((2*hill_met_env*env(x,y,z+1)) + (2*M_met_env*(w_met-w_env))))

           sqr(linindex) = sqr(linindex)*point_met_vec(linindex)*point_met_vec(linindex)
        end do
     end do
  end do

  call VecRestoreArrayF90(met_vec,point_met_vec,ierr)

  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call MatDestroy(lhs_met_mat,ierr)

  return

end subroutine MetFunction







































subroutine PhtFunction(snes,pht_vec,ret_vec,dummy,ierr)
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
  Vec pht_vec,ret_vec
  Mat lhs_pht_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pht_vec(:)

  real*8, dimension(psx,psy,psz+2) :: D_pht

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


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)


           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D_pht(x,y,z) = (M_pht_met*sigma_pht_met*met(x,y,z)) + (M_pht_pyr*sigma_pht_pyr*pyr(x,y,z)) + (M_pht_env*sigma_pht_env*env(x,y,z))
        end do
     end do
  end do



  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)-1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y-1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x-1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_pht(x,y,(z+1)+1)+D_pht(x,y,(z+1)-1)+&
                &D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,wrap(y-1,psy),z+1)+&
                &D_pht(wrap(x+1,psx),y,z+1)+D_pht(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pht(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(wrap(x+1,psx),y,z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,wrap(y+1,psy),z+1)+D_pht(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pht(x,y,(z+1)+1)+D_pht(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6








  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

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


  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_pht_mat,ierr)
  call MatMult(lhs_pht_mat,pht_vec,ret_vec,ierr)


  call VecGetArrayF90(pht_vec,point_pht_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

       !!! Calculate linear prefactor
           lnr(linindex) = (M_pht_met*sigma_pht_met*del2met(x,y,z+1)) + (2*M_pht_met*hill_pht_met*met(x,y,z+1)*met(x,y,z+1)) + (4*M_pht_met*(w_pht-w_met)*met(x,y,z+1)) + &
                & (M_pht_pyr*sigma_pht_pyr*del2pyr(x,y,z+1)) + (2*M_pht_pyr*hill_pht_pyr*pyr(x,y,z+1)*pyr(x,y,z+1)) + (4*M_pht_pyr*(w_pht-w_pyr)*pyr(x,y,z+1)) + &
                & (M_pht_env*sigma_pht_env*del2env(x,y,z+1)) + (2*M_pht_env*hill_pht_env*env(x,y,z+1)*env(x,y,z+1)) + (4*M_pht_env*(w_pht-w_env)*env(x,y,z+1)) 

           lnr(linindex) = lnr(linindex)*point_pht_vec(linindex)

       !!! Calculate quadratic prefactor
           sqr(linindex) = -((2*M_pht_met*hill_pht_met*met(x,y,z+1)) + (2*M_pht_met*(w_pht-w_met))) - &
                & ((2*M_pht_pyr*hill_pht_pyr*pyr(x,y,z+1)) + (2*M_pht_pyr*(w_pht-w_pyr))) - &
                & ((2*M_pht_env*hill_pht_env*env(x,y,z+1)) + (2*M_pht_env*(w_pht-w_env)))

           sqr(linindex) = sqr(linindex)*point_pht_vec(linindex)*point_pht_vec(linindex)
        end do
     end do
  end do

  call VecRestoreArrayF90(pht_vec,point_pht_vec,ierr)

  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call MatDestroy(lhs_pht_mat,ierr)

  return

end subroutine PhtFunction





































subroutine PyrFunction(snes,pyr_vec,ret_vec,dummy,ierr)
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
  Vec pyr_vec,ret_vec
  Mat lhs_pyr_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_pyr_vec(:)

  real*8, dimension(psx,psy,psz+2) :: D_pyr

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


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)


           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D_pyr(x,y,z) = (M_pyr_met*sigma_pyr_met*met(x,y,z)) + (M_pyr_pht*sigma_pyr_pht*pht(x,y,z)) + (M_pyr_env*sigma_pyr_env*env(x,y,z))
        end do
     end do
  end do



  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)-1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y-1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x-1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_pyr(x,y,(z+1)+1)+D_pyr(x,y,(z+1)-1)+&
                &D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,wrap(y-1,psy),z+1)+&
                &D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_pyr(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(wrap(x+1,psx),y,z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,wrap(y+1,psy),z+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_pyr(x,y,(z+1)+1)+D_pyr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6








  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

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


  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_pyr_mat,ierr)
  call MatMult(lhs_pyr_mat,pyr_vec,ret_vec,ierr)


  call VecGetArrayF90(pyr_vec,point_pyr_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

       !!! Calculate linear prefactor
           lnr(linindex) = (M_pyr_met*sigma_pyr_met*del2met(x,y,z+1)) + (2*M_pyr_met*hill_pyr_met*met(x,y,z+1)*met(x,y,z+1)) + (4*M_pyr_met*(w_pyr-w_met)*met(x,y,z+1)) + &
                & (M_pyr_pht*sigma_pyr_pht*del2pht(x,y,z+1)) + (2*M_pyr_pht*hill_pyr_pht*pht(x,y,z+1)*pht(x,y,z+1)) + (4*M_pyr_pht*(w_pyr-w_pht)*pht(x,y,z+1)) + &
                & (M_pyr_env*sigma_pyr_env*del2env(x,y,z+1)) + (2*M_pyr_env*hill_pyr_env*env(x,y,z+1)*env(x,y,z+1)) + (4*M_pyr_env*(w_pyr-w_env)*env(x,y,z+1)) 

           lnr(linindex) = lnr(linindex)*point_pyr_vec(linindex)

       !!! Calculate quadratic prefactor
           sqr(linindex) = -((2*M_pyr_met*hill_pyr_met*met(x,y,z+1)) + (2*M_pyr_met*(w_pyr-w_met))) - &
                & ((2*M_pyr_pht*hill_pyr_pht*pht(x,y,z+1)) + (2*M_pyr_pht*(w_pyr-w_pht))) - &
                & ((2*M_pyr_env*hill_pyr_env*env(x,y,z+1)) + (2*M_pyr_env*(w_pyr-w_env)))

           sqr(linindex) = sqr(linindex)*point_pyr_vec(linindex)*point_pyr_vec(linindex)
        end do
     end do
  end do

  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)

  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call MatDestroy(lhs_pyr_mat,ierr)

  return

end subroutine PyrFunction





































































subroutine EnvFunction(snes,env_vec,ret_vec,dummy,ierr)
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
  Vec env_vec,ret_vec
  Mat lhs_env_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_env_vec(:)

  real*8, dimension(psx,psy,psz+2) :: D_env

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


  real*8 :: sigma_pht_pyr
  real*8 :: sigma_pyr_pht
  real*8 :: sigma_met_pyr
  real*8 :: sigma_pyr_met
  real*8 :: sigma_env_pyr
  real*8 :: sigma_pyr_env


  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1

           !! Volume per mole of Pht = 19.130 cc assuming a density of 4.6 g/cc
           !! Number of moles per m^3 = 52275
           f_pht = ((mu(x,y,z)*mu(x,y,z)*(1E-6)) + 20.53*T - 72050)*52275
           w_pht = f_pht - (mu(x,y,z)*52275)

           !! Volume per mole of air = 22.4 * (T/298) * 1E3 cc
           !! Number of moles per m^3 = 13303/T
           f_env = mu(x,y,z)*(exp(mu(x,y,z)/(R*T))*(13303/T))
           w_env = 0.0d0

           !! Volume per mole of Fe = 7.122 cc assuming a density of 7.84 g/cc
           !! Number of moles per m^3 = 140401
           f_met = 0.0d0*140401

!           w_met = f_met - (mu(x,y,z)*140401*0.0015d0)
           w_met = 0.0d0*140401


           !! Volumer per mole of Pyrite = 24 cc assuming a density of 5 g/cc
           !! Number of moles per m^3 = 41667
           f_pyr = ((mu(x,y,z)*mu(x,y,z)*(1E-9)) + 50.355*T - 98710)*41667
           w_pyr = f_pyr - (mu(x,y,z)*2*41667)


           sigma_pyr_met = sigma_pyr_met_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_pht = sigma_pyr_pht_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))
           sigma_pyr_env = sigma_pyr_env_0*(1+0.5*cos(4*(0.185+atan(delypyr(x,y,z)/(delzpyr(x,y,z)+1E-14)))-opyr(x,y,z)))

           sigma_met_pyr = sigma_pyr_met
           sigma_pht_pyr = sigma_pyr_pht
           sigma_env_pyr = sigma_pyr_env

        end do
     end do
  end do



  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D_env(x,y,z) = (M_env_met*sigma_env_met*met(x,y,z)) + (M_env_pyr*sigma_env_pyr*pyr(x,y,z)) + (M_env_pht*sigma_env_pht*pht(x,y,z))
        end do
     end do
  end do



  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)-1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y-1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x-1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_env(x,y,(z+1)+1)+D_env(x,y,(z+1)-1)+&
                &D_env(x,wrap(y+1,psy),z+1)+D_env(x,wrap(y-1,psy),z+1)+&
                &D_env(wrap(x+1,psx),y,z+1)+D_env(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D_env(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(wrap(x+1,psx),y,z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D_env(x,wrap(y+1,psy),z+1)+D_env(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D_env(x,y,(z+1)+1)+D_env(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6








  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

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


  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_env_mat,ierr)
  call MatMult(lhs_env_mat,env_vec,ret_vec,ierr)


  call VecGetArrayF90(env_vec,point_env_vec,ierr)
  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

       !!! Calculate linear prefactor
           lnr(linindex) = (M_env_met*sigma_env_met*del2met(x,y,z+1)) + (2*M_env_met*hill_env_met*met(x,y,z+1)*met(x,y,z+1)) + (4*M_env_met*(w_env-w_met)*met(x,y,z+1)) + &
                & (M_env_pyr*sigma_env_pyr*del2pyr(x,y,z+1)) + (2*M_env_pyr*hill_env_pyr*pyr(x,y,z+1)*pyr(x,y,z+1)) + (4*M_env_pyr*(w_env-w_pyr)*pyr(x,y,z+1)) + &
                & (M_env_pht*sigma_env_pht*del2pht(x,y,z+1)) + (2*M_env_pht*hill_env_pht*pht(x,y,z+1)*pht(x,y,z+1)) + (4*M_env_pht*(w_env-w_pht)*pht(x,y,z+1)) 

           lnr(linindex) = lnr(linindex)*point_env_vec(linindex)

       !!! Calculate quadratic prefactor
           sqr(linindex) = -((2*M_env_met*hill_env_met*met(x,y,z+1)) + (2*M_env_met*(w_env-w_met))) - &
                & ((2*M_env_pyr*hill_env_pyr*pyr(x,y,z+1)) + (2*M_env_pyr*(w_env-w_pyr))) - &
                & ((2*M_env_pht*hill_env_pht*pht(x,y,z+1)) + (2*M_env_pht*(w_env-w_pht)))

           sqr(linindex) = sqr(linindex)*point_env_vec(linindex)*point_env_vec(linindex)
        end do
     end do
  end do

  call VecRestoreArrayF90(env_vec,point_env_vec,ierr)

  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call MatDestroy(lhs_env_mat,ierr)

  return

end subroutine EnvFunction






