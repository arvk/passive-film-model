subroutine electroFunction(snes,elpot_vec,ret_vec,dummy,ierr)
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
  Vec elpot_vec,ret_vec
  Mat lhs_mat
  integer dummy(*)
  PetscErrorCode ierr
  PetscScalar, pointer :: point_elpot_vec(:)

  ! A/B/JA matrices for implicit solver
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: vector_locator
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: const,lnr,sqr

  real*8, dimension(:), allocatable :: A
  integer, dimension(:), allocatable :: JA, IA

  integer :: x, y, z   ! Loop variables
  integer :: linindex, contindex
  integer :: wrap

  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  real*8 :: c0 = 1.0d0 !! Moles/m^3
  real*8 :: el_charg = 1.60217657E-19*6.022E23

  real*8 :: epsilon0, epsilon_met, epsilon_mkw, epsilon_pht, epsilon_pyr, epsilon_env

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: loc_elpot

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: epsilonr

  real*8 :: exponent
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: nonlin

  call VecGetArrayF90(elpot_vec,point_elpot_vec,ierr)
  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           loc_elpot(x,y,z+1) = point_elpot_vec(linindex)

        end do
     end do
  end do
  call VecRestoreArrayF90(elpot_vec,point_elpot_vec,ierr)

  do z = 1,ghost_width
     do y = 1,psy
        do x = 1,psx
           loc_elpot(x,y,z) = loc_elpot(x,y,1+ghost_width) ; loc_elpot(x,y,psz+z+ghost_width) = loc_elpot(x,y,psz+ghost_width)
        end do
     end do
  end do


  epsilon_met = 500.0d0
  epsilon_mkw = 500.0d0
  epsilon_pht = 2.0d0
  epsilon_pyr = 11.0d0
  epsilon_env = 80.0d0

  epsilon0 = 8.854187817E-12 !! Define vacuum permittivity

  do z = 1,psz+(2*ghost_width)
     do y = 1,psy
        do x = 1,psx
           epsilonr(x,y,z) = epsilon_met*met(x,y,z) + epsilon_mkw*mkw(x,y,z) + epsilon_pht*pht(x,y,z) + epsilon_pyr*pyr(x,y,z) + epsilon_env*env(x,y,z)
           epsilonr(x,y,z) = epsilonr(x,y,z)*epsilon0*(1.0d0-voids(x,y,z))
        end do
     end do
  end do



  allocate(A((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(JA((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(IA((psx*psy*(psz+(2*ghost_width)-2))+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0; linindex = 0

  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           vector_locator(linindex) = linindex-1

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = (0.5d0*(epsilonr(x,y,(z+1)-1)+epsilonr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = (0.5d0*(epsilonr(x,wrap(y-1,psy),z+1)+epsilonr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = (0.5d0*(epsilonr(wrap(x-1,psx),y,z+1)+epsilonr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = epsilonr(x,y,(z+1)+1)+epsilonr(x,y,(z+1)-1)+&
                &epsilonr(x,wrap(y+1,psy),z+1)+epsilonr(x,wrap(y-1,psy),z+1)+&
                &epsilonr(wrap(x+1,psx),y,z+1)+epsilonr(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*epsilonr(x,y,z+1)
           A(contindex) = A(contindex)*(-0.5d0)/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = (0.5d0*(epsilonr(wrap(x+1,psx),y,z+1)+epsilonr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = (0.5d0*(epsilonr(x,wrap(y+1,psy),z+1)+epsilonr(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz+(2*ghost_width)-2) then
              contindex = contindex + 1
              A(contindex) = (0.5d0*(epsilonr(x,y,(z+1)+1)+epsilonr(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz+(2*ghost_width)-2)-1)*psx*psy) + ((y-1)*psx) + x
           end if

        end do
     end do
  end do

  IA((psx*psy*(psz+(2*ghost_width)-2))+1)= IA(psx*psy*(psz+(2*ghost_width)-2))+6


  do linindex = 1,psx*psy*(psz+(2*ghost_width)-2)
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

  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*(psz+(2*ghost_width)-2),psx*psy*(psz+(2*ghost_width)-2),IA,JA,A,lhs_mat,ierr)
  call MatMult(lhs_mat,elpot_vec,ret_vec,ierr)

  lnr = 0.0d0
  sqr = 0.0d0

  do z = 1,psz+(2*ghost_width)-2
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           exponent = (loc_elpot(x,y,z+1)*96485)/(R*T)
           nonlin(linindex) = c0*el_charg*(exp(0.0d0-exponent)-exp(0.0d0+exponent))*env(x,y,z+1)

        end do
     end do
  end do

  call VecSetValues(ret_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,nonlin,ADD_VALUES,ierr)

  call VecAssemblyBegin(ret_vec,ierr)
  call VecAssemblyEnd(ret_vec,ierr)

  call MatDestroy(lhs_mat,ierr)

  return

end subroutine electroFunction
