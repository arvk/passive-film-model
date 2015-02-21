subroutine F_Pyr_Met(snes,pyr_vec,ret_vec,dummy,ierr)
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

  integer dummy(*)
  SNES snes
  Vec pyr_vec,ret_vec
  Mat lhs_pyr_mat
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
  real*8 :: f_met, f_pyr, f_pht, f_env
  real*8 :: w_met, w_pyr, w_pht, w_env

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D_pyr(x,y,z) = M_pyr_met*sigma_pyr_met*met(x,y,z)
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

           include 'define_f_w.user'

       !!! Calculate linear prefactor
           lnr(linindex) = (4*w_pyr*met(x,y,z+1)*M_pyr_met)
           lnr(linindex) = lnr(linindex)*point_pyr_vec(linindex)

       !!! Calculate quadratic prefactor
           sqr(linindex) = (2*M_pyr_met*w_pyr) - (2*M_pyr_met*hill_pyr_met*met(x,y,z+1))
           sqr(linindex) = sqr(linindex)*point_pyr_vec(linindex)*point_pyr_vec(linindex)

        end do
     end do
  end do

  call VecRestoreArrayF90(pyr_vec,point_pyr_vec,ierr)

  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,lnr,ADD_VALUES,ierr)
  call VecSetValues(ret_vec,psx*psy*psz,vector_locator,sqr,ADD_VALUES,ierr)

  call VecAssemblyBegin(ret_vec,ierr)
  call VecAssemblyEnd(ret_vec,ierr)

  call MatDestroy(lhs_pyr_mat,ierr)

  return

end subroutine F_Pyr_Met
