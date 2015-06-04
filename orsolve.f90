subroutine orsolve(iter)
  use commondata
  use fields
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  integer, intent(in) :: iter                 ! Iteration count
  integer :: x, y, z                          ! Index for x-, y-, and z-direction (Loop)
  integer :: wrap                             ! Wrapping along the x-, y-, and z-direction
  integer :: odiff                            ! Calcualte difference in orientation between two points
  real*8, parameter :: infinitesimal = 1E-15  ! A hard-coded 'small' number

  ! Orientation field mobilities and diffusivities
  real*8, dimension(psx,psy,(psz+(2*ghost_width)-2)) :: M_opyr
  real*8 :: M_opyr_max = 1.0E-16
  real*8 :: M_opyr_min = 1.0E-19
  real*8, dimension(psx,psy,(psz+(2*ghost_width)-2)) :: D_opyr
  real*8, dimension(psx,psy,(psz+(2*ghost_width)-2)) :: del_opyr
  real*8 :: D_opyr_max = 1E1   ! Truncation for the orientation field

  ! Equation and solution matrices and vectors
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: B, approxsol
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: vector_locator
  real*8, dimension(:), allocatable :: A
  integer, dimension(:), allocatable :: JA, IA
  integer :: linindex, contindex
  real*8 :: orc,orc1,orc2,orc3,orc4,orc5,orc6

  ! Sorting of IA/JA matrices
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  PetscErrorCode ierr
  Vec or_vec,rhs_vec
  Mat lhs_mat
  KSP ksp_or
  PetscScalar, pointer :: point_or_vec(:)
  KSPConvergedReason or_converged_reason

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call swap_or()

  M_opyr = 0.0d0   !! Initialize to 0

  do x = 1,psx
     do y = 1,psy
        do z = 1,(psz+(2*ghost_width)-2)

           M_opyr(x,y,z) = min(M_opyr_min/(((pyr(x,y,z+1))**4)+infinitesimal),M_opyr_max)
           D_opyr(x,y,z) = D_opyr_max

        end do
     end do
  end do


  approxsol = 0.0d0

  !! Calculate the RHS (B matrix in Ax=B)
  !! B = C_{i,j} + o_pyr/dt

  orc = 0.0d0

  do z = 1,(psz+(2*ghost_width)-2)
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           orc = odiff(opyr(x,y,z+1),opyr(wrap(x-1,psx),y,z+1))-(opyr(x,y,z+1)-opyr(wrap(x-1,psx),y,z+1))
           orc1 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(wrap(x-1,psx),y,z))/(dpf*dpf)

           orc = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(x,y,z+1))-(opyr(wrap(x+1,psx),y,z+1)-opyr(x,y,z+1))
           orc2 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1),opyr(x,wrap(y-1,psy),z+1))-(opyr(x,y,z+1)-opyr(x,wrap(y-1,psy),z+1))
           orc3 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,wrap(y-1,psy),z))/(dpf*dpf)

           orc = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,y,z+1))-(opyr(x,wrap(y+1,psy),z+1)-opyr(x,y,z+1))
           orc4 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1),opyr(x,y,z+1-1))-(opyr(x,y,z+1)-opyr(x,y,z+1-1))
           orc5 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1))-(opyr(x,y,z+1+1)-opyr(x,y,z+1))
           orc6 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,z))/(dpf*dpf)

           orc = (orc1+orc3+orc5)-(orc2+orc4+orc6)

           orc = 0.0d0

           B(linindex) = (opyr(x,y,z+1)/dt) + orc

           if (z.eq.1) then
              B(linindex) = B(linindex) + ((M_opyr(x,y,z)*0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf))*opyr(x,y,z+1-1))
           elseif (z.eq.(psz+(2*ghost_width)-2)) then
              B(linindex) = B(linindex) + ((M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,z))/(dpf*dpf))*opyr(x,y,z+1+1))
           end if

           approxsol(linindex) = opyr(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do



  call VecCreate(PETSC_COMM_SELF,or_vec,ierr)
  call VecSetSizes(or_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2),ierr)
  call VecSetFromOptions(or_vec,ierr)
  call VecSetUp(or_vec,ierr)
  call VecSetValues(or_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,approxsol,INSERT_VALUES,ierr)

  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2),ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,B,INSERT_VALUES,ierr)



  !! Calculate the LHS (A matrix in Ax=B)
  allocate(A((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(JA((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(IA((psx*psy*(psz+(2*ghost_width)-2))+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0

  do z = 1,(psz+(2*ghost_width)-2)
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z-1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,z)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,(psz+(2*ghost_width)-2))-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y-1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x-1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_opyr(x,y,wrap(z+1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,wrap(z-1,(psz+(2*ghost_width)-2)))+&
                &D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,wrap(y-1,psy),z)+&
                &D_opyr(wrap(x+1,psx),y,z)+D_opyr(wrap(x-1,psx),y,z)
           A(contindex) = A(contindex) + 6*D_opyr(x,y,z)
           A(contindex) = A(contindex)*0.5d0*M_opyr(x,y,z)/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z+1,(psz+(2*ghost_width)-2)))+D_opyr(x,y,z)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,(psz+(2*ghost_width)-2))-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*(psz+(2*ghost_width)-2))+1)= IA(psx*psy*(psz+(2*ghost_width)-2))+6

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

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

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*(psz+(2*ghost_width)-2),psx*psy*(psz+(2*ghost_width)-2),IA,JA,A,lhs_mat,ierr)

  call KSPCreate(PETSC_COMM_SELF,ksp_or,ierr)
  call KSPSetOperators(ksp_or,lhs_mat,lhs_mat,ierr)
  call KSPSetInitialGuessNonzero(ksp_or,PETSC_TRUE,ierr)
  call KSPSetFromOptions(ksp_or,ierr)
  call KSPSetUp(ksp_or,ierr)
  call KSPSolve(ksp_or,rhs_vec,or_vec,ierr)

  call KSPGetConvergedReason(ksp_or,or_converged_reason,ierr)

  if (or_converged_reason.gt.0) then     ! If Solver is converged

  call VecGetArrayF90(or_vec,point_or_vec,ierr)

  do z = 1,(psz+(2*ghost_width)-2)
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           opyr(x,y,z+1) = (pyr(x,y,z+1)*point_or_vec(linindex)) + ((1.0d0-pyr(x,y,z+1))*opyr(x,y,z+1))

        end do
     end do
  end do

  call VecRestoreArrayF90(or_vec,point_or_vec,ierr)

  end if

  call VecDestroy(or_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call MatDestroy(lhs_mat,ierr)
  call KSPDestroy(ksp_or,ierr)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  do x = 1,psx
     do y = 1,psy
        do z = 3,psz+(2*ghost_width)-2
           if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then  ! If it is iron(sulfide)/environment boundary, set the gradient of the OPYR field to zero there
              opyr(x,y,z-2) = opyr(x,y,z)
              opyr(x,y,z-1) = opyr(x,y,z)
              opyr(x,y,z+1) = opyr(x,y,z)
              opyr(x,y,z+2) = opyr(x,y,z)
              exit
           end if
        end do
     end do
  end do



end subroutine orsolve





double precision function odiff(or1,or2)
  implicit none
  real*8, intent(in) :: or1,or2
  real*8 :: Pi = 3.14159265d0

  if (or1>or2) then

     if ((or1-or2).lt.(or2+(Pi/2.0d0)-or1)) then
        odiff = -(or1-or2)
     else
        odiff = (or2+(Pi/2.0d0)-or1)
     end if

  else

     if ((or2-or1).lt.(or1+(Pi/2.0d0)-or2)) then
        odiff = or2-or1
     else
        odiff = -(or1+(Pi/2.0d0)-or2)
     end if

  end if

  return
end function odiff




integer function wrap(a,lim)
  implicit none
  integer, intent(in) :: a,lim
  wrap = modulo(a-1,lim)+1
  return
end function wrap
