subroutine pHsolve(iter)
  
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  integer, intent(in) :: iter  ! Iteration count
  integer :: x, y, z           ! Index for x-, y-, and z-direction (Loop)
  integer :: wrap              ! Wrapping along the x-, y-, and z-direction
  integer, dimension(psx,psy) :: interface_loc
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: newpH

  !! Diffusivities
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D
  real*8 :: interface_multiply ! Flag to control interfacial diffusion

  ! Equation and solution matrices and vectors
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: B, approxsol
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: vector_locator
  real*8, dimension(:), allocatable :: A
  integer, dimension(:), allocatable :: JA, IA
  integer :: linindex, contindex
  integer :: iterations, solver_info

  ! Sorting of IA/JA matrices
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  PetscErrorCode ierr
  Vec pH_vec,rhs_vec
  Mat lhs_mat
  KSP ksp_pH
  PetscScalar, pointer :: point_pH_vec(:)
  KSPConvergedReason pH_converged_reason

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  call swap_pH()

  newpH = pH

  D = 0.0d0

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+(2*ghost_width)

           D(x,y,z) = max(min(met(x,y,z)-0.005d0,1.0d0),0.0d0)*D_H_met
           D(x,y,z) = D(x,y,z) + max(min(mkw(x,y,z)-0.005d0,1.0d0),0.0d0)*D_H_mkw
           D(x,y,z) = D(x,y,z) + max(min(pht(x,y,z)-0.005d0,1.0d0),0.0d0)*D_H_pht
           D(x,y,z) = D(x,y,z) + max(min(pyr(x,y,z)-0.005d0,1.0d0),0.0d0)*D_H_pyr
           D(x,y,z) = D(x,y,z) + max(min(env(x,y,z)-0.005d0,1.0d0),0.0d0)*D_H_env
           D(x,y,z) = D(x,y,z)*(1.0d0-voids(x,y,z))

        end do
     end do
  end do


  !! Identify the iron(sulfide)/environment boundary
  interface_loc = 0
  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then
              interface_loc(x,y) = z
              exit
           end if
        end do
     end do
  end do


  approxsol = 0.0d0

  do z = 1,(psz+(2*ghost_width)-2)
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           B(linindex) =  (pH(x,y,z+1)/dt) 

           if (z.eq.1) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1)+D(x,y,z+1-1))/(dpf*dpf))*pH(x,y,z+1-1))
           elseif (z.eq.(psz+(2*ghost_width)-2)) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1+1)+D(x,y,z+1))/(dpf*dpf))*pH(x,y,z+1+1))
           end if

           approxsol(linindex) = pH(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do




  call VecCreate(PETSC_COMM_SELF,pH_vec,ierr)
  call VecSetSizes(pH_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2),ierr)
  call VecSetFromOptions(pH_vec,ierr)
  call VecSetUp(pH_vec,ierr)
  call VecSetValues(pH_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,approxsol,INSERT_VALUES,ierr)

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

           interface_multiply = 1.0d0
           if ((z.gt.interface_loc(x,y)-2).and.(z.lt.interface_loc(x,y)+1)) then
              interface_multiply = 0.0d0
           end if

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D(x,y,(z+1)-1)+D(x,y,z+1)))/(dpf*dpf)
              A(contindex) = A(contindex)*interface_multiply
              JA(contindex) = ((wrap(z-1,(psz+(2*ghost_width)-2))-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(x,wrap(y-1,psy),z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(wrap(x-1,psx),y,z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = (D(x,y,(z+1)+1) + D(x,y,(z+1)-1) + 2*D(x,y,z+1))*interface_multiply
           A(contindex) = A(contindex)+&
                &D(x,wrap(y+1,psy),z+1)+D(x,wrap(y-1,psy),z+1)+&
                &D(wrap(x+1,psx),y,z+1)+D(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 4*D(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(wrap(x+1,psx),y,z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(x,wrap(y+1,psy),z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. (psz+(2*ghost_width)-2)) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D(x,y,(z+1)+1)+D(x,y,z+1)))/(dpf*dpf)
              A(contindex) = A(contindex)*interface_multiply
              JA(contindex) = ((wrap(z+1,(psz+(2*ghost_width)-2))-1)*psx*psy) + ((y-1)*psx) + x
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

  call KSPCreate(PETSC_COMM_SELF,ksp_pH,ierr)
  call KSPSetOperators(ksp_pH,lhs_mat,lhs_mat,ierr)
  call KSPSetInitialGuessNonzero(ksp_pH,PETSC_TRUE,ierr)
  call KSPSetFromOptions(ksp_pH,ierr)
  call KSPSetUp(ksp_pH,ierr)
  call KSPSolve(ksp_pH,rhs_vec,pH_vec,ierr)

  call KSPGetConvergedReason(ksp_pH,pH_converged_reason,ierr)

  if (pH_converged_reason.gt.0) then  ! If Solver is converged

     call VecGetArrayF90(pH_vec,point_pH_vec,ierr)

     do z = 1,(psz+(2*ghost_width)-2)
        do y = 1,psy
           do x = 1,psx
              linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

              newpH(x,y,z+1) =  point_pH_vec(linindex)

           end do
        end do
     end do

     call VecRestoreArrayF90(pH_vec,point_pH_vec,ierr)

  else                                ! If Solver is not converged

     write(6,'(A,I2,A,I6,A,I2)') " ERROR: pH field evolution did not converge in rank ", rank, " and iteration ", iter, ". ERROR CODE: ", pH_converged_reason

     newpH = pH

  end if


  call VecDestroy(pH_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call MatDestroy(lhs_mat,ierr)
  call KSPDestroy(ksp_pH,ierr)



  if (rank.eq.0) then          ! Apply global boundary conditions to process 0 (Metal end)
     do x = 1,psx
        do y = 1,psy
           newpH(x,y,1+ghost_width) = newpH(x,y,2+ghost_width)
        end do
     end do
  elseif(rank.eq.procs-1) then ! Apply global boundary conditions to process procs-1 (Environment end)
     do x = 1,psx
        do y = 1,psy
           newpH(x,y,psz+ghost_width) = pH(x,y,psz+ghost_width)
        end do
     end do
  end if


!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           pH(x,y,z) = newpH(x,y,z)
        end do
     end do
  end do


end subroutine pHsolve
