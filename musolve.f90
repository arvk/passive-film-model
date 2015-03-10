subroutine musolve(iter)
  use commondata
  use fields
  use laplacians
  use thermo_constants
  use diffusion_constants
  implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>


  PetscErrorCode ierr
  Vec mus_vec,rhs_vec
  Mat lhs_mat
  KSP ksp_mu
  PetscScalar, pointer :: point_mu_vec(:)

  integer, intent(in) :: iter

  integer :: x, y, z   ! Loop variables
  real*8 :: asd ! Anti-Surface-Diffusion

  !! Diffusivities
  real*8 :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  real*8, dimension(psx,psy,psz+2) :: D

  !! Derivative of sulfur density with chemical potential
  real*8 :: drho_dmu_met, drho_dmu_mkw, drho_dmu_pht, drho_dmu_pyr, drho_dmu_env, Chi

  integer, dimension(psx,psy) :: interface_loc
  real*8 :: noise

  real*8, dimension(psx,psy,psz+2) :: newmu
  integer :: wrap


  ! A/B/JA matrices for implicit solver
  real*8, dimension(psx*psy*psz) :: B
  real*8, dimension(psx*psy*psz) :: approxsol
  integer, dimension(psx*psy*psz) :: vector_locator
  real*8, dimension(psx*psy*psz) :: vecread
  real*8, dimension(psx*psy*psz) :: scratch1,scratch2,scratch3

  real*8, dimension(:), allocatable :: A, LU
  integer, dimension(:), allocatable :: JA, IA

  integer :: linindex, contindex
  integer :: iterations, solver_info


  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  newmu = 0.0d0

  if (mod(iter,swap_freq_pf).eq.1) then
     call swap_mu()
  end if

  open(unit = 55, file = '/dev/null')

  D_inter_met = max(D_Fe_met,D_S_met)
  D_inter_mkw = max(D_Fe_mkw,D_S_mkw)
  D_inter_pht = max(D_Fe_pht,D_S_pht)
  D_inter_pyr = max(D_Fe_pyr,D_S_pyr)
  D_inter_env = max(D_Fe_env,D_S_env)

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+2
           D(x,y,z) = ((met(x,y,z)*D_inter_met)+(mkw(x,y,z)*D_inter_mkw)+(pht(x,y,z)*D_inter_pht)+(pyr(x,y,z)*D_inter_pyr)+(env(x,y,z)*D_inter_env))
        end do
     end do
  end do


  approxsol = 0.0d0


  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           !! Calculate derivative of sulfur concentration with chemical potential
           drho_dmu_met = (0.0015d0*140401)/(R*T)
           drho_dmu_mkw = 0.95d0*48683.0d0/(2*250896.0d0)
           drho_dmu_pht = 52275.0d0/(2*250896.0d0)

           if (iter.lt.nomc/100) then
              drho_dmu_env = (0.0015d0*140401)/(R*T)
           else
              drho_dmu_env = (0.0015d0*(13303/T))/(R*T)
           end if

!           drho_dmu_pyr = 2*41667.0d0/(2*25089600.0d0)
           drho_dmu_pyr = (2*41667.0d0)/(2*250896.0d0)

           !! Calculate chemical 'specific heat'
           Chi = (met(x,y,z+1)*drho_dmu_met)+(mkw(x,y,z+1)*drho_dmu_mkw)+(pht(x,y,z+1)*drho_dmu_pht)+(pyr(x,y,z+1)*drho_dmu_pyr)+(env(x,y,z+1)*drho_dmu_env)

           B(linindex) =  (mu(x,y,z+1)/dt) - (((dpht_dt(x,y,z+1)*rho_pht) + (dmet_dt(x,y,z+1)*rho_met) + (dmkw_dt(x,y,z+1)*rho_mkw) + (denv_dt(x,y,z+1)*rho_env) + (dpyr_dt(x,y,z+1)*rho_pyr))/Chi)

           if (z.eq.1) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1)+D(x,y,z+1-1))/(dpf*dpf))*mu(x,y,z+1-1))
           elseif (z.eq.psz) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1+1)+D(x,y,z+1))/(dpf*dpf))*mu(x,y,z+1+1))
           end if

           approxsol(linindex) = mu(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do




  call VecCreate(PETSC_COMM_SELF,mus_vec,ierr)
  call VecSetSizes(mus_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(mus_vec,ierr)
  call VecSetUp(mus_vec,ierr)


  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*psz,ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*psz,vector_locator,B,INSERT_VALUES,ierr)


  !! Calculate the LHS (A matrix in Ax=B)
  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(LU((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))


  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0




  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D(x,y,(z+1)-1)+D(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(x,wrap(y-1,psy),z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(wrap(x-1,psx),y,z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D(x,y,(z+1)+1)+D(x,y,(z+1)-1)+&
                &D(x,wrap(y+1,psy),z+1)+D(x,wrap(y-1,psy),z+1)+&
                &D(wrap(x+1,psx),y,z+1)+D(wrap(x-1,psx),y,z+1)
           A(contindex) = A(contindex) + 6*D(x,y,z+1)
           A(contindex) = A(contindex)*0.5d0/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(wrap(x+1,psx),y,z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (0.5d0*(D(x,wrap(y+1,psy),z+1)+D(x,y,z+1)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (0.5d0*(D(x,y,(z+1)+1)+D(x,y,z+1)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6











  do linindex = 1,psx*psy*psz
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

  call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,psx*psy*psz,psx*psy*psz,IA,JA,A,lhs_mat,ierr)


  call KSPCreate(PETSC_COMM_SELF,ksp_mu,ierr)
  call KSPSetOperators(ksp_mu,lhs_mat,lhs_mat,ierr)
  call KSPSolve(ksp_mu,rhs_vec,mus_vec,ierr)



  write(55,*) "RANK", rank,iterations
  

  call VecGetArrayF90(mus_vec,point_mu_vec,ierr)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

              if (point_mu_vec(linindex).eq.point_mu_vec(linindex)) then
                 newmu(x,y,z+1) =  point_mu_vec(linindex)
                 else
                    write(6,*) 'Uh-oh' 
              end if

        end do
     end do
  end do

  call VecRestoreArrayF90(mus_vec,point_mu_vec,ierr)


  call VecDestroy(mus_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call MatDestroy(lhs_mat,ierr)
  call KSPDestroy(ksp_mu,ierr)


  !! Apply boundary conditions to chemical-potential-field update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           newmu(x,y,2) = mu(x,y,2)
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           newmu(x,y,psz+1) = mu(x,y,psz+1)
        end do
     end do
  end if

  !! Normalize chemical potential in bulk phases
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           if (env(x,y,z).gt.0.97) then
              newmu(x,y,z) = avg_mu_env
           end if
        end do
     end do
  end do


  !! Impose boundary counditions on the composition field (sulfidation rate)
  
  if (dt.gt.0) then

     interface_loc = 0

     rho_pht = 52275.0d0
     rho_met = 0.0015d0*140401       


     do x = 1,psx
        do y = 1,psy
           do z = psz+1,2,-1
              if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then
                 interface_loc(x,y) = z
                 newmu(x,y,interface_loc(x,y)) = mu(x,y,interface_loc(x,y)) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt)/(dpf))
                 newmu(x,y,interface_loc(x,y)) = min(newmu(x,y,interface_loc(x,y)),max_mu)
                 newmu(x,y,interface_loc(x,y)+1) = newmu(x,y,interface_loc(x,y)-1)
                 exit
              end if
           end do
        end do
     end do

  end if

  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+1
           mu(x,y,z) = newmu(x,y,z)
        end do
     end do
  end do



close(55)


end subroutine musolve
