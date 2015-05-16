subroutine musolve(iter)
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


  PetscErrorCode ierr
  Vec mus_vec,rhs_vec
  Mat lhs_mat
  KSP ksp_mu
  PetscScalar, pointer :: point_mu_vec(:)
  KSPConvergedReason mu_converged_reason

  integer, intent(in) :: iter

  integer :: x, y, z   ! Loop variables
  real*8 :: interface_multiply ! Flag to control interfacial diffusion
  real*8 :: asd ! Anti-Surface-Diffusion

  !! Diffusivities
  real*8 :: D_inter_met, D_inter_mkw, D_inter_pht, D_inter_pyr, D_inter_env
  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: D

  !! Derivative of sulfur density with chemical potential
  real*8 :: drho_dmu_met, drho_dmu_mkw, drho_dmu_pht, drho_dmu_pyr, drho_dmu_env, Chi

  integer, dimension(psx,psy) :: interface_loc
  real*8 :: noise

  real*8, dimension(psx,psy,psz+(2*ghost_width)) :: newmu
  integer :: wrap

  real*8 :: sulf_rate_gas_met, sulf_rate_gas_pht, sulf_rate_gas_pyr
  real*8 :: sulf_rate_liq_met, sulf_rate_liq_mkw, sulf_rate_liq_pht, sulf_rate_liq_pyr

  ! A/B/JA matrices for implicit solver
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: B
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: approxsol
  integer, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: vector_locator
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: vecread
  real*8, dimension(psx*psy*(psz+(2*ghost_width)-2)) :: scratch1,scratch2,scratch3

  real*8, dimension(:), allocatable :: A, LU
  integer, dimension(:), allocatable :: JA, IA

  integer :: linindex, contindex
  integer :: iterations, solver_info


  ! Sorting for A/B/JA
  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap


  newmu = 0.0d0

  D_inter_met = D_S_met
  D_inter_mkw = max(D_Fe_mkw,D_S_mkw)
  D_inter_pht = max(D_Fe_pht,D_S_pht)
  D_inter_pyr = max(D_Fe_pyr,D_S_pyr)
  D_inter_env = D_S_env

  !! Calculate derivative of sulfur concentration with chemical potential
  drho_dmu_met = (0.0015d0*140401)/(R*T)
  drho_dmu_mkw = 0.95d0*48683.0d0/(2*250896.0d0)
  drho_dmu_pht = 52275.0d0/(2*250896.0d0)

  if (iter.lt.nomc/100) then
     drho_dmu_env = (0.0015d0*140401)/(R*T)
  else
     drho_dmu_env = (0.0015d0*(13303/T))/(R*T)
  end if

  drho_dmu_pyr = (2*41667.0d0)/(2*250896.0d0)

  D = 0.0d0
  Chi = 0.0d0

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz+(2*ghost_width)

           !! Calculate chemical 'specific heat'
           Chi = max(min(met(x,y,z)-0.005d0,1.0d0),0.0d0)*drho_dmu_met
           Chi = Chi + max(min(mkw(x,y,z)-0.005d0,1.0d0),0.0d0)*drho_dmu_mkw
           Chi = Chi + max(min(pht(x,y,z)-0.005d0,1.0d0),0.0d0)*drho_dmu_pht
           Chi = Chi + max(min(pyr(x,y,z)-0.005d0,1.0d0),0.0d0)*drho_dmu_pyr
           Chi = Chi + max(min(env(x,y,z)-0.005d0,1.0d0),0.0d0)*drho_dmu_env

           D(x,y,z) = max(min(met(x,y,z)-0.005d0,1.0d0),0.0d0)*D_inter_met
           D(x,y,z) = D(x,y,z) + max(min(mkw(x,y,z)-0.005d0,1.0d0),0.0d0)*D_inter_mkw
           D(x,y,z) = D(x,y,z) + max(min(pht(x,y,z)-0.005d0,1.0d0),0.0d0)*D_inter_pht
           D(x,y,z) = D(x,y,z) + max(min(pyr(x,y,z)-0.005d0,1.0d0),0.0d0)*D_inter_pyr
           D(x,y,z) = D(x,y,z) + max(min(env(x,y,z)-0.005d0,1.0d0),0.0d0)*D_inter_env
           D(x,y,z) = D(x,y,z)/Chi
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
           !! Calculate chemical 'specific heat'
           Chi = max(min(met(x,y,z+1)-0.005d0,1.0d0),0.0d0)*drho_dmu_met
           Chi = Chi + max(min(mkw(x,y,z+1)-0.005d0,1.0d0),0.0d0)*drho_dmu_mkw
           Chi = Chi + max(min(pht(x,y,z+1)-0.005d0,1.0d0),0.0d0)*drho_dmu_pht
           Chi = Chi + max(min(pyr(x,y,z+1)-0.005d0,1.0d0),0.0d0)*drho_dmu_pyr
           Chi = Chi + max(min(env(x,y,z+1)-0.005d0,1.0d0),0.0d0)*drho_dmu_env

           B(linindex) =  (mu(x,y,z+1)/dt) - (((dpht_dt(x,y,z+1)*rho_pht) + (dmet_dt(x,y,z+1)*rho_met) + (dmkw_dt(x,y,z+1)*rho_mkw) + (denv_dt(x,y,z+1)*rho_env) + (dpyr_dt(x,y,z+1)*rho_pyr))/Chi)

           if (z.eq.1) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1)+D(x,y,z+1-1))/(dpf*dpf))*mu(x,y,z+1-1))
           elseif (z.eq.(psz+(2*ghost_width)-2)) then
              B(linindex) = B(linindex) + ((0.5d0*(D(x,y,z+1+1)+D(x,y,z+1))/(dpf*dpf))*mu(x,y,z+1+1))
           end if

           approxsol(linindex) = mu(x,y,z+1)

           vector_locator(linindex) = linindex-1

        end do
     end do
  end do




  call VecCreate(PETSC_COMM_SELF,mus_vec,ierr)
  call VecSetSizes(mus_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2),ierr)
  call VecSetFromOptions(mus_vec,ierr)
  call VecSetUp(mus_vec,ierr)
  call VecSetValues(mus_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,approxsol,INSERT_VALUES,ierr)


  call VecCreate(PETSC_COMM_SELF,rhs_vec,ierr)
  call VecSetSizes(rhs_vec,PETSC_DECIDE,psx*psy*(psz+(2*ghost_width)-2),ierr)
  call VecSetFromOptions(rhs_vec,ierr)
  call VecSetUp(rhs_vec,ierr)
  call VecSetValues(rhs_vec,psx*psy*(psz+(2*ghost_width)-2),vector_locator,B,INSERT_VALUES,ierr)


  !! Calculate the LHS (A matrix in Ax=B)
  allocate(A((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(JA((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
  allocate(LU((7*psx*psy*(psz+(2*ghost_width)-2))-(2*psx*psy)))
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


  call KSPCreate(PETSC_COMM_SELF,ksp_mu,ierr)
  call KSPSetOperators(ksp_mu,lhs_mat,lhs_mat,ierr)
  call KSPSetInitialGuessNonzero(ksp_mu,PETSC_TRUE,ierr)
  call KSPSetFromOptions(ksp_mu,ierr)
  call KSPSetUp(ksp_mu,ierr)
  call KSPSolve(ksp_mu,rhs_vec,mus_vec,ierr)

  call KSPGetConvergedReason(ksp_mu,mu_converged_reason,ierr)

  if (mu_converged_reason.gt.0) then

     call VecGetArrayF90(mus_vec,point_mu_vec,ierr)

     do z = 1,(psz+(2*ghost_width)-2)
        do y = 1,psy
           do x = 1,psx
              linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

              newmu(x,y,z+1) =  point_mu_vec(linindex)

           end do
        end do
     end do

     call VecRestoreArrayF90(mus_vec,point_mu_vec,ierr)

  else  ! if mu_converged_reason < 0

     newmu = mu

  end if


  call VecDestroy(mus_vec,ierr)
  call VecDestroy(rhs_vec,ierr)
  call MatDestroy(lhs_mat,ierr)
  call KSPDestroy(ksp_mu,ierr)


  !! Apply boundary conditions to chemical-potential-field update
  if (rank.eq.0) then
     do x = 1,psx
        do y = 1,psy
           newmu(x,y,1+ghost_width) = newmu(x,y,2+ghost_width)
        end do
     end do
  elseif(rank.eq.procs-1) then
     do x = 1,psx
        do y = 1,psy
           newmu(x,y,psz+ghost_width) = mu(x,y,psz+ghost_width)       
        end do
     end do
  end if

  !! Normalize chemical potential in bulk phases
  do x = 1,psx
     do y = 1,psy
        do z = 2,psz+(2*ghost_width)-1
           if (env(x,y,z).gt.0.97) then
              newmu(x,y,z) = avg_mu_env
           end if
        end do
     end do
  end do


!!! EVERYTHING BELOW IS IN nm/s

  sulf_rate_gas_met = 10**((0.00473*T)-5.645+((avg_mu_env+63562)/(R*T))) !! Ref = Assessing Corrosion in Oil Refining and Petrochemical Processing, Materials Research, Vol 7, No 1, pp. 163-173, 2004
  sulf_rate_gas_met = max(sulf_rate_gas_met,0.0d0) 

  sulf_rate_gas_pht = exp(-(11766/T)-0.6478)*1E9 !! Ref = Mechanisms Governing the Growth, Reactivity and Stability of Iron Sulfides, Ph.D Thesis William Herbert, MIT
  sulf_rate_gas_pht = max(sulf_rate_gas_pht,0.0d0) 

  sulf_rate_gas_pyr = 7.45E8 * exp(-(98400/(R*T))) !! Ref = Kinetics of sulfidation of chalcopyrite with gaseous sulfur. Padilla R. et. al., Met Trans B, Vol 34B, Feb 2003, 61-68, 
  sulf_rate_gas_pyr = max(sulf_rate_gas_pyr,0.0d0) 

  sulf_rate_liq_met = 0.0666*(1 + ((avg_mu_env+63562)/(R*T))) !! Ref = Kinetics of iron sulfiede and mixed iron sulfide/carbonate scale precipitation in CO2/H2S corrosion, Corrosion 2006, Paper 06644
  sulf_rate_liq_met = max(sulf_rate_liq_met,0.0d0) 

  sulf_rate_liq_mkw = 0.1332*(1 + (2*(avg_mu_env+63562)/(R*T))) !! Ref = Mechanistic model of H2S corrosion of mild steel
  sulf_rate_liq_mkw = max(sulf_rate_liq_mkw,0.0d0) 

  sulf_rate_liq_pht = 0.01372 + 0.04356*(exp(avg_mu_env/(R*T))) !! Ref = Corrosion, January 1990, Vol. 46, No. 1, pp. 66-74
  sulf_rate_liq_pht = max(sulf_rate_liq_pht,0.0d0) 

  sulf_rate_liq_pyr = 0.003543 !! Ref = Crystal growth of pyrite in Aqueous solutions. Inhibition by organophosphorous compounds, Harmandas NG. et. al., Langmuir 14, 1250-1255, 1998.
  sulf_rate_liq_pyr = max(sulf_rate_liq_pyr,0.0d0) 

  rho_pht = 52275.0d0
  rho_met = 0.0015d0*140401

  if (include_dissolve) then
     sulfidation_rate = sulf_rate_liq_pht
  else
     sulfidation_rate = sulf_rate_gas_met
  end if

  sulfidation_rate = sulfidation_rate*1E-1

  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           if ((env(x,y,z) .lt. 5.0E-1).and.(env(x,y,z+1) .gt. 5.0E-1)) then

              Chi = max(min(met(x,y,z-1)-0.005d0,1.0d0),0.0d0)*drho_dmu_met
              Chi = Chi + max(min(mkw(x,y,z-1)-0.005d0,1.0d0),0.0d0)*drho_dmu_mkw
              Chi = Chi + max(min(pht(x,y,z-1)-0.005d0,1.0d0),0.0d0)*drho_dmu_pht
              Chi = Chi + max(min(pyr(x,y,z-1)-0.005d0,1.0d0),0.0d0)*drho_dmu_pyr
              Chi = Chi + max(min(env(x,y,z-1)-0.005d0,1.0d0),0.0d0)*drho_dmu_env


              newmu(x,y,z) = mu(x,y,z-1) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt)/(dpf)) - max((dt*D(x,y,z-2)*(mu(x,y,z-1)-mu(x,y,z-2))/dpf),0.0d0) - (((dpht_dt(x,y,z-1)*rho_pht) + (dmet_dt(x,y,z-1)*rho_met) + (dmkw_dt(x,y,z-1)*rho_mkw) + (denv_dt(x,y,z-1)*rho_env) + (dpyr_dt(x,y,z-1)*rho_pyr))/Chi)
              newmu(x,y,z) = min(newmu(x,y,z),avg_mu_env)

              newmu(x,y,z-1) = mu(x,y,z-1) + ((((rho_pht-rho_met)/drho_dmu_pht)*sulfidation_rate*dt)/(dpf)) - max((dt*D(x,y,z-2)*(mu(x,y,z-1)-mu(x,y,z-2))/dpf),0.0d0) - (((dpht_dt(x,y,z-1)*rho_pht) + (dmet_dt(x,y,z-1)*rho_met) + (dmkw_dt(x,y,z-1)*rho_mkw) + (denv_dt(x,y,z-1)*rho_env) + (dpyr_dt(x,y,z-1)*rho_pyr))/Chi)
              newmu(x,y,z-1) = min(newmu(x,y,z-1),avg_mu_env)

              exit
           end if
        end do
     end do
  end do


 
  do x = 1,psx
     do y = 1,psy
        do z = 1+ghost_width,psz+ghost_width
           mu(x,y,z) = newmu(x,y,z)
        end do
     end do
  end do



end subroutine musolve
