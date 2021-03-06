subroutine initialize_geometry(simstate)
  use commondata
  use fields
  use thermo_constants
  implicit none
#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscmat.h>
#include <finclude/petscpc.h>
#include <finclude/petscksp.h>
#include <finclude/petscsnes.h>
#include <finclude/petscdm.h>
#include <finclude/petscdmda.h>
#include <finclude/petscdmda.h90>
  !! **Initialize the simulation system containing a metal surface in contact with the given environment**

  PetscErrorCode :: ierr !! MPI error flag
  PetscInt :: x,y,z                               !! Coordinates inside the simulation system
  PetscInt :: met_z_end                           !! Location of boundary between metal and environment
  PetscScalar :: avg_mu_met                           !! Chemical potential in the metal region
  PetscScalar, parameter :: infinitesimal = 0.00001d0 !! A hard-coded 'small' number
  PetscScalar, pointer :: statepointer(:,:,:,:) !! Pointer array referenced to individual gridpoints inside the simulation cell
  PetscRandom :: random_orientation_context     !! Context to seed and generate random numbers for the orientation field
  PetscReal :: random_number                    !! Pseudo random number generated from a PETSc context
  type(context), intent(inout) :: simstate      !! Field variables stored in PETSc vectors and DMDA objects

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  met_z_end = 15*(psz_g/16)

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end         ! Initialize metal region
           met_g(x,y,z) = 1.0d0-infinitesimal
           mkw_g(x,y,z) = 0.0d0+infinitesimal
           pht_g(x,y,z) = 0.0d0+infinitesimal
           pyr_g(x,y,z) = 0.0d0+infinitesimal
           env_g(x,y,z) = 0.0d0+infinitesimal
        end do

        do z = met_z_end+1,psz_g   ! Initialize environmental region
           met_g(x,y,z) = 0.0d0+infinitesimal
           mkw_g(x,y,z) = 0.0d0+infinitesimal
           pht_g(x,y,z) = 0.0d0+infinitesimal
           pyr_g(x,y,z) = 0.0d0+infinitesimal
           env_g(x,y,z) = 1.0d0-infinitesimal
        end do

     end do
  end do

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  avg_mu_met = mus_met_mkw_eqb - (R*T*0.5d0)
  avg_mu_env = mus_mkw_pht_eqb - (R*T*2.5d0)

  do x = 1,psx_g
     do y = 1,psy_g

        do z = 1,met_z_end
           mu_g(x,y,z) = avg_mu_met    ! Initialize metal region
        end do

        do z = met_z_end+1,psz_g
           mu_g(x,y,z) = avg_mu_env    ! Initialize environmental region
        end do

     end do
  end do

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  elpot_g = metal_potential            ! Initialize global electrical potential

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  call PetscRandomCreate(MPI_COMM_WORLD,random_orientation_context,ierr)
  call PetscRandomSetType(random_orientation_context,PETSCRAND,ierr)

  call DMDAVecGetArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)
  do x = simstate%startx , simstate%startx + simstate%widthx-1
     do y = simstate%starty , simstate%starty + simstate%widthy-1
        do z = simstate%startz , simstate%startz + simstate%widthz-1
           call PetscRandomGetValueReal(random_orientation_context,random_number,ierr)
           statepointer(nmet,x,y,z) = met_g(x+1,y+1,z+1)
           statepointer(nmkw,x,y,z) = mkw_g(x+1,y+1,z+1)
           statepointer(npht,x,y,z) = pht_g(x+1,y+1,z+1)
           statepointer(npyr,x,y,z) = pyr_g(x+1,y+1,z+1)
           statepointer(nenv,x,y,z) = env_g(x+1,y+1,z+1)
           statepointer(nmus,x,y,z) = mu_g(x+1,y+1,z+1)
           statepointer(npH,x,y,z) = 10**(0.0d0-pH_in)
           statepointer(nang,x,y,z) = random_number * (3.14159265d0/2.0d0)   ! Populate the global orientation field with random numbers
           statepointer(npot,x,y,z) = elpot_g(x+1,y+1,z+1)
           statepointer(nvoi,x,y,z) = 1.0d0
        end do
     end do
  end do
  call DMDAVecRestoreArrayF90(simstate%lattval,simstate%slice,statepointer,ierr)

  call PetscRandomDestroy(random_orientation_context,ierr)

  call DMDASetFieldName(simstate%lattval,nmet,"Metal",ierr)
  call DMDASetFieldName(simstate%lattval,nmkw,"Mackinawite",ierr)
  call DMDASetFieldName(simstate%lattval,npht,"Pyrrhotite",ierr)
  call DMDASetFieldName(simstate%lattval,npyr,"Pyrite",ierr)
  call DMDASetFieldName(simstate%lattval,nenv,"Environment",ierr)
  call DMDASetFieldName(simstate%lattval,nmus,"Sulfur chemical potential",ierr)
  call DMDASetFieldName(simstate%lattval,npH,"H+ concentration",ierr)
  call DMDASetFieldName(simstate%lattval,nang,"Grain orientation",ierr)
  call DMDASetFieldName(simstate%lattval,npot,"Electric potential",ierr)
  call DMDASetFieldName(simstate%lattval,nvoi,"Voids",ierr)

end subroutine initialize_geometry

