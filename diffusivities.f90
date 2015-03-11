subroutine diffusivities
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  
! in m^2/s

  D_Fe_met = 7.87E-7 * exp((96500*(0.0-0.60d0))/(R*T))   ! Ref = PHYSICAL REVIEW B 80, 144111-4 2009
  D_S_met = 0                                            ! Assumed to be negligible

  D_Fe_mkw = 3.37E-16 * exp((0.0d0-15500)/(R*T))         ! Ref = SIMULATION OF SOLID-STATE GROWTH OF IRON SULFIDES IN SOUR CORROSION CONDITIONS, Jamel Amri and Jon Kvarekvål, NACE 2011, Paper 11078
  D_S_mkw = D_Fe_mkw                                     ! Ref = Assumptions. To be later taken from other models

  D_Fe_pht = 10**((0-7056/T)-2.679)    ! Ref = William's data
  D_S_pht = D_Fe_pht/10                ! Ref = Assumption of one order of magnitude lesser

  D_S_pyr = 1.75E-14 * exp((0.0d0-132100.0d0)/(R*T))   ! Ref = Geochimica et Cosmochimica Acta 73 (2009) 4792–4802
  D_Fe_pyr = 2.5E-16 * exp((0.0d0-41840.0d0)/(R*T))    ! Ref = METALLURGICAL TRANSACTIONS B VOLUME 6B, JUNE 1975-333

  D_S_env = 1.73E-09  ! Ref = J. Chem. Eng. Data 1994, 39, 330-332 and http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
  D_Fe_env =7.19E-10 ! Ref = http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf


!! Hydrogen diffusivity

!  D_H_env = 9.31E-9  ! Ref = http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
!  D_H_met = 1.5E-7 * exp(0-(8492/(R*T))) ! Ref = PHYSICAL REVIEW B 70 , 064102-7 (2004)
!  D_H_pht = 1.75E-9 ! Ref = Int. J. Electrochem. Sci.,8, (2013) 2880-2891
!  D_H_pyr = D_H_met ! Assumed. Cannot find reference


! D_Fe_met = 7.87E-7 * exp((96500*(0.0-0.60d0))/(R*T))   ! Ref = PHYSICAL REVIEW B 80, 144111-4 2009
! D_S_env = 1.70E-9  ! Ref = J. Chem. Eng. Data 1994, 39, 330-332 and http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf

end subroutine diffusivities
