subroutine diffusivities
  use commondata
  use fields
  use diffusion_constants

  D_Fe_pht = 10**((0-7056/T)-3.679)    ! Ref = William's data
  D_S_pht = D_Fe_pht/10                ! Ref = Assumption of one order of magnitude lesser

  D_Fe_met = 7.87E-9 * exp((96500*(0.0-0.60d0))/(R*T))   ! Ref = PHYSICAL REVIEW B 80, 144111-4 2009
  D_S_met = 0                                            ! Assumed to be negligible

  D_S_env = 1.70E-15  ! Ref = J. Chem. Eng. Data 1994, 39, 330-332 and http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
  D_Fe_env = 7.0E-10 ! Ref = http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf

  D_S_pyr = 1.75E-14 * exp((0.0d0-132100.0d0)/(R*T))   ! Ref = Geochimica et Cosmochimica Acta 73 (2009) 4792–4802
  D_Fe_pyr = 2.5E-16 * exp((0.0d0-41840.0d0)/(R*T))    ! Ref = METALLURGICAL TRANSACTIONS B VOLUME 6B, JUNE 1975-333

end subroutine diffusivities
