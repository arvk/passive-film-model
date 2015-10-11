subroutine diffusivities
  use commondata
  use fields
  use thermo_constants
  use diffusion_constants
  implicit none
  !! **Diffusivity of \(Fe\) and \(S\) atoms in different \(FeS\) phases**


  !!---

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!

  !!#### Sulfur diffusivity in \(\frac{m^2}{s}\)

  !! $$D_{Fe}^{metal} = 7.87 \times 10^{-7} \times e^{-96500 \times \frac{0.60}{R T}}$$
  !! Reference : Physical Review B 80, 144111-4 2009
  D_Fe_met = 7.87E-7 * exp((96500*(0.0-0.60d0))/(R*T))

  !! $$D_{S}^{metal} = 1.78 \times e^{-\frac{205000}{R T}}$$
  !! Reference : Scripta Metallurgica, Col 16, pp 537-540, 1982. (Fillastre C. et. al. - Diffusion of sulfur in iron-chromium alloys)
  D_S_met = 1.78d0*exp((0.0d0-205000)/(R*T))

  !!----------------------------------------

  !! $$D_{Fe}^{mackinawite} = 3.37 \times 10^{-16} \times e^{-\frac{15500}{R T}}$$
  !! Reference : Simulation Of Solid-State Growth Of Iron Sulfides In Sour Corrosion Conditions, Jamel Amri and Jon Kvarekvål, NACE 2011, Paper 11078
  D_Fe_mkw = 3.37E-16 * exp((0.0d0-15500)/(R*T))

  !! $$D_{S}^{mackinawite} = D_{Fe}^{mackinawite} $$
  !! Reference : Assumption
  D_S_mkw = D_Fe_mkw

  !!----------------------------------------

  !! $$D_{Fe}^{pyrrhotite} = 10^{-\frac{7056}{R T}-2.679}$$
  !! Reference :  Mechanisms Governing the Growth, Reactivity and Stability of Iron Sulfides, Ph.D Thesis William Herbert, MIT
  if (T.lt.(273.15+315)) then
     D_Fe_pht = 10**((0-7056/T)-2.679)
  else
     D_Fe_pht = 10**((0-4220/T)-7.311)
  end if

  !! $$D_{S}^{pyrrhotite} = D_{Fe}^{pyrrhotite}/10$$
  !! Reference : Assumption of one order of magnitude lesser
  D_S_pht = D_Fe_pht/10

  !!----------------------------------------

  !! $$D_{Fe}^{pyrite} = 1.75 \times 10^{-14} \times e^{-\frac{132100}{R T}}$$
  !! Reference : Geochimica et Cosmochimica Acta 73 (2009) 4792–4802
  D_S_pyr = 1.75E-14 * exp((0.0d0-132100.0d0)/(R*T))

  !! $$D_{S}^{pyrite} = 2.5 \times 10^{-16} \times e^{-\frac{41840}{R T}}$$
  !! Reference: Metallurgical Transactions B Volume 6B, JUNE 1975-333
  D_Fe_pyr = 2.5E-16 * exp((0.0d0-41840.0d0)/(R*T))

  !!----------------------------------------

  !! $$D_{Fe}^{environment} = 7.19 \times 10^{-10}$$
  !! Reference : http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
  D_Fe_env =7.19E-10

  !! $$D_{S}^{environment} = 1.73 \times 10^{-9}$$
  !! Reference : J. Chem. Eng. Data 1994, 39, 330-332
  !! Reference : http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
  D_S_env = 1.73E-09

  !!----------------------------------------


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!


  !!#### Hydrogen diffusivity in \(\frac{m^2}{s}\)

  !! $$D_{H}^{environment} = 9.31 \times 10^{-9}$$
  !! Reference : http://mail.sci.ccny.cuny.edu/~pzhang/EAS44600/EAS446lec16.pdf
   D_H_env = 9.31E-9

  !!----------------------------------------

   !! $$D_{H}^{metal} = 1.57 \times 10^{-7} \times e^{-\frac{8492}{R T}}$$
   !! Reference : Physical Review B 70 , 064102-7 (2004)
   D_H_met = 1.5E-7 * exp(0-(8492/(R*T)))

  !!----------------------------------------

   !! $$D_{H}^{pyrrhotite} = 1.75 \times 10^{-9}$$
   !! Reference : Int. J. Electrochem. Sci.,8, (2013) 2880-2891
   D_H_pht = 1.75E-9

  !!----------------------------------------

   !! $$D_{H}^{pyrite} = D_{H}^{metal}$$
   !! Reference : Assumed. Cannot find a suitable reference
   D_H_pyr = D_H_met

end subroutine diffusivities
