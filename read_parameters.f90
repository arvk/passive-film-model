subroutine read_parameters()
  use commondata
  use fields
  use kmc_data
  implicit none

  integer :: error_temp, error_nomc, error_ph, error_dissolve, error_electro, error_metpotl
  character*1 :: isdissolve
  character*1 :: useelectro

  !! Is the simulation being started from scratch?
  call system("cat param.in | grep ^RESTART | sed 's/RESTART//g'| sed 's/=//g' > .isrestart.readin")
  open(unit = 7000, file = ".isrestart.readin", status = 'old') ; read(7000,*) isrestart ; close(7000)
  call system ("rm .isrestart.readin")


  !! If it IS a restarted run
  if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then
     write (6,*) 'Restarted run'

     !! If it IS NOT a restarted run
  else
     write(6,*) '-------------------------------------'
     write(6,*) 'INFO: Simulation started from scratch'
     write(6,*) '-------------------------------------'

  end if  !! End of isrestart loop

  !! What is the size of the phase-field simulation box?
  call system("cat param.in | grep ^PFSIZE | sed 's/PFSIZE//g'| sed 's/=//g' > .pfsize.readin")
  open(unit = 7101, file = ".pfsize.readin", status = 'old') ; read(7101,*) psx_g, psy_g, psz_g ; close(7101)
  call system ("rm .pfsize.readin")

  if ((psx_g .le. 0) .or. (psy_g .le. 0) .or. (psz_g .le. 0)) then
     stop "ERROR: Phase field grid sizes are non-positive. Exiting"
  end if

  psx = psx_g ; psy = psy_g ; psz = psz_g/procs
  ksx_g = psx_g*kg_scale ; ksy_g = psy_g*kg_scale
  ksx = ksx_g ; ksy = ksy_g/procs

  !! What is the temperature of the simulation box?
  call system("cat param.in | grep ^TEMP | sed 's/TEMP//g'| sed 's/=//g' > .temp.readin")
  open(unit = 7102, file = ".temp.readin", status = 'old');   read(7102,*,IOSTAT=error_temp) T ; close(7102)
  call system ("rm .temp.readin")

  if (error_temp .lt. 0) then
     T = 573
     write(6,*) "----------------------------------------------------------------"
     write(6,*) "WARNING: No temperature input. Simulation will continue at 573K."    
     write(6,*) "----------------------------------------------------------------"
  end if

  if (T .le. 0) then
     T = 573
     write(6,*) "WARNING: Error in reading temperature. Simulation will continue at 573K."    
  end if



  !! What is the pH of the simulation box?
  call system("cat param.in | grep ^PH | sed 's/PH//g'| sed 's/=//g' > .ph.readin")
  open(unit = 7102, file = ".ph.readin", status = 'old');   read(7102,*,IOSTAT=error_ph) pH_in ; close(7102)
  call system ("rm .ph.readin")

  if (error_ph .lt. 0) then
     pH_in = 3.0d0
     write(6,*) "------------------------------------------------------------"
     write(6,*) "WARNING: No pH input. Simulation will continue at a pH of 3."    
     write(6,*) "------------------------------------------------------------"
  end if


  !! What is the length of the simulation?
  call system("cat param.in | grep ^SIMLEN | sed 's/SIMLEN//g'| sed 's/=//g' > .nomc.readin")
  open(unit = 7103, file = ".nomc.readin", status = 'old');   read(7103,*,IOSTAT=error_nomc) nomc ; close(7103)
  call system ("rm .nomc.readin")

  if (error_nomc .lt. 0) then
     nomc = 1000000
     write(6,*) "-----------------------------------------------------------------"
     write(6,*) "WARNING: No SIMLEN specified. Simulation will run for 10^6 steps."    
     write(6,*) "-----------------------------------------------------------------"
  end if

  if (nomc .le. 0) then
     nomc = 1000000
     write(6,*) "---------------------------------------------------------------------"
     write(6,*) "WARNING: Error in reading SIMLEN. Simulation will run for 10^6 steps."
     write(6,*) "---------------------------------------------------------------------"   
  end if



  !! Should dissolution be included in the simulation?
  call system("cat param.in | grep ^DISSOLVE | sed 's/DISSOLVE//g'| sed 's/=//g' > .isdissolve.readin")
  open(unit = 7000, file = ".isdissolve.readin", status = 'old') ; read(7000,*,IOSTAT=error_dissolve) isdissolve ; close(7000)
  call system ("rm .isdissolve.readin")

  if (error_dissolve .lt. 0) then

     include_dissolve = .FALSE.
     write(6,*) "-----------------------------------------------------"
     write(6,*) "WARNING: No dissolution flag set.No film dissolution."
     write(6,*) "-----------------------------------------------------"

  else

     !! If dissolution SHOULD be included in the simulation
     if ((isdissolve .eq. 'Y') .or. (isdissolve .eq. 'y')) then
        include_dissolve = .TRUE.
        write (6,*) 'Liquid environment. Film dissolution included.'
        !! If dissolution SHOULD NOT be included in the simulation
     else
        include_dissolve = .FALSE.
        write(6,*) 'Gaseous environment. No film dissolution.'

     end if  !! End of isdissolve loop


  end if



  !! Should the electrochemistry module be included?
  call system("cat param.in | grep ^ELECTRO | sed 's/ELECTRO//g'| sed 's/=//g' > .useelectro.readin")
  open(unit = 7000, file = ".useelectro.readin", status = 'old') ; read(7000,*,IOSTAT=error_electro) useelectro ; close(7000)
  call system ("rm .useelectro.readin")


  if (error_electro .lt. 0) then
     include_electro = .FALSE.
     write(6,*) "-----------------------------------------------------------------------------"
     write(6,*) "WARNING: Electrochemistry flag not set. Electrochemistry module not included."
     write(6,*) "-----------------------------------------------------------------------------"

  else

     !! If dissolution SHOULD be included in the simulation
     if ((useelectro .eq. 'Y') .or. (useelectro .eq. 'y')) then
        include_electro = .TRUE.
        write (6,*) 'Electrochemistry module included.'


        call system("cat param.in | grep ^METPOTL | sed 's/METPOTL//g'| sed 's/=//g' > .metpotl.readin")
        open(unit = 7000, file = ".metpotl.readin", status = 'old') ; read(7000,*,IOSTAT=error_metpotl) metal_potential ; close(7000)
        call system ("rm .metpotl.readin")

        if (error_metpotl .lt. 0) then
           metal_potential = -0.50d0
           write(6,*) "----------------------------------------------------"
           write(6,*) "WARNING: Metal potential not set. -0.5V SHE assumed."
           write(6,*) "----------------------------------------------------"
        end if


        !! If dissolution SHOULD NOT be included in the simulation
     else
        include_electro = .FALSE.
        write(6,*) 'Electrochemistry module not included.'

     end if  !! End of useelectro loop

  end if

end subroutine read_parameters
