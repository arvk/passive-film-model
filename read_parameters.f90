subroutine read_parameters()
  use commondata
  use fields
  use kmc_data
  implicit none

  integer :: error_temp, error_nomc, error_ph, error_dissolve, error_electro, error_metpotl, error_snapshots
  character*1 :: isdissolve
  character*1 :: useelectro

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= RESTART | Question= Is the simulation being started from scratch? | Output_var= isrestart !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *RESTART' | sed 's/RESTART//g'| sed 's/=//g' > .isrestart.readin")
  open(unit = 7000, file = ".isrestart.readin", status = 'old') ; read(7000,*) isrestart ; close(7000)
  call system ("rm .isrestart.readin")

  if ((isrestart .eq. 'Y') .or. (isrestart .eq. 'y')) then ! If it IS a restarted run
     write (6,*) 'INFO: Simulation restarted from a previous run.'
  else                                                     ! If it IS NOT a restarted run
     write(6,*) 'INFO: Simulation started from scratch.'
  end if
  !! FINISH | Tag= RESTART !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= PFSIZE | Question= What is the size of the phase-field simulation box? | Output_var= psx_g, psy_g, psz_g  !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *PFSIZE' | sed 's/PFSIZE//g'| sed 's/=//g' > .pfsize.readin")
  open(unit = 7101, file = ".pfsize.readin", status = 'old') ; read(7101,*) psx_g, psy_g, psz_g ; close(7101)
  call system ("rm .pfsize.readin")

  if ((psx_g .le. 0) .or. (psy_g .le. 0) .or. (psz_g .le. 0)) then ! If negative values are read in
     stop "ERROR: Phase field grid sizes are non-positive. Exiting."
  end if

  psx = psx_g ; psy = psy_g ; psz = psz_g/procs
  ksx_g = psx_g*kg_scale ; ksy_g = psy_g*kg_scale
  ksx = ksx_g ; ksy = ksy_g/procs
  !! FINISH | Tag= PFSIZE !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= TEMP | Question= What is the temperature of the simulation box? | Output_var= T !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *TEMP' | sed 's/TEMP//g'| sed 's/=//g' > .temp.readin")
  open(unit = 7102, file = ".temp.readin", status = 'old');   read(7102,*,IOSTAT=error_temp) T ; close(7102)
  call system ("rm .temp.readin")

  if (error_temp .lt. 0) then ! If no tempertaure values are read in
     T = 573
     write(6,*) "WARNING: No temperature input. Simulation will continue at 573K."
  end if

  if (T .le. 0) then ! If negative temperatures are read in
     T = 573
     write(6,*) "WARNING: Error in reading temperature. Simulation will continue at 573K."
  end if
  !! FINISH | Tag= TEMP !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= PH | What is the pH of the simulation box? | Output_var= pH_in !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *PH' | sed 's/PH//g'| sed 's/=//g' > .ph.readin")
  open(unit = 7102, file = ".ph.readin", status = 'old');   read(7102,*,IOSTAT=error_ph) pH_in ; close(7102)
  call system ("rm .ph.readin")

  if (error_ph .lt. 0) then ! If no pH values are read in
     pH_in = 3.0d0
     write(6,*) "WARNING: No pH input. Simulation will continue at a pH of 3."
  end if
  !! FINISH | Tag= PH !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= SIMLEN | Question= What is the length of the simulation? | Output_var= nomc !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *SIMLEN' | sed 's/SIMLEN//g'| sed 's/=//g' > .nomc.readin")
  open(unit = 7103, file = ".nomc.readin", status = 'old');   read(7103,*,IOSTAT=error_nomc) nomc ; close(7103)
  call system ("rm .nomc.readin")

  if (error_nomc .lt. 0) then ! If no simulation length is read in
     nomc = 1000000
     write(6,*) "WARNING: No SIMLEN specified. Simulation will run for 10^6 steps."
  end if

  if (nomc .le. 0) then ! If simulation length is negative
     nomc = 1000000
     write(6,*) "WARNING: Error in reading SIMLEN. Simulation will run for 10^6 steps."
  end if
  !! FINISH | Tag= SIMLEN !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= DISSOLVE | Question= Should dissolution be included in the simulation? | Output_var= include_dissolve !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *DISSOLVE' | sed 's/DISSOLVE//g'| sed 's/=//g' > .isdissolve.readin")
  open(unit = 7000, file = ".isdissolve.readin", status = 'old') ; read(7000,*,IOSTAT=error_dissolve) isdissolve ; close(7000)
  call system ("rm .isdissolve.readin")

  if (error_dissolve .lt. 0) then ! If dissolution tag is not read in

     include_dissolve = .FALSE.
     write(6,*) "WARNING: No dissolution flag set. No film dissolution."

  else                            ! If dissolution tag is correctly read in

     if ((isdissolve .eq. 'Y') .or. (isdissolve .eq. 'y')) then ! If dissolution SHOULD be included in the simulation
        include_dissolve = .TRUE.
        write (6,*) 'INFO: Liquid environment. Film dissolution included.'
     else                                                       ! If dissolution SHOULD NOT be included in the simulation
        include_dissolve = .FALSE.
        write(6,*) 'INFO: Gaseous environment. No film dissolution.'
     end if

  end if
  !! FINISH | Tag= SIMLEN !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= ELECTRO | Question= Should the electrochemistry module be included? | Output_var= include_electro !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *ELECTRO' | sed 's/ELECTRO//g'| sed 's/=//g' > .useelectro.readin")
  open(unit = 7000, file = ".useelectro.readin", status = 'old') ; read(7000,*,IOSTAT=error_electro) useelectro ; close(7000)
  call system ("rm .useelectro.readin")


  if (error_electro .lt. 0) then
     include_electro = .FALSE.
     write(6,*) "WARNING: Electrochemistry flag not set. Electrochemistry module not included."

  else

     if ((useelectro .eq. 'Y') .or. (useelectro .eq. 'y')) then ! If dissolution SHOULD be included in the simulation
        include_electro = .TRUE.
        write (6,*) 'INFO: Electrochemistry module included.'

        !!----Start reading Metal Potential-----!!
        call system("cat param.in | sed 's/!.*//g' | grep '^ *METPOTL' | sed 's/METPOTL//g'| sed 's/=//g' > .metpotl.readin")
        open(unit = 7000, file = ".metpotl.readin", status = 'old') ; read(7000,*,IOSTAT=error_metpotl) metal_potential ; close(7000)
        call system ("rm .metpotl.readin")

        if (error_metpotl .lt. 0) then ! If no metal potential is read in
           metal_potential = -0.50d0
           write(6,*) "WARNING: Metal potential not set. -0.5V SHE assumed."
        end if
        !!----End reading Metal Potential-----!!

     else                                                      ! If dissolution SHOULD NOT be included in the simulation
        include_electro = .FALSE.
        write(6,*) 'INFO: Electrochemistry module not included.'
     end if

  end if
  !! FINISH | Tag= ELECTRO !!

!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!

  !! START | Tag= SNAPSHOTS | Question= How many output snapshots? | Output_var= noimg !!
  call system("cat param.in | sed 's/!.*//g' | grep '^ *SNAPSHOTS' | sed 's/SNAPSHOTS//g'| sed 's/=//g' > .nosnapshots.readin")
  open(unit = 7109, file = ".nosnapshots.readin", status = 'old');   read(7109,*,IOSTAT=error_snapshots) noimg ; close(7109)
  call system ("rm .nosnapshots.readin")

  if (error_snapshots .lt. 0) then ! If no snapshots count is read in
     noimg = min(nomc,100)
     write(6,*) "WARNING: No SNAPSHOTS specified. 100 output files will be created."
  end if

  if (noimg .le. 0) then ! If number of snapshots is negative
     noimg = min(nomc,100)
     write(6,*) "WARNING: Error in reading SNAPSHOTS. 100 output files will be created."
  end if
  !! FINISH | Tag= SNAPSHOTS !!

end subroutine read_parameters
