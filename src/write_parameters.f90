subroutine write_parameters()
  use commondata
  use fields
  implicit none
  !!####Write-out current simulation parameters to the PARAMS.out file. This can be used to restart a simulation
  !!---

  call system("rm -rf PARAMS.out")
  open(unit = 8000, file = "PARAMS.out", status = 'new')

  write(8000,'(A10, A1)') "RESTART = ", isrestart

  write(8000,'(A9,I4.4,1X,I4.4,1X,I4.4)') "PFSIZE = ", psx_g, psy_g, psz_g

  write(8000,'(A7,I4.4)') "TEMP = ", T

  write(8000,'(A5,F5.2)') "PH = ", pH_in

  write(8000,'(A9,I8.8)') "SIMLEN = ", nomc

  if (include_dissolve) then
     write(8000,'(A12)') "DISSOLVE = Y"
  else
     write(8000,'(A12)') "DISSOLVE = N"
  end if

  if (include_electro) then
     write(8000,'(A11)') "ELECTRO = Y"
     write(8000,'(A10,F6.2)') "METPOTL = ",metal_potential
  else
     write(8000,'(A11)') "ELECTRO = N"
  end if

  write(8000,'(A12,I5.5)') "SNAPSHOTS = ", num_images

  close(8000)

end subroutine write_parameters
