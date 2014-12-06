subroutine write_parameters()
  use commondata
  implicit none

  call system("rm -rf PARAMS.out")
  open(unit = 8000, file = "PARAMS.out", status = 'new') ; 
  write(8000,'(A10, A1)') "RESTART = ", isrestart
  write(8000,'(A10,I5,I5,I5)') "KMCSIZE = ", ksx_g, ksy_g
  write(8000,'(A9,I5,I5,I5)') "PFSIZE = ", psx_g, psy_g, psz_g
  write(8000,'(A7,I4)') "TEMP = ", T
  write(8000,'(A5,F6.2)') "pH = ", pH_in
  write(8000,'(A9,I8)') "SIMLEN = ", nomc
  close(8000)

end subroutine write_parameters
