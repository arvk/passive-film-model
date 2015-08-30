subroutine write_fields(iter)
  use commondata
  use fields
  use kmc_data
  implicit none

  integer :: x, y, z          ! Index along x-, y- and z-directions (Loop)
  integer, intent(in) :: iter ! Iteration number
  character*5 :: img_id       ! Index of current image (derived from current iteration number)

  call system("rm -rf MET.out MKW.out PHT.out PYR.out ENV.out MUS.out OPYR.out POT.out PH.out")

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Write phase-fraction fields; Units == 1XX
  open (unit=101, file="MET.out", status="new")
  write(101,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(101,'(F9.6)')  met_g(x,y,z)         ! Write metal phase fraction
        end do
     end do
  end do
  write(101,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(101)

  open (unit=102, file="MKW.out", status="new")
  write(102,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(102,'(F9.6)')  mkw_g(x,y,z)         ! Write mackinawite phase fraction
        end do
     end do
  end do
  write(102,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(102)

  open (unit=103, file="PHT.out", status="new")
  write(103,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(103,'(F9.6)')  pht_g(x,y,z)         ! Write pyrrhotite phase fraction
        end do
     end do
  end do
  write(103,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(103)

  open (unit=104, file="PYR.out", status="new")
  write(104,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(104,'(F9.6)')  pyr_g(x,y,z)         ! Write pyrite phase fraction
        end do
     end do
  end do
  write(104,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(104)

  open (unit=105, file="ENV.out", status="new")
  write(105,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(105,'(F9.6)')  env_g(x,y,z)         ! Write environment phase fraction
        end do
     end do
  end do
  write(105,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(105)

  write(6,'(A,F9.0,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2)') " INFO: TIME= ", iter*dt, " s. MET= ", sum(met_g), " MKW= ", sum(mkw_g), " PHT= ", sum(pht_g), " PYR= ", sum(pyr_g), " ENV= ", sum(env_g)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Write chemical-potential field; Units == 2XX
  open (unit=200, file="MUS.out", status="new")
  write(200,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(200,'(F9.0)')  mu_g(x,y,z)          ! Write chemical potential field
        end do
     end do
  end do
  write(200,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(200)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Write orientation field; Units == 3XX
  open (unit=304, file="OPYR.out", status="new")
  write(304,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(304,'(F9.6)')  opyr_g(x,y,z)        ! Write orientation field
        end do
     end do
  end do
  write(304,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(304)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Write electri potential field; Units == 4XX
  open (unit=401, file="POT.out", status="new")
  write(401,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(401,'(F6.3)')  elpot_g(x,y,z)       ! Write electric potential field
        end do
     end do
  end do
  write(401,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(401)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!

  ! Write pH field; Units == 5XX
  open (unit=500, file="PH.out", status="new")
  write(500,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(500,'(F9.6)')  pH_g(x,y,z)          ! Write pH field
        end do
     end do
  end do
  write(500,'(A,F9.0,A)') "# TIME: ", iter*dt, " s"
  close(500)

!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!!!!


  write(img_id,'(I5.5)') iter/max(floor(real(nomc/noimg)),1)
  call system("cp MET.out MET_"//img_id//".out")
  call system("cp MKW.out MKW_"//img_id//".out")
  call system("cp PHT.out PHT_"//img_id//".out")
  call system("cp PYR.out PYR_"//img_id//".out")
  call system("cp ENV.out ENV_"//img_id//".out")

  call system("cp MUS.out MUS_"//img_id//".out")

  call system("cp PH.out PH_"//img_id//".out")

  call system("cp OPYR.out OPYR_"//img_id//".out")

  call system("cp POT.out POT_"//img_id//".out")

end subroutine write_fields

