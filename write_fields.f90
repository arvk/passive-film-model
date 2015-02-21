subroutine write_fields(iter)
  use commondata
  use fields
  use kmc_data
  implicit none

  integer, intent(in) :: iter  
  integer :: x, y, z
  character*5 :: img_id

  call system("rm -rf PHT.out ENV.out MET.out PYR.out MUS.out OPYR.out")


  !! Write phase fraction fields
  open (unit=101, file="MET.out", status="new")
  write(101,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(101,*)  met_g(x,y,z)
        end do
     end do
  end do
  close(101)  

  open (unit=102, file="MKW.out", status="new")
  write(102,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(102,*)  mkw_g(x,y,z)
        end do
     end do
  end do
  close(102)  

  open (unit=103, file="PHT.out", status="new")
  write(103,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(103,*)  pht_g(x,y,z)
        end do
     end do
  end do
  close(103)  

  open (unit=104, file="PYR.out", status="new")
  write(104,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(104,*)  pyr_g(x,y,z)
        end do
     end do
  end do
  close(104)  

  open (unit=105, file="ENV.out", status="new")
  write(105,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(105,*)  env_g(x,y,z)
        end do
     end do
  end do
  close(105)  



  !! Write chemical potential fields
  open (unit=200, file="MUS.out", status="new")
  write(200,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(200,*)  mu_g(x,y,z)
        end do
     end do
  end do
  close(200)  



  !! Write orientation fields
  open (unit=304, file="OPYR.out", status="new")
  write(304,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(304,*)  opyr_g(x,y,z)
        end do
     end do
  end do
  close(304)  




  write(img_id,'(I5.5)') iter/((nomc/noimg)-1)
  call system("cp MET.out MET_"//img_id//".out")
  call system("cp PHT.out PHT_"//img_id//".out")
  call system("cp PYR.out PYR_"//img_id//".out")
  call system("cp ENV.out ENV_"//img_id//".out")

  call system("cp MUS.out MUS_"//img_id//".out")

  call system("cp OPYR.out OPYR_"//img_id//".out")

end subroutine write_fields

