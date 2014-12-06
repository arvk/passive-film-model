subroutine write_fields(iter)
  use commondata
  implicit none

  integer, intent(in) :: iter  
  integer :: x, y, z
  character*5 :: img_id

  call system("rm -rf PHT.out ENV.out MET.out MUS.out Hplus.out ALL.out")

  open (unit=101, file="PHT.out", status="new")
  write(101,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(101,*)  pht_g(x,y,z)
        end do
     end do
  end do
  close(101)  

  open (unit=102, file="ENV.out", status="new")
  write(102,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(102,*)  env_g(x,y,z)
        end do
     end do
  end do
  close(102)  

  open (unit=103, file="MET.out", status="new")
  write(103,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g
     do y = 1,psy_g
        do z = 1,psz_g
           write(103,*)  met_g(x,y,z)
        end do
     end do
  end do
  close(103)  

  open (unit=105, file="MUS.out", status="new")
  write(105,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(105,*)  mu_g(x,y,z)
        end do
     end do
  end do
  close(105)  

  open (unit=106, file="Hplus.out", status="new")
  write(106,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(106,*)  ph_g(x,y,z)
        end do
     end do
  end do
  close(106)  


  open (unit=107, file="ALL.out", status="new")
  write(107,*) "#", psx_g, psy_g, psz_g, ksx_g, ksy_g
  do x = 1,psx_g 
     do y = 1,psy_g
        do z = 1,psz_g
           write(107,*)  (met_g(x,y,z)*1.0d0)+(pht_g(x,y,z)*2.0d0)+(env_g(x,y,z)*4.0d0)
        end do
     end do
  end do
  close(107)  


  write(img_id,'(I5.5)') iter/((nomc/noimg)-1)
  call system("cp PHT.out PHT_"//img_id//".out")
  call system("cp ENV.out ENV_"//img_id//".out")
  call system("cp MET.out MET_"//img_id//".out")
  call system("cp MUS.out MUS_"//img_id//".out")
  call system("cp Hplus.out Hplus_"//img_id//".out")
  call system("cp ALL.out ALL_"//img_id//".out")

end subroutine write_fields

