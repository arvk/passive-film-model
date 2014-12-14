subroutine seed_prng
  use commondata
  implicit none

  integer :: n,istat
  integer, dimension(8) :: datetime

   open(89, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(89) seed
     close(89)
  else
     call date_and_time(values=datetime)
     seed(n) = datetime(8); seed(1) = datetime(8)*datetime(7)*datetime(6)
  end if

  call random_seed(put=seed)

end subroutine seed_prng
