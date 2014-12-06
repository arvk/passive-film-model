subroutine draw_progress_bar(iter,nomc,dt)
  implicit none

  integer, intent(in) :: iter, nomc
  real*8, intent(in) :: dt

  character(len=68)::bar="????????s, ???% |                                                  |"

  integer :: i            ! Loop variable


  !! Calculate the timecount and percentage
  write(unit=bar(12:14),fmt="(i3)") int((100.0d0*iter)/nomc)
  write(unit=bar(1:8),fmt="(i8)") int(iter*dt)

  !! Define length of progress bar
  do i=1, int((50.0d0*iter)/nomc)
     bar(17+i:17+i)="="
  enddo

  !! (Over)write the progress bar
  write(unit=6,fmt="(a1,a1,a68)", advance="no") '+',char(13), bar

  return

end subroutine draw_progress_bar

