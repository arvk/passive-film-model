subroutine kmcsolve(iter)
  use commondata
  use fields
  use kmc_data
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer :: i, j, k   ! Loop variables
  integer :: iii, jjj, kkk   ! Loop variables
  integer :: mc_step   ! KMC time loop

  real*8 :: tripper, cutoff
  integer :: fromx,fromy,tox,toy
  real*8 :: sum_kg

  character*10 :: curr_time

  integer, intent(in) :: iter

  integer :: ierr,status(MPI_STATUS_SIZE)

!!!!-------------------------------------------!!!!
!!!! PSEUDORANDOM NUMBER GENERATOR VARIABLES   !!!!
  character(LEN=10) :: date, time              !!!!
  integer :: dates(8)                          !!!!
  integer :: rand_seed_size                    !!!!
  integer, dimension(12) :: rand_seed          !!!!
  real :: random                               !!!!
!!!!-------------------------------------------!!!!


!!! Initialize PRNG
  call DATE_AND_TIME(date,time,VALUES=dates)
  rand_seed=(dates(1)+dates(2)+dates(3)+dates(5)+dates(6)+dates(7)+dates(8))*(rank+1)
!  write(6,*) "RANDSEED",rank,rand_seed
  rand_seed_size = 1
  call random_seed(size = rand_seed_size)
  call random_seed(put = rand_seed)
  !------------------------------

!!! Boundary conditions and swaps

  call swap_kmc()

  !! Initialize process list
  do x = 1,ksx+2
     do y = 1,ksy+2
        do i = 1,12
           plist(((y-1)*(ksx+2))+x,i)%fx = 0
           plist(((y-1)*(ksx+2))+x,i)%fy = 0
           plist(((y-1)*(ksx+2))+x,i)%from = 0
           plist(((y-1)*(ksx+2))+x,i)%tx = 0
           plist(((y-1)*(ksx+2))+x,i)%ty = 0
           plist(((y-1)*(ksx+2))+x,i)%to = 0
           plist(((y-1)*(ksx+2))+x,i)%prob = 0.0d0
        end do
     end do
  end do

  do x = 2,ksx+1
     do y = 2,ksy+1
        call add_to_plist(x,y)
     end do
  end do

  do mc_step = 1,60000

!     write(6,*) 'PER',mc_step,rank,sum(vfe_f_g)

     call mpi_barrier(MPI_COMM_WORLD,ierr)

     if (mod(mc_step,swap_freq_kmc).eq.0) then
        call swap_kmc()
     end if

     call random_number(random)

     cutoff = random*sum(plist%prob)

     tripper = 0.0d0

     do i = 1,(ksx+2)*(ksy+2)
        do j = 1,12
           tripper = tripper + plist(i,j)%prob
           if (tripper .ge. cutoff) then
              fromx = plist(i,j)%fx; fromy = plist(i,j)%fy
              tox = plist(i,j)%tx; toy = plist(i,j)%ty
              kg(plist(i,j)%fx,plist(i,j)%fy) = plist(i,j)%from
              kg(plist(i,j)%tx,plist(i,j)%ty) = plist(i,j)%to
              call upnl(fromx,fromy) ; call upnl(tox,toy)
              exit
           end if
        end do
        if (tripper .ge. cutoff) then
           exit
        end if
     end do

  end do




  call gather_kmc()

  if(rank.eq.0)then

     call system("rm -rf kgrid pgrid")
     open(unit=500,file='kgrid', status = 'new')
     do j = 1,ksy_g
        write(500,'(<ksx_g>(1I2))') (kg_g(i,j),i=1,ksx_g)
     end do
     close(500)

  end if

end subroutine kmcsolve





subroutine add_to_plist(x,y)
  use commondata
  use fields
  use kmc_data
  implicit none

  integer, intent(in) :: x,y
  integer :: i,j    ! Loop variables
  logical :: isfe
  integer :: nhbd   ! Neighborhood count


  nhbd = 0

  do i = x-1,x+1
     do j = y-1,y+1
        nhbd = nhbd + kg(i,j)
     end do
  end do

  ! Is it an Fe ion?
  isfe = XOR((mod(x,2)==0),(mod(y,2)==0))
  !------------------------------

  if (isfe) then   !! If it is a Fe cell

     if (kg(x,y).lt.1) then  !! If the Fe cell is a vacancy

       call make_ptl(x,y,x,y,1,1,vfe_a(x,y)*exp(9*nhbd/8.0d0))

        do i = x-1,x+1,2
           do j = y-1,y+1,2

              if (kg(i,j).gt.0) then
                 call make_ptl(x,y,i,j,1,0,v_diff(x,y)*exp(9*nhbd/8.0d0))
              end if

           end do
        end do


     else   !! If the Fe cell is not a vacancy

        call make_ptl(x,y,x,y,0,0,vfe_f(x,y)*exp(9*(9-nhbd)/8.0d0))

        do i = x-1,x+1,2
           j = y

           if (kg(i,j).gt.0) then
              call make_ptl(x,y,i,j,0,0,fes_diss(x,y)*exp(9*(9-nhbd)/7.0d0))
           end if

        end do


        do j = y-1,y+1,2
           i = x

           if (kg(i,j).gt.0) then
              call make_ptl(x,y,i,j,0,0,fes_diss(x,y)*exp(9*(9-nhbd)/7.0d0))
           end if

        end do


     end if

  else   !! If it is a S cell


     if (kg(x,y).lt.1) then  !! If the S cell is a vacancy

        call make_ptl(x,y,x,y,1,1,vs_a(x,y)*exp(9*nhbd/8.0d0))

        do i = x-1,x+1,2
           do j = y-1,y+1,2

              if (kg(i,j).gt.0) then
                 call make_ptl(x,y,i,j,1,0,v_diff(x,y)*exp(9*nhbd/8.0d0))
              end if

           end do
        end do


     else !! If the S cell is occupied

        call make_ptl(x,y,x,y,0,0,vs_f(x,y)*exp(9*(9-nhbd)/8.0d0))

     end if

  end if

end subroutine add_to_plist




subroutine make_ptl(fromx,fromy,tox,toy,from,to,prob)
  use commondata
  use fields
  use kmc_data

  implicit none

  integer, intent(in) :: fromx,fromy,tox,toy,from,to
  real*8, intent(in) :: prob
  integer :: lp
  integer :: atom_id
  logical :: ispending

  atom_id = ((fromy-1)*(ksx+2))+fromx

  ispending = .TRUE.

  do lp = 1, 12
     if (plist(atom_id, lp)%fx == 0) then
        plist(atom_id,lp)%fx = fromx ; plist(atom_id,lp)%fy = fromy
        plist(atom_id,lp)%tx = tox ; plist(atom_id,lp)%ty = toy
        plist(atom_id,lp)%from = from ; plist(atom_id,lp)%to = to
        plist(atom_id,lp)%prob = prob
        ispending = .FALSE.
        exit
     end if

  end do


  if (ispending) then
     write(6,*) 'You should no be seeing this'
  end if


end subroutine make_ptl




subroutine upnl(x,y)
  use commondata
  use fields
  use kmc_data
  implicit none

  integer, intent(in) :: x,y
  integer :: i,j

  do i = x-1,x+1
     do j = y-1,y+1

        if ((i.gt.1).and.(i.lt.ksx+2).and.(j.gt.1).and.(j.lt.ksy+2)) then
           call rfpls(i,j)
           call add_to_plist(i,j)
        end if

     end do
  end do


end subroutine upnl




subroutine rfpls(x,y)
  use commondata
  use fields
  use kmc_data
  implicit none

  integer, intent(in) :: x,y
  integer :: atom_id, proc_id

  atom_id = ((y-1)*(ksx+2))+x

  do proc_id = 1,12
     plist(atom_id,proc_id)%fx = 0 ; plist(atom_id,proc_id)%fy = 0
     plist(atom_id,proc_id)%tx = 0 ; plist(atom_id,proc_id)%ty = 0
     plist(atom_id,proc_id)%from = 0 ; plist(atom_id,proc_id)%to = 0
     plist(atom_id,proc_id)%prob = 0.0d0
  end do

end subroutine rfpls


