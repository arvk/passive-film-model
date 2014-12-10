subroutine initialize_kmc()
  use commondata
  use fields
  use kmc_data
  implicit none

!!!!-------------------------------------------!!!!
!!!! PSEUDORANDOM NUMBER GENERATOR VARIABLES   !!!!
  character(LEN=10) :: date, time              !!!!
  integer :: dates(8)                          !!!!
  integer :: rand_seed_size                    !!!!       
  integer, dimension(12) :: rand_seed          !!!!
  real :: random                               !!!!
!!!!-------------------------------------------!!!!

  integer, dimension(psx_g,psy_g) :: interface_loc
  integer :: x, y, z   ! Loop variables
  integer :: i, j, k   ! Loop variables
  integer :: ii, jj    ! Loop variables
  real*8 :: sum_kg
  real*8 :: my_ph  ! Avg pH used to calculate dissolutipon probability

!!! Initialize PRNG
  call DATE_AND_TIME(date,time,VALUES=dates)
  rand_seed=dates(1)+dates(2)+dates(3)+dates(5)+dates(6)+dates(7)+dates(8)
  rand_seed_size = 1
  call random_seed(size = rand_seed_size)  
  call random_seed(put = rand_seed)
  !------------------------------

!!! Identify location of the interface
  interface_loc = 0

  do x = 1,psx_g
     do y = 1,psy_g
        do z = psz_g-1,1,-1
           if ((pht_g(x,y,z) .gt. 5.0E-1).and.(pht_g(x,y,z+1) .lt. 5.0E-1)) then
              interface_loc(x,y) = z
              exit
           end if
        end do
     end do
  end do


  do x = 1,psx_g
     do y = 1,psy_g

        sum_kg = 0.0d0

        do i = 1,kg_scale
           do j = 1,kg_scale
              sum_kg = sum_kg + kg_g(((x-1)*kg_scale)+i,((y-1)*kg_scale)+j)
           end do
        end do
        sum_kg = sum_kg/(kg_scale*kg_scale)



        do while((sum_kg.lt.((pht_g(x,y,interface_loc(x,y))-0.10d0))).or.(sum_kg.gt.((pht_g(x,y,interface_loc(x,y))+0.10d0))))


           if (sum_kg.lt.(pht_g(x,y,interface_loc(x,y)))) then

              do i = 1,kg_scale
                 do j = 1,kg_scale

                    if (kg_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) .eq. 0) then
                       call random_number(random)
                       if (random .lt. ((pht_g(x,y,interface_loc(x,y))-sum_kg)/(1.0d0-sum_kg))) then
                          kg_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 1                      
                       end if
                    end if

                 end do
              end do

           else

              do i = 1,kg_scale
                 do j = 1,kg_scale

                    if (kg_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) .eq. 1) then
                       call random_number(random)
                       if (random .lt. ((sum_kg-pht_g(x,y,interface_loc(x,y)))/sum_kg)) then
                          kg_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0
                       end if
                    end if

                 end do
              end do


           end if






           sum_kg = 0.0d0

           do i = 1,kg_scale
              do j = 1,kg_scale
                 sum_kg = sum_kg + kg_g(((x-1)*kg_scale)+i,((y-1)*kg_scale)+j)
              end do
           end do
           sum_kg = sum_kg/(kg_scale*kg_scale)

        end do








     end do
  end do




  !! Initialize baseline probabilities
  do x = 1,psx_g
     do y = 1,psy_g

        my_ph = 0.0d0

        do z = interface_loc(x,y),psz_g
           my_ph = my_ph + (14.0d0 - pH_in)
        end do

        my_ph = my_ph/(psz_g-interface_loc(x,y))       

        do i = 1,kg_scale
           do j = 1,kg_scale

              vfe_f_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0.001d0
              vfe_a_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0.001d0
              vs_f_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0.001d0
              vs_a_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0.001d0
              fes_diss_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 0.01d0*env_g(x,y,interface_loc(x,y)+1)*exp(my_ph/10)
              v_diff_g(((x-1)*kg_scale + i), ((y-1)*kg_scale + j)) = 20.0d0*pht(x,y,interface_loc(x,y))

           end do
        end do


     end do
  end do






end subroutine initialize_kmc
