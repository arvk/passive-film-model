subroutine orsolve(iter)
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer, intent(in) :: iter

  real*8, dimension(psx,psy,psz) :: M_opyr
  real*8 :: M_opyr_max = 1E-9
  real*8 :: M_opyr_min = 1E-14

  real*8, dimension(psx,psy,psz) :: D_opyr
  real*8, dimension(psx,psy,psz) :: del_opyr
  real*8 :: D_opyr_max = 1E-8   !! Truncation for the orientation field

  real*8 :: delx,dely,delz
  real*8 :: epsilon = 1E-12

  real*8 :: odiff
  integer :: wrap

  integer :: linindex

  real*8 :: orc,orc1,orc2,orc3,orc4,orc5,orc6

  real*8, dimension(psx*psy*psz) :: B
  real*8, dimension(psx*psy*psz) :: approxsol
  real*8, dimension(psx*psy*psz) :: scratch1,scratch2,scratch3

  real*8, dimension(:), allocatable :: A, LU
  integer, dimension(:), allocatable :: JA, IA

  integer :: contindex
  integer :: iterations, solver_info
  logical :: solve_for_theta


  solve_for_theta = .FALSE.

  if (mod(iter,swap_freq_pf).eq.1) then
     call swap_or()
  end if

  M_opyr = 0.0d0   !! Initialize to 0

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz

           M_opyr(x,y,z) = M_opyr_min + (M_opyr_max-M_opyr_min)*pyr(x,y,z+1) !! Modify as necessary

           delx = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(wrap(x-1,psx),y,z+1))
           dely = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,wrap(y-1,psy),z+1))
           delz = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1-1))

           del_opyr(x,y,z) = sqrt((delx*delx)+(dely*dely)+(delz*delz))/dpf

           D_opyr(x,y,z) = ((pyr(x,y,z+1)*pyr(x,y,z+1))/del_opyr(x,y,z)) + epsilon
           D_opyr(x,y,z) = min(D_opyr(x,y,z),D_opyr_max)

        end do
     end do
  end do


  approxsol = 1.0d0

  !! Calculate the RHS (B matrix in Ax=B)
  !! B = C_{i,j} + o_pyr/dt

  orc = 0.0d0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           orc = odiff(opyr(x,y,z+1),opyr(wrap(x-1,psx),y,z+1))-(opyr(x,y,z+1)-opyr(wrap(x-1,psx),y,z+1))
           orc1 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(wrap(x-1,psx),y,z))/(dpf*dpf)

           orc = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(x,y,z+1))-(opyr(wrap(x+1,psx),y,z+1)-opyr(x,y,z+1))
           orc2 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1),opyr(x,wrap(y-1,psy),z+1))-(opyr(x,y,z+1)-opyr(x,wrap(y-1,psy),z+1))
           orc3 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,wrap(y-1,psy),z))/(dpf*dpf)

           orc = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,y,z+1))-(opyr(x,wrap(y+1,psx),z+1)-opyr(x,y,z+1))
           orc4 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1),opyr(x,y,z+1-1))-(opyr(x,y,z+1)-opyr(x,y,z+1-1))
           orc5 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1))-(opyr(x,y,z+1+1)-opyr(x,y,z+1))
           orc6 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,psz))+D_opyr(x,y,z))/(dpf*dpf)

           orc = (orc2+orc4+orc6)-(orc1+orc3+orc5)


           approxsol(linindex) = opyr(x,y,z+1)

           if (pyr(x,y,z+1).gt.0.1d0) then
              solve_for_theta = .TRUE.
              B(linindex) = (opyr(x,y,z+1)/dt) + orc
           else
              B(linindex) = 1E-1
           end if

        end do
     end do
  end do


  !! Calculate the LHS (A matrix in Ax=B)
  allocate(A(7*psx*psy*psz))
  allocate(JA(7*psx*psy*psz))
  allocate(LU(7*psx*psy*psz))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           IA(linindex) = contindex + 1

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z-1,psz))+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y-1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + (wrap(y-1,psy)-1)*psx + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x-1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + (y-1)*psx + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_opyr(x,y,wrap(z+1,psz))+D_opyr(x,y,wrap(z+1,psz))+&
                &D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,wrap(y-1,psy),z)+&
                &D_opyr(wrap(x+1,psx),y,z)+D_opyr(wrap(x-1,psx),y,z)
           A(contindex) = A(contindex) + 6*D_opyr(x,y,z)
           A(contindex) = A(contindex)*0.5d0*M_opyr(x,y,z)/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + (y-1)*psx + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + (y-1)*psx + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + (wrap(y+1,psy)-1)*psx + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z+1,psz))+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x

        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+7

  if (solve_for_theta) then

     call iccglu(A,psx*psy*psz,IA,JA,LU,B,approxsol,scratch1,scratch2,scratch3,1E-2,20,iterations,0,solver_info)

     do z = 1,psz
        do y = 1,psy
           do x = 1,psx
              linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
              if (pyr(x,y,z+1).gt.0.1d0) then

                 if (approxsol(linindex).eq.approxsol(linindex)) then
                    opyr(x,y,z+1) = min(max(approxsol(linindex),3.14d0/4.0d0),0.0d0)
                 end if

              end if
           end do
        end do
     end do

  end if


end subroutine orsolve


real function odiff(or1,or2)
  implicit none
  real*8, intent(in) :: or1,or2
  real*8 :: Pi = 3.14159265d0

  odiff = mod(abs(or1-or2),(2*Pi)/4)
  odiff = min(odiff,abs(odiff+(2*Pi/4)))

  return
end function odiff




integer function wrap(a,lim)
  implicit none
  integer, intent(in) :: a,lim
  wrap = modulo(a-1,lim)+1
  return
end function wrap

