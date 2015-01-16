subroutine orsolve(iter)
  use commondata
  use fields
  implicit none
  include 'mpif.h'

  integer :: x, y, z   ! Loop variables
  integer, intent(in) :: iter
  integer :: ierr

  real*8, dimension(psx,psy,psz) :: M_opyr
  real*8 :: M_opyr_max = 1.5E-10
  real*8 :: M_opyr_min = 1.5E-13

  real*8, dimension(psx,psy,psz) :: D_opyr
  real*8, dimension(psx,psy,psz) :: del_opyr
  real*8 :: D_opyr_max = 1E-6   !! Truncation for the orientation field

  real*8 :: delx,dely,delz
  real*8 :: epsilon = 1E-18

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

  logical :: is_sorted
  integer :: rowindex, JAleft, JAright, JAswap
  real*8 :: Aswap

  if (mod(iter,swap_freq_pf).eq.1) then
     call swap_or()
  end if

  M_opyr = 0.0d0   !! Initialize to 0

  do x = 1,psx
     do y = 1,psy
        do z = 1,psz

           M_opyr(x,y,z) = min(M_opyr_min/((pyr(x,y,z+1)+0.01d0)**10),M_opyr_max)

           delx = odiff(opyr(wrap(x+1,psx),y,z+1),opyr(wrap(x-1,psx),y,z+1))
           dely = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,wrap(y-1,psy),z+1))
           delz = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1-1))

           del_opyr(x,y,z) = (sqrt((delx*delx)+(dely*dely)+(delz*delz))/dpf) + 1E-15

           D_opyr(x,y,z) = ((pyr(x,y,z+1)*pyr(x,y,z+1))/del_opyr(x,y,z)) + epsilon
           D_opyr(x,y,z) = min(D_opyr(x,y,z),D_opyr_max)

        end do
     end do
  end do


  approxsol = 0.0d0

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

           orc = odiff(opyr(x,wrap(y+1,psy),z+1),opyr(x,y,z+1))-(opyr(x,wrap(y+1,psy),z+1)-opyr(x,y,z+1))
           orc4 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1),opyr(x,y,z+1-1))-(opyr(x,y,z+1)-opyr(x,y,z+1-1))
           orc5 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf)

           orc = odiff(opyr(x,y,z+1+1),opyr(x,y,z+1))-(opyr(x,y,z+1+1)-opyr(x,y,z+1))
           orc6 = orc * M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,psz))+D_opyr(x,y,z))/(dpf*dpf)

           orc = (orc1+orc3+orc5)-(orc2+orc4+orc6)

           B(linindex) = (opyr(x,y,z+1)/dt) + orc

           if (z.eq.1) then
              B(linindex) = B(linindex) + ((M_opyr(x,y,z)*0.5d0*(D_opyr(x,y,z)+D_opyr(x,y,max(z-1,1)))/(dpf*dpf))*opyr(x,y,z+1-1))
           elseif (z.eq.psz) then
              B(linindex) = B(linindex) + ((M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,min(z+1,psz))+D_opyr(x,y,z))/(dpf*dpf))*opyr(x,y,z+1+1))
           end if

           approxsol(linindex) = opyr(x,y,z+1)

        end do
     end do
  end do


  !! Calculate the LHS (A matrix in Ax=B)
  allocate(A((7*psx*psy*psz)-(2*psx*psy)))
  allocate(JA((7*psx*psy*psz)-(2*psx*psy)))
  allocate(LU((7*psx*psy*psz)-(2*psx*psy)))
  allocate(IA((psx*psy*psz)+1))

  IA=0 ; JA=0 ; A=0
  contindex = 0
  linindex = 0

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx

           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x

           IA(linindex) = contindex + 1

           if (z .gt. 1) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z-1,psz))+D_opyr(x,y,z)))/(dpf*dpf)
              JA(contindex) = ((wrap(z-1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y-1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y-1,psy)-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x-1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x-1,psx)

           contindex = contindex + 1
           A(contindex) = D_opyr(x,y,wrap(z+1,psz))+D_opyr(x,y,wrap(z-1,psz))+&
                &D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,wrap(y-1,psy),z)+&
                &D_opyr(wrap(x+1,psx),y,z)+D_opyr(wrap(x-1,psx),y,z)
           A(contindex) = A(contindex) + 6*D_opyr(x,y,z)
           A(contindex) = A(contindex)*0.5d0*M_opyr(x,y,z)/(dpf*dpf)
           A(contindex) = A(contindex) + (1.0d0/dt)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + x

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(wrap(x+1,psx),y,z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((y-1)*psx) + wrap(x+1,psx)

           contindex = contindex + 1
           A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,wrap(y+1,psy),z)+D_opyr(x,y,z)))/(dpf*dpf)
           JA(contindex) = ((z-1)*psx*psy) + ((wrap(y+1,psy)-1)*psx) + x

           if (z .lt. psz) then
              contindex = contindex + 1
              A(contindex) = 0.0d0 - (M_opyr(x,y,z) * 0.5d0*(D_opyr(x,y,wrap(z+1,psz))+D_opyr(x,y,z)))/(dpf*dpf)
              JA(contindex) = ((wrap(z+1,psz)-1)*psx*psy) + ((y-1)*psx) + x
           end if


        end do
     end do
  end do

  IA((psx*psy*psz)+1)= IA(psx*psy*psz)+6



  do linindex = 1,psx*psy*psz
     is_sorted = .FALSE.

     do while (is_sorted .eq. .FALSE.)

        do rowindex = IA(linindex),IA(linindex+1)-2

           is_sorted = .TRUE.

           JAleft = JA(rowindex)
           JAright = JA(rowindex+1)

           if (JAleft .gt. JAright) then

              JAswap = JA(rowindex)
              JA(rowindex) = JA(rowindex+1)
              JA(rowindex+1) = JAswap

              Aswap = A(rowindex)
              A(rowindex) = A(rowindex+1)
              A(rowindex+1) = Aswap

              is_sorted = .FALSE.
              exit
           end if

        end do


     end do

  end do


  call iccglu(A,int(psx*psy*psz),IA,JA,LU,B,approxsol,scratch1,scratch2,scratch3,1E-3,100,iterations,0,solver_info)

  do z = 1,psz
     do y = 1,psy
        do x = 1,psx
           linindex = ((z-1)*psx*psy) + ((y-1)*psx) + x
           if (pyr(x,y,z+1).gt.0.1d0) then

              if (approxsol(linindex).eq.approxsol(linindex)) then
                 opyr(x,y,z+1) = approxsol(linindex)
              end if

           end if
        end do
     end do
  end do

end subroutine orsolve



double precision function odiff(or1,or2)
  implicit none
  real*8, intent(in) :: or1,or2
  real*8 :: Pi = 3.14159265d0

  if (or1>or2) then

     if ((or1-or2).lt.(or2+(Pi/2.0d0)-or1)) then
        odiff = -(or1-or2)
     else 
        odiff = (or2+(Pi/2.0d0)-or1)
     end if

  else

     if ((or2-or1).lt.(or1+(Pi/2.0d0)-or2)) then
        odiff = or2-or1
     else
        odiff = -(or1+(Pi/2.0d0)-or2)
     end if

  end if

  return
end function odiff




integer function wrap(a,lim)
  implicit none
  integer, intent(in) :: a,lim
  wrap = modulo(a-1,lim)+1
  return
end function wrap
