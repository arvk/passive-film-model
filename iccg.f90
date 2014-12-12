      subroutine matgen(ifli,iflm,iflo,irot,nonz,ia,ja,a,f)
      implicit real*8(a-h,o-z)
      integer ia(1),ja(1),ifli,iflo
      real*8 a(1)
      integer bwi,bwo
      real*8 f(28000),acof0(28000),acof1(28000),acof2(28000),acof3(28000),acof4(28000),acof5(28000),acof6(28000),gl(1000),gu(1000)
       real*8 cvx(30,30,30),cvy(30,30,30),cvz(30,30,30),omg(30,30,30)
       common/coef/acof0,acof1,acof2,acof3,acof4,acof5,acof6
       common /order/ iorflg
!!!!     ineum5=1 neuman on bot;=0 dir on bot.
!!!!     ineum6=1 neuman on top;=0 dir on top.
!!!!     irot = 1, rotational flow field
!!!!     irot = 0, non-rotational flow field
      flx=1.e0
      fly=1.e0
      flz=1.e0
!!!!     write(6,*)' call prbtyp'
       call prbtyp(ineum5,ineum6,ising,imx,jmx,kmx,1)
!!!!     write(6,*)' out of prbtyp in matgen'
      id=1
      idp1=id+1
       iptflg=0
       go to (3,5,7),iorflg
 3     continue
!!!!     write(6,*)' continue 3'
       bwo=imx*jmx
       bwi=imx
       ni=1
       nj=bwi
       nk=bwo
       nc=-bwo-bwi
       nib=1
       njb=imx
       ncb=-imx
!!!!      (ijk) (5310246)
       go to 9
 5     continue
!!!!     write(6,*)' continue 5'
!!!!      (kij) (3150624)
       bwo=kmx*imx
       bwi=kmx
       ni=bwi
       nj=bwo
       nk=1
       nc=-bwo-bwi
       nib=1
       njb=imx
       ncb=-imx
       go to 9
 7     continue
!!!!     write(6,*)' continue 7'
!!!!      (jki) (1530462)
       bwo=jmx*kmx
       bwi=jmx
       ni=bwo
       nj=1
       nk=bwi
       nc=-bwo-bwi
       nib=jmx
       njb=1
       ncb=-jmx
 9     continue
!!!!     write(6,*)' continue 9'
!!!!     write(6,*)'imx,jmx,kmx',imx,jmx,kmx
      dx=flx/float(imx)
      dy=fly/float(jmx)
      dz=flz/float(kmx)
      rdx2=1.e0/dx**2
      rdy2=1.e0/dy**2
      rdz2=1.e0/dz**2
      mx=imx*jmx*kmx
      imxm1=imx-1
      jmxm1=jmx-1
      kmxm1=kmx-1
      nonz=3*mx-2+2*(mx-bwi+ifli)+2*(mx-bwo+iflo)
      nonzt=nonz
!!!!     write(6,*)' call inits'
      call inits(a,mx,ia,ja)
      if(ifli .ge. bwi-1) go to 220
      if(iflo .ge. bwo-bwi) go to 220
      do 10 m=1,mx
      f(m)=0.e0
  10  continue
!!!!     write(6,*)' continue 10'
      do 11 i=1,imx
      do 11 j=1,jmx
      do 11 k=1,kmx
      xi=(float(i)-0.5e0)*dx
      yj=(float(j)-0.5e0)*dy
      zk=(float(k)-0.5e0)*dz
      facx = 1.0e0
      facy = 1.0e0
      if( irot .eq. 0 ) go to 12
      facx = xi - .05e0
      facy = yj - .05e0
   12 continue
      dum=32.0e0*xi*(1.e0-xi)*yj*(1.e0-yj)*zk
      dum1=2.0e0*xi*yj*zk**2
       dum=25.e0*dum
       dum1=2.0e0*dum1
       cvx(i,j,k)=dum*facx
       cvy(i,j,k)=dum*facy
       cvz(i,j,k)=dum1
       omg(i,j,k)=0.e0
       dum2 = (xi*xi)*yj*zk
       m = i*ni + j*nj + k*nk + nc
       f(m) = dum2
11    continue
!!!!     write(6,*)' continue 11'
      do 15 j=1,jmx
      do 15 i=1,imx
      mb=i*nib+j*njb+ncb
      gl(mb)=1.e0
      gu(mb)=2.0e0
  15  continue
!!!!     write(6,*)' continue 15'
      do 20 m=1,mx
      acof0(m)=0.e0
      acof1(m)=0.e0
      acof2(m)=0.e0
      acof3(m)=0.e0
      acof4(m)=0.e0
      acof5(m)=0.e0
      acof6(m)=0.e0
  20  continue
!!!!     write(6,*)' continue 20'
      if(imx .lt. 3) go to 45
      if(jmx .lt. 3) go to 45
      if(kmx .lt. 3) go to 45
      do 30 j=2,jmxm1
      do 30 i=2,imxm1
      do 30 k=2,kmxm1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=2.0e0*(rdx2+rdy2+rdz2)+omg(i,j,k)
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
  30  continue
!!!!     write(6,*)' continue 30'
  45  continue
!!!!     write(6,*)' continue 45'
      if((jmx.lt.3).or.(kmx.lt.3)) go to 55
      do 50 j=2,jmxm1
      do 50 k=2,kmxm1
      i=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+2.0e0*(rdy2+rdz2)+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      i=imx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+2.0e0*(rdy2+rdz2)+omg(i,j,k) +0.5e0*cvx(i,j,k)/dx
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
  50  continue
!!!!     write(6,*)' continue 50'
  55  continue
!!!!     write(6,*)' continue 55'
      if((imx.lt.3).or.(kmx.lt.3)) go to 65
      do 60 i=2,imxm1
      do 60 k=2,kmxm1
      j=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdy2+ 2.0e0*(rdx2+rdz2)+omg(i,j,k) -0.5e0*cvy(i,j,k)/dy
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      j=jmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdy2+2.0e0*(rdx2+rdz2)+omg(i,j,k)+0.5e0*cvy(i,j,k)/dy
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
  60  continue
!!!!     write(6,*)' continue 60'
  65  continue
!!!!     write(6,*)' continue 65'
      if((imx.lt.3).or.(jmx.lt.3)) go to 75
      do 70 i=2,imxm1
      do 70 j=2,jmxm1
      k=1
      m=(j-1)*bwo+(i-1)*bwi +k
      acof0(m)=2.0e0*(rdx2+rdy2)+3.0e0*rdz2 +omg(i,j,k)+cvz(i,j,k)/(2.0e0*dz)
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 71
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
   71 continue
!!!!     write(6,*)' continue 71'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=2.0e0*(rdx2+rdy2)+3.0e0*rdz2+omg(i,j,k)-cvz(i,j,k)/(2.0e0*dz)
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 72
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
   72 continue
!!!!     write(6,*)' continue 72'
  70  continue
!!!!     write(6,*)' continue 70'
  75  continue
!!!!     write(6,*)' continue 75'
      if(kmx.lt.3) go to 85
      do 80 k=2,kmxm1
      i=1
      j=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+2.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      j=jmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+2.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      i=imx
      m=(j-1)*bwo+(i-1)*bwi +k
      acof0(m)=rdx2+rdy2+2.0e0*rdz2+omg(i,j,k) +0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+ 0.5e0*cvz(i,j,k)/dz
      j=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+2.0e0*rdz2+omg(i,j,k)+0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
  80  continue
!!!!     write(6,*)' continue 80'
  85  continue
!!!!     write(6,*)' continue 85'
      if(imx.lt.3)go to 95
      do 90 i=2,imxm1
      k=1
      j=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=2.0e0*rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvy(i,j,k)/dy+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 91
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
   91 continue
!!!!     write(6,*)' continue 91'
      j=jmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=2.0e0*rdx2 +rdy2 +3.0e0*rdz2+omg(i,j,k)+0.5e0*cvy(i,j,k)/dy+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 92
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+ cvz(i,j,k)/dz)*gl(mb)
   92 continue
!!!!     write(6,*)' continue 92'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=2.0e0*rdx2+rdy2+3.0e0*rdz2+omg(i,j,k) +0.5e0*cvy(i,j,k)/dy-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 93
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
   93 continue
!!!!     write(6,*)' continue 93'
      j=1
      m=(j-1)*bwo+(i-1)*bwi +k
      acof0(m)=2.0e0*rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvy(i,j,k)/dy-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 94
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
   94 continue
!!!!     write(6,*)' continue 94'
  90  continue
!!!!     write(6,*)' continue 90'
  95  continue
!!!!     write(6,*)' continue 95'
      if(jmx.lt.3)go to 105
      do 100 j=2,jmxm1
      k=1
      i=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+ 2.0e0*rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 101
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  101 continue
!!!!     write(6,*)' continue 101'
      i=imx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+2.0e0*rdy2+3.0e0*rdz2+omg(i,j,k)+0.5e0*cvx(i,j,k)/dx+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 102
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  102 continue
!!!!     write(6,*)' continue 102'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+2.0e0*rdy2+3.0e0*rdz2+omg(i,j,k)+0.5e0*cvx(i,j,k)/dx-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 103
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  103 continue
!!!!     write(6,*)' continue 103'
      i=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+2.0e0*rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 104
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  104 continue
!!!!     write(6,*)' continue 104'
 100  continue
!!!!     write(6,*)' continue 100'
 105  continue
!!!!     write(6,*)' continue 105'
      i=1
      j=1
      k=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k) -0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy +0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dx
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 106
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  106 continue
!!!!     write(6,*)' continue 106'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy -0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 107
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  107 continue
!!!!     write(6,*)' continue 107'
      i=1
      j=jmx
      k=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)-0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy +0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 108
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  108 continue
!!!!     write(6,*)' continue 108'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k) -0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy -0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof2(m)=-rdx2+0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 109
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  109 continue
!!!!     write(6,*)' continue 109'
      i=imx
      j=1
      k=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k) +0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy +0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 111
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  111 continue
!!!!     write(6,*)' continue 111'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)+0.5e0*cvx(i,j,k)/dx-0.5e0*cvy(i,j,k)/dy -0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof4(m)=-rdy2+0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 112
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  112 continue
!!!!     write(6,*)' continue 112'
      i=imx
      j=jmx
      k=1
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k) +0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy +0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 -cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof6(m)=-rdz2+0.5e0*cvz(i,j,k)/dz
      if (ineum5 .eq. 1) go to 113
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2+cvz(i,j,k)/dz)*gl(mb)
  113 continue
!!!!     write(6,*)' continue 113'
      k=kmx
      m=i*ni+j*nj+k*nk+nc
      acof0(m)=rdx2+rdy2+3.0e0*rdz2+omg(i,j,k)  +0.5e0*cvx(i,j,k)/dx+0.5e0*cvy(i,j,k)/dy -0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) acof0(m) = acof0(m) - 2.0e0*rdz2 +cvz(i,j,k)/dz
      acof1(m)=-rdx2-0.5e0*cvx(i,j,k)/dx
      acof3(m)=-rdy2-0.5e0*cvy(i,j,k)/dy
      acof5(m)=-rdz2-0.5e0*cvz(i,j,k)/dz
      if (ineum6 .eq. 1) go to 114
      mb=i*nib+j*njb+ncb
      f(m)=f(m)+(2.0e0*rdz2-cvz(i,j,k)/dz)*gu(mb)
  114 continue
!!!!     write(6,*)' continue 114'
       if(ineum5*ineum6 .eq. 0) go to 117
       if( ising .eq. 1 ) go to 117
       i=1
       j=1
       k=1
       m=i*ni+j*nj+k*nk+nc
       acof2(m)=0.e0
       acof4(m)=0.e0
       acof6(m)=0.e0
       f(m)=0.e0
       i=2
       m=i*ni+j*nj+k*nk+nc
       acof1(m)=0.e0
       i=1
       j=2
       m=i*ni+j*nj+k*nk+nc
       acof3(m)=0.e0
       j=1
       k=2
       m=i*ni+j*nj+k*nk+nc
       acof5(m)=0.e0
 117    continue
!!!!     write(6,*)' continue 117'
!!!!  load coef
!!!!     write(6,*)'load matrix'
      do 150 i=1+bwo,mx
             j=i-bwo
             t=acof3(i)
             call put(t,a,mx,ia,ja,i,j)
 150   continue
!!!!     write(6,*)' continue 150'
       do 155 i=1+bwi,mx
              j=i-bwi
              t=acof1(i)
              call put(t,a,mx,ia,ja,i,j)
 155   continue
!!!!     write(6,*)' continue 155'
      do 160 i=2,mx
         j=i-1
         t=acof5(i)
         call put(t,a,mx,ia,ja,i,j)
 160   continue
!!!!     write(6,*)' continue 160'
       do 165 i=1,mx
        j=i
        t=acof0(i)
       call put(t,a,mx,ia,ja,i,j)
 165   continue
!!!!     write(6,*)' continue 165'
       do 170 i=1,mx-1
          j=i+1
          t=acof6(i)
          call put(t,a,mx,ia,ja,i,j)
 170   continue
!!!!     write(6,*)' continue 170'
       do 175 i=1,mx-bwi
          j=i+bwi
          t=acof2(i)
          call put(t,a,mx,ia,ja,i,j)
 175   continue
!!!!     write(6,*)' continue 175'
       do 180 i=1,mx-bwo
          j=i+bwo
          t=acof4(i)
          call put(t,a,mx,ia,ja,i,j)
 180   continue
!!!!     write(6,*)' continue 180'
       if(iflo .eq. 0) go to 195
       do 192 k = 1,iflo
       do 185 i=1+bwo-k,mx
          j=i-bwo+k
          t=0.0e0
          call put(t,a,mx,ia,ja,i,j)
 185   continue
!!!!     write(6,*)' continue 185'
      do 190 i=1,mx-bwo+k
         j=i+bwo-k
         t=0.0e0
         call put(t,a,mx,ia,ja,i,j)
 190   continue
 192   continue
!!!!     write(6,*)' continue 190'
 195   continue
!!!!     write(6,*)' continue 195'
       if(ifli .eq. 0) go to 210
       do 208 k = 1,ifli
       do 200 i=1+bwi-k,mx
       j=i-bwi+k
       t=0.0e0
       call put(t,a,mx,ia,ja,i,j)
 200   continue
!!!!     write(6,*)' continue 200'
       do 205 i=1,mx-bwi+k
          j=i+bwi-k
          t=0.0e0
          call put(t,a,mx,ia,ja,i,j)
 205   continue
 208   continue
!!!!     write(6,*)' continue 205'
 210   continue
!!!!     write(6,*)' continue 210'
       if( iflm .eq. 0 ) go to 219
       do 217 k = 1,iflm
          do 212 i = 1+bwi+k,mx
             j = i - bwi - k
             call put(0.0e0,a,mx,ia,ja,i,j)
  212     continue
          do 213 i = 1,mx-bwi-k
             j = i + bwi + k
             call put(0.0e0,a,mx,ia,ja,i,j)
  213     continue
  217 continue
  219 continue

!!!!     write(6,*) 'formula,actual',nonzt,ia(mx+1)
      go to 215
 220  write(6,*) 'ifli or iflo too large'
 215  continue
!!!!     write(6,*)' continue 215'
      return
      end
      subroutine prbtyp(ineum5,ineum6,ising,imx,imy,imz,iflg)
       integer ineum5,ineum6,ising,imx,imy,imz,iflg
       integer n5,n6,isg,ix,iy,iz
      if(iflg .eq. 0) go to 10
!!!!     write(6,*)' in prbtyp ix,iy,iz recall',ix,iy,iz
      ineum5=n5
      ineum6=n6
      ising=isg
      imx=ix
      imy=iy
      imz=iz
!!!!     write(6,*)' in prbtyp recall imx,imy,imz',imx,imy,imz
      go to 20
 10   continue
!!!!     write(6,*)' continue #'
!!!!     write(6,*)' in prbtyp init imx,imy,imz',imx,imy,imz
      n5=ineum5
      n6=ineum6
      isg=ising
      ix=imx
      iy=imy
      iz=imz
!!!!     write(6,*)   ' *********** problem type *************'
!!!!     if( ineum5 .eq. 0 ) then
!!!!        write(6,*)' *****       dir on bottom        *****'
!!!!     else
!!!!        write(6,*)' *****     neuman on bottom       *****' 
!!!!     endif
!!!!     if( ineum6 .eq. 0 ) then
!!!!        write(6,*)' *****         dir on top         *****'
!!!!     else
!!!!        write(6,*)' *****      neuman on top         *****'
!!!!     endif
!!!!     if( ising .eq. 0 ) then
!!!!        write(6,*)' *****         regular            *****'
!!!!     else
!!!!        write(6,*)' ***** singular if neuman on top and bottom *****'
!!!!     endif
!!!!     write(6,*)   ' ***** mesh cells in x dir',imx,'     *****'
!!!!     write(6,*)   ' ***** mesh cells in y dir',imy,'     *****'
!!!!     write(6,*)   ' ***** mesh cells in z dir',imz,'     *****'
 20   continue
!!!!     write(6,*)' continue #'
      return
      end
      subroutine iccglu(a,n,ia,ja,af,b,x,r,p,s,tol,maxits,its,job,info)

      integer n,ia(1),ja(1),maxits,its,job,info
      real*8 a(1),af(1),x(1),r(1),p(1),s(1),b(1),tol
!!!!    this routine performs preconditioned conjugate gradient on a
!!!!    sparse matrix. the preconditioner is an incomplete lu of
!!!!    the matrix.
!!!!    on entry:
!!!!       a  real*8()
!!!!          contains the elements of the matrix.
!!!!       n  integer
!!!!          is the order of the matrix.
!!!!       ia integer(n+1)
!!!!          contains pointers to the start of the rows in the arrays
!!!!          a and ja.
!!!!       ja integer()
!!!!          contains the column location of the corresponding elements
!!!!          in the array a.
!!!!       b real*8 (n)
!!!!          contains the right hand side.
!!!!       x real*8 (n)
!!!!          contains an estimate of the solution, the closer it
!!!!          is to the true solution the faster tthe method will
!!!!          converge.
!!!!       tol real*8
!!!!          is the accuracy desired in the solution.
!!!!       maxits integer
!!!!          is the maximun number of iterations to be taken
!!!!          by the routine.
!!!!       job integer
!!!!          is a flag to signal if incomplete factorization
!!!!          aready exits in array af.
!!!!          job = 0  perform incomplete factorization
!!!!          job = 1  skip incomplete factorization
!!!!    on output
!!!!       af real*8 ()
!!!!          contains the incomplete factorization of the matrix
!!!!          contained in the array a.
!!!!       x real*8 (n)
!!!!         contains the solution.
!!!!       r,p,s real*8 (n)
!!!!         these are scratch work arrays.
!!!!       its integer
!!!!         contains the number of iterations need to converge.
!!!!       info integer
!!!!         signals if normal termination.
!!!!         info = 0 method converged in its iterations
!!!!         info = 1 method converged, but exit occurred because
!!!!                  residual norm was less than sqrt(rtxmin).
!!!!         info = -999 method didnot converge in maxits iterations.
!!!!    the algorithm has the following form.
!!!!       form incomplete factors l and u
!!!!       x(0) <- initial estiate
!!!!       r(0) <- b - a*x(0)
!!!!       r(0) <- trans(a)*invtrans(l)*invtrans(u)*inv(u)*inv(l)*r(0)
!!!!       p(0) <- r(0)
!!!!       i    <- 0
!!!!       while r(i) > tol do
!!!!          s      <- inv(u)*inv(l)*a*p(i)
!!!!          a(i)   <- trans(r(i))*r(i)/(trans(s)*s)
!!!!          x(i+1) <- x(i) + a(i)*p(i)
!!!!          r(i+1) <- r(i) - a(i)*trans(a)*invtrans(l)*invtrans(u)*s
!!!!          b(i)   <- trans(r(i+1))*r(i+1)/(trans(r(i))*r(i))
!!!!          p(i+1) <- r(i+1) + b(i)*p(i)
!!!!          i      <- i + 1
!!!!       end
      real*8 ai,bi,rowold,rownew,xnrm,anrm
      real*8 sdot
      real*8 rtxmax,rtxmin
      common /gear14/ rtxmax,rtxmin
      data rtxmin/1.0e-16/
      info = 0
      anrm = abs(a(isamax(ia(n+1)-1,a,1)))
!!!!    form incomplete factors l and u
      if( job .ne. 0 ) go to 5
         call scopy(ia(n+1)-1,a,1,af,1)
         call lu(af,n,ia,ja)
    5 continue
!!!!    r(0) <- b - a*x(0)
      call cgres(a,n,ia,ja,x,b,r)
!!!!    r(0) <- trans(a)*invtrans(l)*invtrans(u)*inv(u)*inv(l)*r(0)
!!!!    inv(u)*inv(l)*r(0)
      call scopy(n,r,1,s,1)
      call ssol(af,n,ia,ja,s,1)
      call ssol(af,n,ia,ja,s,2)
!!!!    invtrans(l)*invtrans(u)*above
      call ssol(af,n,ia,ja,s,-2)
      call ssol(af,n,ia,ja,s,-1)
!!!!    trans(a)*above
      call mmult(a,n,ia,ja,r,s,-1)
!!!!    p(0) <- r(0)
      call scopy(n,r,1,p,1)
      rowold = sdot(n,r,1,r,1)
      i = 0
!!!!    while r(i) > tol do
      ai = 1.0d0
   10 continue
      xnrm = abs(x(isamax(n,x,1)))
!!!!    write(6,*) ' iter residual xnrm ai',i,sqrt(rowold)/(anrm*xnrm),
!!!!  $            xnrm,ai
      if( sqrt(rowold) .le. tol*(anrm*xnrm) ) go to 12
!!!!    if (rowold .le. rtxmin) go to 25
   13 continue
!!!!       s      <- inv(u)*inv(l)*a*p(i)
         call mmult(a,n,ia,ja,s,p,1)
         call ssol(af,n,ia,ja,s,1)
         call ssol(af,n,ia,ja,s,2)
!!!!       a(i)   <- trans(r(i))*r(i)/(trans(s)*s)
         ai = rowold/sdot(n,s,1,s,1)
!!!!       x(i+1) <- x(i) + a(i)*p(i)
         call saxpy(n,ai,p,1,x,1)
!!!!       r(i+1) <- r(i) - a(i)*trans(a)*invtrans(l)*invtrans(u)*s
         call ssol(af,n,ia,ja,s,-2)
         call ssol(af,n,ia,ja,s,-1)
         call mmult(a,n,ia,ja,b,s,-1)
         call saxpy(n,-ai,b,1,r,1)
!!!!       b(i)   <- trans(r(i+1))*r(i+1)/(trans(r(i))*r(i))
         rownew = sdot(n,r,1,r,1)
         bi = rownew/rowold
!!!!       p(i+1) <- r(i+1) + b(i)*p(i)
         call saxpy2(n,bi,p,1,r,1)
         rowold = rownew
!!!!       i      <- i + 1
         i = i + 1
         if( i .gt. maxits ) go to 20
         go to 10
   12 continue
   15 continue
      its = i
      return
   20 continue
      info = -999
      its = maxits
      return
!!!!    info = 1
!!!!    its = i
!!!!    return
      end
      subroutine axpy(a,n,ia,ja,i,js,je,t,y)
      integer n,ia(1),ja(1),i,js,je
      real*8 a(1),y(1)
      real*8 t
!!!!    this routine computes an axpy for a row of a sparse matrix
!!!!    with a vector.
!!!!    an axpy is a multiple of a row of the matrix is added to the
!!!!    vector y.
      is = ia(i)
      ie = ia(i+1) - 1
      do 10 ir = is,ie
         j = ja(ir)
         if( j .gt. je ) go to 20
         if( j .lt. js ) go to 10
         y(j) = y(j) + t*a(ir)
   10 continue
   20 continue
      return
      end
      subroutine saxpy2(n,da,dx,incx,dy,incy)
!!!!    constant times a vector plus a vector.
!!!!    uses unrolled loops for increments equal to one.
!!!!    jack dongarra, linpack, 3/11/78.
      real*8 dx(1),dy(1),da
      integer i,incx,incy,m,mp1,n  
      if(n.le.0)return
      if (da .eq. 0.0e0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!!!!       code for unequal increments or equal increments
!!!!         not equal to 1
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dx(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!!!!       code for both increments equal to 1
!!!!       clean-up loop
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dx(i) = dy(i) + da*dx(i)
        dx(i + 1) = dy(i + 1) + da*dx(i + 1)
        dx(i + 2) = dy(i + 2) + da*dx(i + 2)
        dx(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      real*8 function dot(a,n,ia,ja,i,js,je,b)
      integer n,ia(1),ja(1),i,js,je
      real*8 a(1),b(1)
      real*8 t
!!!!    this routine computes an inner product for a row of a sparse matri
!!!!    with a vector.
      t = 0.0e0
      is = ia(i)
      ie = ia(i+1) - 1
      do 10 ir = is,ie
         j = ja(ir)
         if( j .gt. je ) go to 20
         if( j .lt. js ) go to 10
         t = t + a(ir)*b(j)
   10 continue
   20 continue
      dot = t
      return
      end
      subroutine inits(a,n,ia,ja)
      integer n,ia(1),ja(1)
      real*8 a(1)
!!!!    this routine initializes storage for the sparse
!!!!    matrix arrays.
!!!!    note: the matrix is initialized to have zeroes
!!!!    on the diagonal.
      do 10 i = 1,n
         a(i) = 0.0e0
         ja(i) = i
         ia(i) = i
   10 continue
      ia(n+1) = n+1
      return
      end
      subroutine insert(t,a,n,ia,ja,i,j,k)
      integer n,ia(1),ja(1),i,j,k
      integer ip1,np1
      real*8 t,a(1)
!!!!    this routine rearranges the elements in arrays a and ja
!!!!    and updates array ia for the new element.
!!!!    write(6,1000)i,j,ia(i),ia(i+1),t
1000  format('  *** from insert i,j,ia(i),ia(i+1),t ',4i5,1x,e15.8)
      l1 = k
      l2 = ia(n+1) - 1
      do 10 lb = l1,l2
         l = l2 - lb + l1
         a(l+1) = a(l)
         ja(l+1) = ja(l)
   10 continue
      a(k) = t
      ja(k) = j
      ip1 = i + 1
      np1 = n + 1
      do 20 l = ip1,np1
         ia(l) = ia(l) + 1
   20 continue
      return
      end
      logical function locate(a,n,ia,ja,i,j,k)
      integer n,ia(1),ja(1),i,j
      real*8 a(1)
!!!!    this routine will locate the i,j-th element in the
!!!!    sparse matrix structure.
      is = ia(i)
      ie = ia(i+1) - 1
!!!!    row i is from is to ie in array a.
      do 10 l = is,ie
         if( j .gt. ja(l) ) go to 10
         if( j .ne. ja(l) ) go to 5
            locate = .true.
            k = l
            go to 20
    5    continue
         if( j .ge. ja(l) ) go to 10
            locate = .false.
            k = 0
         go to 20
   10 continue
!!!!    get here if should be at the end of a row.
      locate = .false.
      k = 0
   20 continue
      return
      end
      subroutine lu(a,n,ia,ja)
      logical locate
      integer n,ia(1),ja(1)
      integer ikp1,kp1
      real*8 a(1)
!!!!    this subroutine does incomplete gaussian elimenation
!!!!    on a sparse matrix. the matrix is stored in a sparse
!!!!    data structure.
!!!!    note: no pivoting is done
!!!!          and the factorization is incomplete,
!!!!          i.e., only where there exists a storage location
!!!!          will the operation take place.
      do 40 k = 1,n
         if( .not. locate(a,n,ia,ja,k,k,l2) ) go to 50
         if( a(l2) .eq. 0.0e0 ) go to 50
         kp1 = k + 1
         do 30 i = kp1,n
            if( .not. locate(a,n,ia,ja,i,k,ik) ) go to 25
               is = ik
               ie = ia(i+1) - 1
               kj = l2
               ke = ia(k+1) - 1
               a(ik) = a(ik)/a(kj)
               if( kj .eq. ke ) go to 30
               kj = kj + 1
               ikp1 = ik + 1
               do 20 j = ikp1,ie
   10             continue
                     if( kj .gt. ke ) go to 30
                     if( ja(kj) .ge. ja(j) ) go to 15
                        kj = kj + 1
                        go to 10
   15                continue
                  if( ja(kj) .gt. ja(j) ) go to 20
                  a(j) = a(j) - a(ik)*a(kj)
   20          continue
   25       continue
   30    continue
   40 continue
      return
   50 continue
!!!!    write(6,*)' value of k and a(k,k)',k,a(l2)
!!!!    write(6,*)' error zero on diagonal'
!!!!    write(6,*)' matrix probably specified incorrectly or'
!!!!    write(6,*)' incomplete factorization produces a singular matrix'
      return
      end
      subroutine mmult(a,n,ia,ja,b,x,job)
      integer n,ia(1),ja(1),job
      real*8 a(1),b(1),x(1)
!!!!    this routine performs matrix vector multiple
!!!!    for a sparse matrix structure
!!!!    job has the following input meanings:
!!!!        job = 1  matrix vector multiple
!!!!            = -1 matrix transpose vector multiple
!!!!            = 2  unit lower matrix vector multiple
!!!!            = -2 unit lower matrix transpose vector multiple
!!!!            = 3  upper matrix vector multiple
!!!!            = -3 upper matrix transpose vector multiple
      real*8 dot
!!!!    a*x
      if( job .ne. 1 ) go to 15
         do 10 i = 1,n
            b(i) = dot(a,n,ia,ja,i,1,n,x)
   10    continue
!!!!    trans(a)*x
   15 continue
      if( job .ne. -1 ) go to 35
         do 20 i = 1,n
            b(i) = 0.0e0
   20    continue
         do 30 i = 1,n
            call axpy(a,n,ia,ja,i,1,n,x(i),b)
   30    continue
!!!!    l*x   when l is unit lower
   35 continue
      if( job .ne. 2 ) go to 45
         do 40 i = 1,n
            b(i) = x(i) + dot(a,n,ia,ja,i,1,i-1,x)
   40    continue
!!!!    trans(l)*x   when l is unit lower
   45 continue
      if( job .ne. -2 ) go to 55
         do 49 i = 1,n
            b(i) = x(i)
   49    continue
         do 50 i = 1,n
            call axpy(a,n,ia,ja,i,i,n,x(i),b)
   50    continue
!!!!    u*x
   55 continue
      if( job .ne. 3 ) go to 65
         do 60 i = 1,n
            b(i) = dot(a,n,ia,ja,i,i,n,x)
   60    continue
!!!!    trans(u)*x
   65 continue
      if( job .ne. -3 ) go to 85
         do 70 i = 1,n
            b(i) = 0.0e0
   70    continue
         do 80 i = 1,n
            call axpy(a,n,ia,ja,i,1,i,x(i),b)
   80    continue
   85 continue
      return
      end
      subroutine put(t,a,n,ia,ja,i,j)
      integer n,ia(1),ja(1),i,j
      real*8 t,a(1)
!!!!    this routine will insert an element into the sparse matrix
!!!!    structure.
      is = ia(i)
      ie = ia(i+1) - 1
!!!!    write(6,100)i,j,ia(i),ia(i+1)
100   format('  *** from put i,j,ia(i),ia(i+1) ',4i7)
!!!!    row i is from is to ie in array a.
      do 10 k = is,ie
         if( j .gt. ja(k) ) go to 10
         if( j .ne. ja(k) ) go to 5
            a(k) = t
            go to 20
    5    continue
         if( j .ge. ja(k) ) go to 12
            call insert(t,a,n,ia,ja,i,j,k)
            go to 20
   12    continue
   10 continue
!!!!    get here if should be at the end of a row.
      k = ie + 1
      call insert(t,a,n,ia,ja,i,j,k)
      go to 30
   20 continue
   30 continue
      return
      end
      subroutine cgres(a,n,ia,ja,x,b,r)
      integer n,ia(1),ja(1)
      real*8 a(1),x(1),b(1),r(1)
!!!!    this routine computes a residual for a*x=b where
!!!!    a is in a sparse structure
      real*8 dot
      do 10 i = 1,n
         r(i) = b(i) - dot(a,n,ia,ja,i,1,n,x)
   10 continue
      return
      end
      subroutine ssol(a,n,ia,ja,b,job)
      integer n,ia(1),ja(1),job
      real*8 a(1),b(1)
!!!!    this routine solves a system of equations based on a sparse
!!!!    matrix date structure.
!!!!    the array b contains the right hand side on input
!!!!        and on output has the solution
!!!!    job has the value 1 if l*x = b is to be solved.
!!!!    job has the value -1 if trans(l)*x = b is to be solved.
!!!!    job has the value 2 if u*x = b is to be solved.
!!!!    job has the value -2 if trans(u)*x = b is to be solved.
      logical locate
      real*8 t
      real*8 dot
!!!!    job = 1  solve  l*x = b
      if( job .ne. 1 ) go to 15
!!!!       solve l*y=b
         do 10 i = 2,n
            b(i) = b(i) - dot(a,n,ia,ja,i,1,i-1,b)
   10    continue
   15 continue
      if( job .ne. 2 ) go to 25
!!!!       solve u*x=y
         do 20 ib = 1,n
            i = n - ib + 1
            t = dot(a,n,ia,ja,i,i+1,n,b)
            if( .not. locate(a,n,ia,ja,i,i,k) ) go to 30
            b(i) = (b(i) - t)/a(k)
   20    continue
!!!!    job = -2  solve  trans(u)*x = b
   25 continue
      if( job .ne. -2 ) go to 35
!!!!       solve trans(u)*y=b
         do 21 i = 1,n
            if( .not. locate(a,n,ia,ja,i,i,k) ) go to 30
            b(i) = b(i)/a(k)
            call axpy(a,n,ia,ja,i,i+1,n,-b(i),b)
   21    continue
!!!!       solve trans(l)*x=y
   35 continue
      if( job .ne. -1 ) go to 45
         do 22 ib = 2,n
            i = n - ib + 2
            call axpy(a,n,ia,ja,i,1,i-1,-b(i),b)
   22    continue
   45 continue
      return
   30 continue
!!!!    write(6,*)' error no diagonal element: from solve'
      return
      end
      subroutine saxpy(n,sa,sx,incx,sy,incy)
!!!!    constant times a vector plus a vector.
!!!!    uses unrolled loop for increments equal to one.
!!!!    jack dongarra, linpack, 3/11/78.
      real*8 sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!!!!       code for unequal increments or equal increments
!!!!         not equal to 1
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!!!!       code for both increments equal to 1
!!!!       clean-up loop
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end
      subroutine scopy(n,sx,incx,sy,incy)
!!!!    copies a vector, x, to a vector, y.
!!!!    uses unrolled loops for increments equal to 1.
!!!!    jack dongarra, linpack, 3/11/78.
      real*8 sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!!!!       code for unequal increments or equal increments
!!!!         not equal to 1
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!!!!       code for both increments equal to 1
!!!!       clean-up loop
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
      real*8 function sdot(n,sx,incx,sy,incy)
!!!!    forms the dot product of two vectors.
!!!!    uses unrolled loops for increments equal to one.
!!!!    jack dongarra, linpack, 3/11/78.
      real*8 sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!!!!       code for unequal increments or equal increments
!!!!         not equal to 1
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
!!!!       code for both increments equal to 1
!!!!       clean-up loop
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end
      integer function isamax(n,sx,incx)
!!!!    finds the index of element having max. absolute value.
!!!!    jack dongarra, linpack, 3/11/78.
      real*8 sx(1),smax
      integer i,incx,ix,n
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!!!!       code for increment not equal to 1
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
!!!!       code for increment equal to 1
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end

