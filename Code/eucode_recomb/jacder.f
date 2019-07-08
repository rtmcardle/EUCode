c     
c     order of species: y(i)=H-,H2+,H2.H3+, HD+,HD,H2D+
c     w(1)=H(n=1), w(2)=e, w(3)=H(n=2)
c     tr=radiation temperature
c     tc=matter temperature
c     a=rate coefficient
c     k=equilibrium constant, dissociation or ionization
c     den = H density
c     nhd=[D]/[H]
c
      SUBROUTINE jacobn(x,y,dfdx,dfdy,n,nmax)
      INTEGER n,nmax,i,j
      DOUBLE PRECISION x,y(*),dfdx(*),dfdy(nmax,nmax),
     . den,nhd,w(3),tc,ye(2)
c      COMMON /prop/ w,nhd
      common /temp/ tc
      external fcnj,estate
c
c      call estate(2,tc,ye)
c      y(24)=y(1)*ye(1)
c      y(22)=y(5)*ye(2)
      do 11 i=1,nmax
        dfdx(i)=0.0d0
11    continue
c
      call fcnj(n,y,dfdy)
c       do 10 j=1,n
c          write(6,*) 'dfdy= ',dfdy(1,1),dfdy(2,1),dfdy(4,1)
c 10    continue
      return
      END

      SUBROUTINE derivs(x,y,dydx)
      integer n
      DOUBLE PRECISION x,y(*),dydx(*),den,nhd,w(3)
     . ,tc,ye(2)
c      COMMON /prop/ w,den,nhd
      common /ndata/ n
      common /temp/ tc
      external fcn1,estate
c
c      call estate(2,tc,ye)
c      y(24)=y(1)*ye(1)
c      y(22)=y(5)*ye(2)
c
       call fcn1(n,y,dydx)
c       write(*,*) 'dydx=',dydx(1),dydx(2),dydx(4)
c       write(*,*) 'y= ',y(1),y(2),y(4)
c
c      write(*,*) 'check', dydx(10), a(24,tr)
c      print *, dydx(5),dydx(6),'tr=',tr,'tc=',tc
      return
      END
