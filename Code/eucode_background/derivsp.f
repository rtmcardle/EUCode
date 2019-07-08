c     
      SUBROUTINE derivsp(z,y)
      integer n
      DOUBLE PRECISION y(*),den,nhd,w(3)
     . ,tc,ye(2),z
c      COMMON /prop/ w,den,nhd
      common /ndata/ n
      common /temp/ tc
      external fcn1p,estate
c
c      call estate(2,tc,ye)
c      y(18)=y(1)*ye(1)
c      y(19)=y(5)*ye(2)
c
       call fcn1p(z,n,y)
c
      return
      END
