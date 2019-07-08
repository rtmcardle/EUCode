      subroutine fcn1p(z,n,x)
c      implicit double precision (a-h,o-z)
      integer r,ihnu,maxr,maxx,idim,ia,ib,nreacs,i,j,n
     . ,rj,pj,ii
      parameter (maxr =260)
      parameter (maxx =100)
      parameter (idim = 8)
      double precision x(maxx)
      double precision xinp(idim),rate(maxr),den,z
      dimension r(maxr,6)
      double precision xk(maxr),xkp(maxr,4),xr(2)
      common /rdata/ r,xk,xkp,nreacs
      common /idata/ den,xinp,ihnu
c
      rj=0
      pj=0
c
      do 108 i=1,nreacs
      ia=r(i,1)
      ib=r(i,2)
      if (ia.gt.0) then
          xr(1)=x(ia)
      elseif(ia.lt.0) then
          xr(1)=xinp(-ia)
c determine H* and D* fractional abundances
          if(ia.eq.-2) xr(1)=xr(1)*xinp(4)
          if(ia.eq.-3) xr(1)=xr(1)*x(2)
      else
          write(6,*) 'reaction # ',i
          stop 'can not happen ia=0  1'
      endif
      if (ib.gt.0) then
          xr(2)=x(ib)
      elseif(ib.lt.0) then
          xr(2)=xinp(-ib)
c determine H* and D* fractional abundances
          if(ia.eq.-2) xr(2)=xr(2)*xinp(4)
          if(ia.eq.-3) xr(2)=xr(2)*x(2)
      else
          stop 'can not happen ib=0 1'
      endif
c check
      if (i.eq.161) then
      write(70,*) den, xr(1), xr(2), xk(i)
      else
      continue
      endif
c
      rate(i)=den*xr(1)*xr(2)*xk(i)
      rj=rj+1
      if (rj.eq.30) then
 10      write(40+pj,999) z,(rate(ii+pj*30), ii=1,30)
         rj=0
         pj=pj+1
         if (i.eq. nreacs) goto 108
      end if 
      if (i.eq. nreacs) goto 10
 108  continue
 999  format(31(e10.4,1x))
      return
      end
