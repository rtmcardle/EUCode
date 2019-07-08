c+++++++++++++++++++++++++++++++++++++++++++++++++
c subroutine to set-up jacobian of rate equations
c
c the number of equations (i.e., gas-species, 1 for each)
c and their fractional abundance (guess) are required.
c the jacobian is returned
c
c various data including rate coefficients are supplied
c in common blocks rdata and idata
c
c----------------------------------------
c written by S. Lepp ca 1984
c modified by P. C. Stancil, 6-25-96
c----------------------------------------
c pd jacobian matrix
c
c drate1,drate2 rate derivatives  
c
c x, xr, ia, and ib are defined in fcn1.f
c
c r, nreacs, xkp, xk, den ihnu, xinp, and idim,
c are defined in input.f
c
c----------------------------------------------
      subroutine fcnj(n,x,pd)
c      implicit double precision (a-h,o-z)
      integer r,ihnu,maxr,maxx,idim,i,j,n,nreacs,ia,ib
      parameter (maxr =260)
      parameter (maxx =100)
      parameter (idim = 8)
      double precision x(n),pd(maxx,maxx),den,
     . drate1,drate2
c  df(i)/dx(j) loaded into pd(i,j).
      double precision xinp(idim),xr(2)
      dimension r(maxr,6)
      double precision xk(maxr),xkp(maxr,4)
      common /rdata/ r,xk,xkp,nreacs
      common /idata/ den,xinp,ihnu
c
      do i = 1,maxx
      do j = 1,maxx
      pd(i,j)=0.0d0
      enddo
      enddo
c
c determine reactant type and correct
c fractional abundance
c

      do 108 i=1,nreacs
      ia=r(i,1)
      ib=r(i,2)
      if (ia.gt.0) then
          xr(1)=x(ia)
      elseif(ia.lt.0) then
          xr(1)=xinp(-ia)
c determine H* and D* fractional abundances
          if(ia.eq.-2) xr(1)=xinp(2)*xinp(4)
          if(ia.eq.-3) xr(1)=xinp(3)*x(2)
c          print*, xinp(2),xinp(3),xinp(4),x(2)
      else
          write(6,*) 'reaction # ',i
          stop 'can not happen ia=0  j'
      endif
      if (ib.gt.0) then
          xr(2)=x(ib)
      elseif(ib.lt.0) then
          xr(2)=xinp(-ib)
c determine H* and D* fractional abundances
          if(ia.eq.-2) xr(2)=xinp(2)*xinp(4)
          if(ia.eq.-3) xr(2)=xinp(3)*x(2)
      else
          stop 'can not happen ib=0 j'
      endif
c
c determine derivatives of reaction rates

      drate1=den*xr(2)*xk(i)
      drate2=den*xr(1)*xk(i)
c
c determine jacobian
c 
      if(ia.gt.0) then
        do 109 j= 1,3
 109       if(r(i,j).gt.0) pd(r(i,j),ia)=pd(r(i,j),ia)-drate1
        do 110 j=4,6
 110       if(r(i,j).gt.0) pd(r(i,j),ia)=pd(r(i,j),ia)+drate1
      endif
      if(ib.gt.0) then
        do 209 j= 1,3
 209       if(r(i,j).gt.0) pd(r(i,j),ib)=pd(r(i,j),ib)-drate2
        do 210 j=4,6
 210       if(r(i,j).gt.0) pd(r(i,j),ib)=pd(r(i,j),ib)+drate2
      endif
 108  continue
      return
      end

