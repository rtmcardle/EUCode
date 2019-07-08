c+++++++++++++++++++++++++++++++++++++++++++++++++++
c subroutine to set-up rate equations (i.e., derivs)
c
c the number of equations (i.e., gas-species, 1 for each)
c and their fractional abundance (guess) are required.
c the derivatives are returned
c
c various data including rate coefficients are supplied
c in common blocks rdata and idata
c 
c----------------------------------------
c written by S. Lepp ca 1984
c modified by P. C. Stancil, 6-25-96
c----------------------------------------
c
c x   fractional abundance of each gas-particle species
c
c xr  temporary storage of x 
c
c ia,ib temporary storage of specied index number
c
c rate  total rate for each reaction
c
c xprime sum of all rates, rate equation
c
c r, nreacs, xkp, xk, den ihnu, xinp, and idim,
c are defined in input.f
c
c------------------------------------------------
      subroutine fcn1(n,x,xprime)
      integer r,ihnu,maxr,maxx,idim,ia,ib,nreacs,i,j,n
      parameter (maxr =260)
      parameter (maxx =100)
      parameter (idim = 8)
      double precision x(maxx),xprime(maxx)
      double precision xinp(idim),rate,den
      dimension r(maxr,6)
      double precision xk(maxr),xkp(maxr,4),xr(2)
      common /rdata/ r,xk,xkp,nreacs
      common /idata/ den,xinp,ihnu
c
      do 107,j=1,maxx
107   xprime(j)=0.0d0
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
      else
          write(6,*) 'reaction # ',i
          stop 'can not happen ia=0  1'
      endif
      if (ib.gt.0) then
          xr(2)=x(ib)
      elseif(ib.lt.0) then
          xr(2)=xinp(-ib)
c determine H* and D* fractional abundances
          if(ia.eq.-2) xr(2)=xinp(2)*xinp(4)
          if(ia.eq.-3) xr(2)=xinp(3)*x(2)
      else
          stop 'can not happen ib=0 1'
      endif
c
c determine reaction rate
c
      rate=den*xr(1)*xr(2)*xk(i)
c
c sum-up reaction rates for all reactions
c
c     destruction processes
c
      do 109 j= 1,3
 109  if(r(i,j).gt.0) xprime(r(i,j))=xprime(r(i,j))-rate
c
c     formation processes
c
      do 110 j=4,6
 110  if(r(i,j).gt.0) xprime(r(i,j))=xprime(r(i,j))+rate
 108  continue
c
      return
      end

