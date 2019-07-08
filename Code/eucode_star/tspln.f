c P.C. Stancil 10-17-94
c creates spline for equilibrium constants
c
c________________________________________
c
      subroutine tspln (n,zm,tm,tm2)
c
      integer n,ie
      double precision tmp1,tmpn
      double precision zm(n),tm(n),tm2(n)
      external spline,splint
c
      open(unit=19, file ='tmata.data',status='old')
c
c reads in constant data 
      do 10 ie=1,n
         read (19,*) zm(ie), tm(ie)
 10   continue
c
      tmp1=0.0d0
      tmpn=0.0d0
c
      call spline(zm,tm,n,tmp1,tmpn,tm2)
c
      return
      end
