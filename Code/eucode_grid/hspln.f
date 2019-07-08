c P.C. Stancil 3-8-95
c creates spline for Dell'Antonio-Rybicki hydrogen density 
c
c________________________________________
c
      subroutine hspln (n,zm,tm,tm2)
c
      integer n,ie,jj
      double precision tmp1,tmpn,dumb
      double precision zm(n),tm(n),tm2(n)
      external spline,splint
c
      open(unit=19, file ='hh+eIII.data')
c
c reads in constant data 
      do 10 ie=1,n
         read (19,*) zm(ie), dumb,tm(ie),dumb
 10   continue
c
      tmp1=0.0d0
      tmpn=0.0d0
c
      call spline(zm,tm,n,tmp1,tmpn,tm2)
c
      return
      end
