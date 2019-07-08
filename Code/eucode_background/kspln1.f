c+++++++++++++++++++++++++++++++++++++++++++++++
c creates spline for LiH+ -> Li + H+ equilibrium
c constants
c
c reads equil1.lih+
c puts spline parameters on common block
c equi1
c________________________________________
c written by P. C. Stancil, 6-25-96
c----------------------------------------
c
c np  number of data points
c
c e   temperature in K
c
c cs1  equilibrium constant
c
c cs12 spline parameters
c------------------------------------------
      subroutine kspln1 
c
      integer ie,np
      parameter(np=30)
      double precision csp1,cspn
      double precision e(np),cs1(np),cs12(np)
      character*70 junk
      common /equi1/ e,cs1,cs12
      external spline,splint
c
      open(unit=18, file ='equil1.lih+',status='old')
c
        read (18,99) junk
        read (18,99) junk
        read (18,99) junk
 99     format(a70)
c reads in constant data 
      do 10 ie=1,np
         read (18,*) e(ie), cs1(ie)
         cs1(ie)=dlog10(cs1(ie))
 10   continue
c
      csp1=0.0d0
      cspn=0.0d0
c
c obtain spline parameters
c
      call spline(e,cs1,np,csp1,cspn,cs12)
c
      return
      end

