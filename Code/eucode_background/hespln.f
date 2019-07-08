c P.C. Stancil 10-28-98
c creates spline for Seager He+,He
c
c________________________________________
c
      subroutine hespln (n,z,hei,hei2,he,he2)
c
      integer n,i
      double precision tmp1,tmpn,dumb
      double precision z(n),hei(n),hei2(n),he(n),he2(n)
      external spline,splint
c
      open(unit=19, file ='model1_file3')
c
c reads in constant data 
      do 10 i=n,1,-1
         read (19,*) z(i),hei(i),he(i),dumb,dumb
 10   continue
c
      tmp1=0.0d0
      tmpn=0.0d0
c
      call spline(z,hei,n,tmp1,tmpn,hei2)
      call spline(z,he,n,tmp1,tmpn,he2)
c
      return
      end
