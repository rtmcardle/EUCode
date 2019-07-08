c+++++++++++++++++++++++++++++++++++++++++++++++++++
c subroutine to calculate relative population of
c hydrogen (n=2) level
c
c output is array ye
c
c ________________________________________
c written by P. C. Stancil, 6-25-96
c----------------------------------------
c
c t    temperature in K
c
c ne   =2 for H and D 
c
c eu   binding energy of n=2 level
c
c gu   statistical weight of n=2 level
c
c hpf  hydrogen partition function
c
c k    boltzmann constant
c
c---------------------------------------
       subroutine estate(ne,t,ye)
       integer ne
       double precision ye(ne),ho,h1,h2,h3,h4,h5,eu
     .   ,lpf,hpf,gu,t,k

c
c  Irwin (1981, ApJ, 45, 624) Spartition function coefficients for h
      parameter (h0=-2.61655891d2,h1=1.63428326d2,h2=-4.06133526d1)
      parameter (h3=5.03282928d0,h4=-3.10998364d-1)
      parameter (h5=7.66654594d-3,eu=0.3748002d0,gu=4.0d0,
     . k=3.166829d-6)

c  hydrogen partition function
         lpf=h0+h1*dlog(t)+h2*(dlog(t))**2
     .    +h3*(dlog(t))**3+h4*(dlog(t))**4+h5*(dlog(t))**5
         hpf=dexp(lpf)
c
         ye(1)=gu/hpf*dexp(-eu/(k*t))
         ye(2)=ye(1)
c
       return
       end

