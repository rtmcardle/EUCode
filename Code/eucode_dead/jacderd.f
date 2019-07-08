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
      INTEGER n,nmax,i
      DOUBLE PRECISION x,y(*),dfdx(*),dfdy(nmax,nmax),
     . a,k,tc,den,tr,nhd,w(3)
      COMMON /prop/ w,den,nhd
      common /temp/ tr,tc
      do 11 i=1,7
        dfdx(i)=0.0d0
11    continue
c  H-
      dfdy(1,1)=-a(1,tr)*k(1,tr)
     . -w(2)*a(2,tc)*den
     . -a(3,tc)*w(1)*den
     . -a(6,tc)*y(2)*den
     . -a(6,tc)*y(5)*den
     . -a(7,tc)*w(2)*den
c     . -a(13,tc)*y(7)*den
c     . -a(26,tc)*y(10)*den
c     . -a(31,tc)*y(8)*den
      dfdy(1,2)=0.0d0
     . -a(6,tc)*y(1)*den
c     . -a(7,tc)*y(1)*den
      dfdy(1,3)=0.0d0
      dfdy(1,4)=0.0d0
      dfdy(1,5)=0.0d0
     . -a(6,tc)*y(1)*den
      dfdy(1,6)=0.0d0
      dfdy(1,7)=0.0d0
c H2+
      dfdy(2,1)=0.0d0
     .    -a(6,tc)*y(2)*den
     .    +a(7,tc)*w(2)*den
c     .    -a(26,tc)*y(10)*den
      dfdy(2,2)=0.0d0
     . -a(4,tr)*k(2,tr)
     . -a(5,tc)*w(1)*den
     . -a(6,tc)*y(1)*den
     . -a(8,tc)*w(2)*den
     . -a(11,tc)*y(3)*den
     . -a(23,tc)*w(1)*nhd*deny
      dfdy(2,3)=0.0d0
     . -a(11,tc)*y(2)*den
      dfdy(2,4)=0.0d0
      dfdy(2,5)=0.0d0
      dfdy(2,6)=0.0d0
      dfdy(2,7)=0.0d0
c H2
      dfdy(3,1)=0.0d0
     .    +a(3,tc)*w(1)*den
     .    +a(6,tc)*y(2)*den
      dfdy(3,2)=0.0d0
     .    +a(5,tc)*w(1)*den
     .    +a(6,tc)*y(1)*den
     .    -a(11,tc)*y(3)*den
     .    +a(32,tc)*w(1)*nhd*den 
      dfdy(3,3)=0.0d0
     .    -a(10,tc)*w(3)*den
     .    -a(11,tc)*y(3)*den
     .    -a(18,tc)*w(2)*nhd*den 
     .    -a(20,tc)*w(1)*nhd*den
     .    -a(27,tc)*y(5)*den
     .    -a(28,tc)*y(5)*den 
     .    -a(30,tc)*w(3)*nhd*den
      dfdy(3,4)=0.0d0
     .    +a(12,tc)*w(2)*den/5.0d0
     .    -a(24,tc)*y(6)*den
      dfdy(3,5)=0.0d0
     .    -a(27,tc)*y(3)*den
     .    -a(28,tc)*y(3)*den
      dfdy(3,6)=0.0d0
     .    +a(19,tc)*w(1)*nhd*den
     .    -a(24,tc)*y(4)*den
      dfdy(3,7)=0.0d0
     .    +a(25,tc)*w(2)*den*0.07d0
c     .    +a(25,tc)*y(10)*den
c      dfdy(9,8)=0.0d0
c      dfdy(9,9)=0.0d0
c      dfdy(9,10)=0.0d0
c     .    +a(25,tc)*y(7)*den
c     .    +a(26,tc)*y(1)*den
c      dfdy(9,11)=0.0d0
c      dfdy(9,12)=0.0d0
c H_3^+
       dfdy(4,1)=0.0d0
       dfdy(4,2)=0.0d0
     . +a(11,tc)*y(3)*den
       dfdy(4,3)=0.0d0
     . +a(10,tc)*w(3)*den
     . +a(11,tc)*y(2)*den
     . +a(28,tc)*y(5)*den
       dfdy(4,4)=0.0d0
     . -a(12,tc)*w(2)*den
     . -a(24,tc)*y(6)*den
     . -a(26,tc)*w(1)*nhd*den
       dfdy(4,5)=0.0d0
     . +a(28,tc)*y(3)*den
       dfdy(4,6)=0.0d0
     . -a(24,tc)*y(4)*den  
       dfdy(4,7)=0.0d0
c HD+
      dfdy(5,1)=0.0d0
     . -a(6,tc)*y(5)*den 
      dfdy(5,2)=0.0d0
      dfdy(5,3)=0.0d0
     . -a(27,tc)*y(5)*den
     . -a(28,tc)*y(5)*den
      dfdy(5,4)=0.0d0
      dfdy(5,5)=0.0d0
     . -a(13,tr)*k(3,tr)
     . -a(14,tr)*k(4,tr)
     . -a(15,tc)*w(2)*den
     . -a(16,tc)*w(1)*den
     . -a(6,tc)*y(1)*den
     . -a(27,tc)*y(3)*den
     . -a(28,tc)*y(3)*den
      dfdy(5,6)=0.0d0
      dfdy(5,7)=0.0d0
c     .    +a(36,tc)*y(8)*nhd*den
c     .    +a(37,tc)*y(8)*nhd*den
c     .    -a(35,tc)*y(12)*den
c      dfdy(12,8)=0.0d0
c     .    +a(36,tc)*y(7)*nhd*den
c     .    +a(37,tc)*y(7)*nhd*den
c     .    -a(38,tc)*y(12)*den
c      dfdy(12,9)=0.0d0
c      dfdy(12,10)=0.0d0
c      dfdy(12,11)=0.0d0
c      dfdy(12,12)=0.0d0
c     .    -a(35,tc)*y(7)*den
c     .    -a(36,tr)*k(10,tr)
c     .    -a(37,tr)*k(11,tr)
c     .    -a(38,tc)*y(8)*den
c HD
      dfdy(6,1)=0.0d0
     .    +a(21,tc)*w(1)*nhd*den
c     .    +a(22,tc)*nhd*w(1)*den
      dfdy(6,2)=0.0d0
      dfdy(6,3)=0.0d0
     . +a(18,tc)*w(2)*nhd*den
     . +a(20,tc)*w(1)*nhd*den
      dfdy(6,4)=0.0d0
     . -a(24,tc)*y(6)*den
      dfdy(6,5)=0.0d0
     . +a(16,tc)*w(1)*den
      dfdy(6,6)=0.0d0
     . -a(17,tr)*k(5,tr)
     . -a(19,tc)*w(2)*den
     . -a(24,tc)*y(4)*den  
     . -a(29,tc)*w(3)*den
      dfdy(6,7)=0.0d0 
     .    +a(25,tc)*w(2)*den*0.20d0
c     .    +a(35,tc)*y(12)*den
c      dfdy(11,8)=0.0d0
c      dfdy(11,9)=0.0d0
c      dfdy(11,10)=0.0d0
c      dfdy(11,11)=0.0d0
c     .    -a(34,tr)*k(9,tr)
c      dfdy(11,12)=0.0d0
c     .    +a(35,tc)*y(7)*den
c
c H2D+
      dfdy(7,1)=0.0d0
      dfdy(7,2)=0.0d0
      dfdy(7,3)=0.0d0
     . +a(27,tc)*y(5)*den
     . +a(30,tc)*w(3)*nhd*den
      dfdy(7,4)=0.0d0
     . +a(24,tc)*y(6)*den
     . +a(26,tc)*w(1)*nhd*deny
      dfdy(7,5)=0.0d0
     . +a(27,tc)*y(3)*den
      dfdy(7,6)=0.0d0
     . +a(24,tc)*y(4)*den
     . +a(29,tc)*w(3)*den
      dfdy(7,7)=0.0d0
     . -a(25,tc)*w(2)*den
      return
      END

      SUBROUTINE derivs(x,y,dydx)
      DOUBLE PRECISION x,y(*),dydx(*),a,k,den,tc,tr,nhd,w(3)
      COMMON /prop/ w,den,nhd
      common /temp/ tr,tc
c H-
      dydx(1)=a(1,tc)*w(1)*w(2)*den
     . -y(1)*a(1,tr)*k(1,tr)
     . -a(2,tc)*y(1)*w(2)*den
     . -a(3,tc)*y(1)*w(1)*den
     . -a(6,tc)*y(1)*y(2)*den
     . -a(6,tc)*y(1)*y(5)*den
     . -a(7,tc)*y(1)*w(2)*den
c     . -a(13,tc)*y(1)*y(7)*den
c     . -a(26,tc)*y(10)*y(1)*den
c     . -a(31,tc)*y(1)*y(8)*den
c H2+
      dydx(2)=a(4,tc)*w(1)*w(2)*den
     .    -a(4,tr)*k(2,tr)*y(2)
     .    -a(5,tc)*w(1)*y(2)*den
     .    -a(6,tc)*y(1)*y(2)*den
     .    +a(7,tc)*w(2)*y(1)*den
     .    -a(8,tc)*y(2)*w(2)*den
c     .    +a(9,tc)*w(1)*w(3)*den
     .    -a(11,tc)*y(2)*y(3)*den
     .    -a(23,tc)*y(2)*w(1)*nhd*den
c H2
      dydx(3)=a(3,tc)*y(1)*w(1)*den
     .    +a(5,tc)*w(1)*y(2)*den
     .    +a(6,tc)*y(1)*y(2)*den
     .    -a(10,tc)*y(3)*w(3)*den
     .    -a(11,tc)*y(2)*y(3)*den
     .    +a(12,tc)*y(4)*w(2)*den/5.0d0
     .    -a(18,tc)*w(2)*nhd*y(3)*den
     .    +a(19,tc)*w(1)*nhd*y(6)*den
     .    -a(20,tc)*w(1)*nhd*y(3)*den
     .    +a(23,tc)*y(2)*w(1)*nhd*den
     .    -a(24,tc)*y(4)*y(6)*den
     .    +a(25,tc)*w(2)*y(7)*den*0.07d0
     .    -a(27,tc)*y(3)*y(5)*den
     .    -a(28,tc)*y(3)*y(5)*den
     .    -a(30,tc)*w(3)*nhd*y(3)*den
c H3+
      dydx(4)=a(10,tc)*y(3)*w(3)*den
     .    +a(11,tc)*y(2)*y(3)*den
     .    -a(12,tc)*y(4)*w(2)*den
     .    -a(24,tc)*y(4)*y(6)*den 
     .    -a(26,tc)*y(4)*w(1)*nhd*den
     .    +a(28,tc)*y(3)*y(5)*den
c HD+
      dydx(5)=a(13,tc)*w(1)*w(2)*nhd*den
     .    +a(14,tc)*w(2)*w(1)*nhd*den
     .    -a(15,tc)*w(2)*y(5)*den
     .    -a(13,tr)*k(3,tr)*y(5)
     .    -a(14,tr)*k(4,tr)*y(5)
     .    -a(16,tc)*w(1)*y(5)*den
     .    -a(6,tc)*y(1)*y(5)*den
     .    -a(27,tc)*y(3)*y(5)*den
     .    -a(28,tc)*y(3)*y(5)*den
c HD
      dydx(6)=a(16,tc)*y(5)*w(1)*den
     .    +a(17,tc)*w(1)*nhd*w(1)*den
     .    -a(17,tr)*k(5,tr)*y(6)
     .    +a(18,tc)*w(2)*nhd*y(3)*den
     .    -a(19,tc)*w(2)*y(6)*den
     .    +a(20,tc)*w(1)*nhd*y(3)*den
     .    +a(21,tc)*y(1)*w(1)*nhd*den
     .    +a(22,tc)*y(1)*nhd*w(1)*den
     .    -a(24,tc)*y(4)*y(6)*den
     .    +a(25,tc)*w(2)*y(7)*den*0.20d0
     .    -a(29,tc)*w(3)*y(6)*den
c
c H2D+
      dydx(7)=a(24,tc)*y(4)*y(6)*den
     .    -a(25,tc)*w(2)*y(7)*den
     .    +a(26,tc)*y(4)*w(1)*nhd*den
     .    +a(27,tc)*y(5)*y(3)*den
     .    +a(29,tc)*w(3)*y(6)*den
     .    +a(30,tc)*w(3)*nhd*y(3)*den
c      write(*,*) 'check', dydx(10), a(24,tr) 
c      print *, dydx(5),dydx(6),'tr=',tr,'tc=',tc
      return
      END

      double precision function a(i,t)
      double precision t
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,
     .  14,15,16,17,18,19,20,21,22,23,24,25,
     .  26,27,28,29,30) i
c
c H- radiative attachment
 1     a=3.0d-16*(t/300.0d0)**0.9*dexp(-t/1.3d4)
c       write(*,*) 'a1= ',a
       return
c
c H- H+ mutual neutralization 
 2     a=4.0d-8*(t/300.0d0)**(-0.5)
c       write(*,*) 'a2= ',a
       return
c
c H- + H -> H2 + e
c from Launay et al.
 3    a=1.5d-9*(t/300.0d0)**(-0.1)*dexp(-t/2.3d4)
      return
c
c H+ + H radiative association
 4     a=1.8d-18*(t/300.0d0)**(1.5)
       return
c
c H2+ + H -> H2 + H+
 5     a=6.4d-10
       return
c
c H2+ + H- -> H2 + H
 6     a=2.3d-7*(t/300.0d0)**(-0.5)
       return
c
c H+ + H- -> H2+ + e
 7     a=8.83d-15*(t/300.0d0)**(-0.32)
       return
c
c H2+ + e -> H + H
 8     a=1.68d-8*(t/300.0d0)**(-0.29)
       return
c
c H(n=2) + H -> H2+ + e, incorrect rate?
 9     a=1.68d-8*(t/300.0d0)**(-0.29)*dexp(-1.27d5/t) 
       return
c
c H2+ +H(n=2) -> H3+ + e
 10    a=1.0d-9
       return
c
c H2 + H2+ -> H3+ + H
 11    a=2.1d-9
       return
c
c H3+ + e -> H2+H
 12    a=2.0d-7*(T/300)**(-0.5)
       return
c
c D+ + H -> HD+ radiative association
 13    a=1.8d-18*(t/300.0d0)**(1.5)
       return
c
c H+ + D ->HD+ radiative association
 14    a=1.8d-18*(t/300.0d0)**(1.5)
       return
c
c HD+ + e -> D + H
 15    a=1.39d-8*(t/300.0d0)**(-0.29)
       return
c
c HD+ + H -> H+ + HD
 16    a=6.4d-10
       return
c
c H + D -> HD + nu
c Lepp & Shull (1984) guess
 17    a=1.0d-25
       return
c
c D+ + H2 -> H+ + HD
 18    a=1.7d-9
       return
c
c H+ + HD -> D+ + H2
 19    a=1.7d-10*dexp(-462.0d0/t)
       return
c
c D + H2 -> HD + H
 20    a=1.09d-18
       return
c
c D + H- -> HD + e
 21    a=1.5d-9*(t/300.0d0)**(-0.1)*dexp(-t/2.3d4)
       return
c
c D- + H -> HD + e
 22    a=1.5d-9*(t/300.0d0)**(-0.1)*dexp(-t/2.3d4)
       return
c
c H2+ + D -> D+ + H2
 23    a=6.4d-10
       return
c
c H3+ + HD -> H2D+ + H2
 24    a=1.7d-9
       return
c
c H2D+ + e -> H + H + D, H2 + D, H + HD
 25    a=6.0d-8
       return
c
c H3+ + D -> H2D+ + H
 26    a=1.8d-9
       return
c
c HD+ + H2 -> H2D+ + H
 27    a=1.1d-9
       return
c
c HD+ + H2 -> H3+ + D
 28    a=1.1d-9
       return
c
c H(n=2) + HD -> H2D+ + e
 29    a=1.0d-9
       return
c
c D(n=2) + H2 -> H2D+ + e
 30    a=1.0d-9
       return
c 
      end

      double precision function k(i,t)
      integer nk,nk1
      parameter(nk=58,nk1=30)
      double precision t,h0,h1,h2,h3,h4,h5,tk2(nk1),
     . ml,mm,ei,em,el,a0cm,pi,h,kb,lpf,hpf,hmpf,lipf,
     . a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,d0,lk,th,kk,
     . c0,c1,c2,c3,c4,doh2p,dh0,d1,d2,d3,d0hd,mu,
     . tk(nk),k1(nk),k2(nk1),k12(nk),k22(nk1),mlm,elm
c  partition function coefficients for h
      parameter (h0=-2.61655891d2,h1=1.63428326d2,h2=-4.06133526d1)
      parameter (h3=5.03282928d0,h4=-3.10998364d-1)
      parameter (h5=7.66654594d-3,ei=0.5d0,mu=0.99945568d0)
c  coefficients for h-
      parameter(em=2.76d-2,mm=0.99945616d0)
c  coefficients for li-
      parameter(elm=2.278d-2,mlm=0.99999218d0)
      parameter (kb=3.166829d-6,a0cm=0.529177249d-8)
c
c  Irwin partition function coefficients for li
      parameter (a0=-2.14567715d3,a1=1.41384065d3,a2=-3.72240759d2)
      parameter (a3=4.89662811d1,a4=-3.21834992d0)
      parameter (a5=8.4554492d-2,ml=0.9999215d0,el=0.1982d0)
c
c  dissociation coefficients for lih
      parameter (b0=9.9506d0,b1=-0.5163d0,b2=0.1833d0)
      parameter (b3=-1.7211d0, b4=1.2805d0, d0=2.4287d0)
c
c  dissociation coefficients for h2+
      parameter (c0=9.9835d0, c1=-0.0664d0, c2=-1.4979d0,
     . c3=-0.0195d0, c4=0.7486d0, doh2p=2.6508d0) 
c
c  dissociation coefficients for h2 (HD)
      parameter (dh0=11.1759d0, d1=-0.8735d0, d2=-0.7470,
     . d3=0.2748d0, d0hd=4.4781d0)
c
      common /equi/ tk,k1,k12
      common /equi1/ tk2,k2,k22
      external splint
c
      pi=4.0*atan(1.0d0)
      h=2.0d0*pi
c  hydrogen partition function
         lpf=h0+h1*dlog(t)+h2*(dlog(t))**2
     .    +h3*(dlog(t))**3+h4*(dlog(t))**4+h5*(dlog(t))**5
         hpf=dexp(lpf)
c lithium partition function
         lpf=a0+a1*dlog(t)+a2*(dlog(t))**2
     .    +a3*(dlog(t))**3+a4*(dlog(t))**4+a5*(dlog(t))**5
         lipf=dexp(lpf)
      go to (1,2,3,4,5) i
c
c  H- photodetachment coefficient
 1     k=(2.0d0*pi*mm*kb*t/h**2)**(1.5)*2.0d0*hpf/a0cm**3
     .  *dexp(-em/(kb*t))
c       write(*,*) 'k1= ',k
       return
c
c  H2+ to H+ + H photodissociation coefficient
 2     th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         k=10**(lk)*10.0d0/t/1.38066d-16
c       write(*,*) 'k2= ',k,t
       return
c
c  HD+ to D+ + H photodissociation coefficient
 3     th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         k=10**(lk)*10.0d0/t/1.38066d-16
c       write(*,*) 'k2= ',k
       return
c
c  HD+ to H+ + D photodissociation coefficient
 4     th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         k=10**(lk)*10.0d0/t/1.38066d-16
c       write(*,*) 'k2= ',k
       return
c  HD to D + H photodissociation coefficient
 5     th=5040.0d0/t
       lk=dh0+d1*dlog10(th)+d2*(dlog10(th))**2
     .      +d3*(dlog10(th))**3-th*d0hd
         k=10**(lk)*10.0d0/t/1.38066d-16
c       write(*,*) 'k2= ',k
       return
 999   format('k=',e16.8, 1x, i1, 1x,e16.8) 
       end
