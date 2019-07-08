c function to determine equilibrium constants
c for photo-destruction processes
c
c processes index i and temperature t are input.
c ke is equilibrium constant returned
c
c LiH and LiH+ (Li + H+, Li+ + H) equilibrium
c constants are spline fits to data by routines
c kspln.f, kspln1.f, and kspln2.f with
c with spline parameters provided in
c common blocks equi, equi1, and equi2
c
c-------------------------------------------------
c written by P. C. Stancil, 6-25-96
c-------------------------------------------------
c
c  mu,mm,mlm,mhe,ml  reduced mass in amu
c
c  ei,em,elm,ehe,el  ionization energy in hartrees
c
c  kb  boltzmanns constant
c
c  a0cm bohr radius in cm
c
c  doh2p,doheh,dohd,ehe2 dissociation energy in eV
c
c------------------------------------------------- 
      double precision function ke(i,t)
      integer nk,nk1,nk2
      parameter(nk=58,nk1=30,nk2=109)
      double precision t,h0,h1,h2,h3,h4,h5,tk2(nk1),
     . ml,mm,ei,em,el,a0cm,pi,h,kb,lpf,hpf,hmpf,lipf,
     . a0,a1,a2,a3,a4,a5,d0,lk,th,kk,
     . c0,c1,c2,c3,c4,doh2p,dh0,d1,d2,d3,d0hd,mu,
     . tk(nk),k1(nk),k2(nk1),k12(nk),k22(nk1),mlm,elm
     . he0,he1,he2,he3,he4,he5,mhe,ep0,ep1,ep2,ep3,ep4
     . ep5,hepf,eppf,eh0,eh1,eh2,eh3,eh4,doheh,eh10,eh11,
     . eh12,eh13,eh14,he20,he21,he22,ehe2
     . ,k3(nk2),k32(nk2),tk3(nk2)
c
c  Irwin (1981, ApJS, 45, 624) partition function coefficients for 
c  atomic hydrogen
c
      parameter (h0=-2.61655891d2,h1=1.63428326d2,h2=-4.06133526d1)
      parameter (h3=5.03282928d0,h4=-3.10998364d-1)
      parameter (h5=7.66654594d-3,ei=0.5d0,mu=0.99945568d0)
c
c  parameters for h-
c  
      parameter(em=2.76d-2,mm=0.99945616d0)
c
c  parameters for li-
c
      parameter(elm=2.278d-2,mlm=0.99999218d0)
      parameter (kb=3.166829d-6,a0cm=0.529177249d-8)
c
c  Irwin partition function coefficients for helium
c
      parameter (he0=-3.76575219d-1,he1=2.33951687d-1,
     . he2=-5.79755525d-2,he3=7.1633316d-3,he4=-4.41302573d-4,
     . he5=1.08442997d-5,mhe=0.999862944d0,ehe=0.9035703d0)
c
c  Irwin partition function coefficients for he+
c
      parameter (ep0=6.93147179d-1,ep1=9.29636701d-10,
     . ep2=-2.30049742d-10,ep3=2.83829746d-11,ep4=-1.74590774d-12,
     . ep5=4.28355287d-14)
c
c  Irwin partition function coefficients for li
c
      parameter (a0=-2.14567715d3,a1=1.41384065d3,a2=-3.72240759d2)
      parameter (a3=4.89662811d1,a4=-3.21834992d0)
      parameter (a5=8.4554492d-2,ml=0.9999215d0,el=0.1982d0)
c
c  Dissociation coefficients for h2+ from Sauval & Tatum
c  (1984, ApJS, 56, 193)
c
      parameter (c0=9.9835d0, c1=-0.0664d0, c2=-1.4979d0,
     . c3=-0.0195d0, c4=0.7486d0, doh2p=2.6508d0) 
c
c  Dissociation coefficients for heh+ -> he + h+
c
      parameter (eh0=10.1667d0,eh1=-0.4814d0,eh2=-0.6444d0,
     . eh3=-0.7916d0,eh4=0.8534d0,doheh=1.845d0)
c
c  Dissociation coefficients for heh+ -> he+ + h from
c  Gaur & Tripathi 1985, J. Quant. Spectrosc. Radiat. Transfer,
c  33, 291. Not used since the fit is unreliable for T<4000 K
c
      parameter (eh10=12.1499222d0,eh11=-3.599778468d0,
     . eh12=2.87783529d0,eh13=-2.0632166d0,eh14=0.52229436d0)
c
c  Dissociation coefficients for h2 (HD)
c
      parameter (dh0=11.1759d0, d1=-0.8735d0, d2=-0.7470,
     . d3=0.2748d0, d0hd=4.4781d0)
c
c  dissociation coefficients for he2+ -> he + he+
c
      parameter (he20=10.2190d0,he21=-0.4011d0,he22=-0.6310d0
     .  ,ehe2=2.365d0)
c
      common /equi/ tk,k1,k12
      common /equi1/ tk2,k2,k22
      common /equi2/ tk3,k3,k32
      external splint
c
      pi=4.0*datan(1.0d0)
      h=2.0d0*pi
c
c  hydrogen partition function
         lpf=h0+h1*dlog(t)+h2*(dlog(t))**2
     .    +h3*(dlog(t))**3+h4*(dlog(t))**4+h5*(dlog(t))**5
         hpf=dexp(lpf)
c
c  helium partition function
         lpf=he0+he1*dlog(t)+he2*(dlog(t))**2
     .    +he3*(dlog(t))**3+he4*(dlog(t))**4+he5*(dlog(t))**5
         hepf=dexp(lpf)
c
c  he+ partition function
         lpf=ep0+ep1*dlog(t)+ep2*(dlog(t))**2
     .    +ep3*(dlog(t))**3+ep4*(dlog(t))**4+ep5*(dlog(t))**5
         eppf=dexp(lpf)
c
c lithium partition function
         lpf=a0+a1*dlog(t)+a2*(dlog(t))**2
     .    +a3*(dlog(t))**3+a4*(dlog(t))**4+a5*(dlog(t))**5
         lipf=dexp(lpf)
c
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) i
c
c  H photoionization coefficient
 1      ke=(2.0d0*pi*mu*kb*t/h**2)**(1.5)*2.0d0/hpf/a0cm**3
     .  *dexp(-ei/(kb*t))
       return
c
c  H- photodetachment coefficient
 2     ke=(2.0d0*pi*mm*kb*t/h**2)**(1.5)*2.0d0*hpf/a0cm**3
     .  *dexp(-em/(kb*t))
       return
c
c  He photoionization coefficient
 3     ke=(2.0d0*pi*mhe*kb*t/h**2)**(1.5)*2.0d0*eppf/hepf/a0cm**3
     .  *dexp(-ehe/(kb*t))
       return
c
c  Li photoionization coefficient
 4      ke=(2.0d0*pi*ml*kb*t/h**2)**(1.5)*2.0d0/lipf/a0cm**3
     .  *dexp(-el/(kb*t))
       return
c
c Li- photodetachment coefficient
 5     ke=(2.0d0*pi*mlm*kb*t/h**2)**(1.5)*2.0d0*lipf/a0cm**3
     .  *dexp(-elm/(kb*t))
       return
c
c  H2+ to H+ + H photodissociation coefficient
 6     th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         ke=10**(lk)*10.0d0/t/1.38066d-16
       return
c
c  HeH+ to H+ + He photodissociation coefficient
 7     th=5040.0d0/t
       lk=eh0+eh1*dlog10(th)+eh2*(dlog10(th))**2
     .      +eh3*(dlog10(th))**3+eh4*(dlog10(th))**4-th*doheh
         ke=10**(lk)*10.0d0/t/1.38066d-16
       return
c
c  HeH+ to H + He+ photodissociation coefficient
 8     th=5040.0d0/t
       lk=eh10+eh11*th+eh12*(th)**2
     .      +eh13*(th)**3+eh14*(th)**4
         ke=10**(lk)/t/1.38066d-16
       return
c
c  HD+ to D+ + H photodissociation coefficient
 9     th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         ke=10**(lk)*10.0d0/t/1.38066d-16
       return
c
c  HD+ to H+ + D photodissociation coefficient
 10    th=5040.0d0/t
       lk=c0+c1*dlog10(th)+c2*(dlog10(th))**2
     .      +c3*(dlog10(th))**3+c4*(dlog10(th))**4-th*doh2p
         ke=10**(lk)*10.0d0/t/1.38066d-16
       return
c
c  HD to D + H photodissociation coefficient
 11    th=5040.0d0/t
       lk=dh0+d1*dlog10(th)+d2*(dlog10(th))**2
     .      +d3*(dlog10(th))**3-th*d0hd
         ke=10**(lk)*10.0d0/t/1.38066d-16
       return
c
c  He2+ to He + He+ photodissociation coefficient
 12    th=5040.0d0/t
       lk=he20+he21*dlog10(th)+he22*(dlog10(th))**2
     .      -th*ehe2
         ke=10**(lk)/t/1.38066d-16
       return
c
c  LiH+ to Li+ + H photodissociation coefficient
 13      call splint(tk,k1,k12,nk,t,kk)
         ke=10**kk
       return
c
c  LiH+ to Li + H+ photodissociation coefficient
 14    if (t .ge. 600.0d0) then
          call splint(tk2,k2,k22,nk1,t,kk)
          ke=10**kk
       else
          ke=0.0d0
       end if
       return
c
c  LiH to Li + H photodissociation coefficient
 15    if (t .ge. 110.0d0) then
          call splint(tk3,k3,k32,nk2,t,kk)
          ke=10**kk
       else
          ke=0.0d0
       end if
       return
c
       end
