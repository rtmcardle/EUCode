c  1-9-2000 PCS cooling/heating switches added
c  8-19-99 PCS Abel 1996 H2 cooling
c  2-9-99 PCS Compton cooling/heating added
c  create jacobian and rate equations
c  add cooling/heating functions by hand to
c    solve for the matter temperature
c
c    molecular line heating included by subtracting cooling
c      function at tr from tc cooling cooling function
c
c     tr=radiation temperature
c     tc=matter temperature
c     den = total H density
c     nhd=[D]/[H]
c
      SUBROUTINE jacobn(x,y,dfdx,dfdy,n,nmax)
      INTEGER n,nmax,i,j,ihnu
      DOUBLE PRECISION x,y(*),dfdx(*),dfdy(nmax,nmax),
     . den,nhd,w(3),tc,ye(2),gamma,mu,xinp(1),tr,dfdyh2,
     . dfdyhd,dfdylih,dfdyhehi,dfdyhdi,dfdyh3i
     .  ,h2lte,lh2low,h2low,dlow,dlte,ddem,h2lteh,
     .  lh2lowh,h2lowh,dlowh,dlteh,ddemh
     .  ,hdlte,lhdlow,hdlow,hdlteh,lhdlowh,hdlowh,
     .  lhdlte,lhdlteh,dfdycomp,kb,dfdy241,
     .  rkh,cvh,cvl,crh,crl,crlt,tch2,thh2, xlt,
     .  drkh,tdcvh,tdcvl,tdcrh,tdcrl,tdcrlt,tdtch2,
     .  dxlt,dcrl,dcvl,dtch2,dthh2,
     .  dfdyhion,dfdyheion,dfdyhepion,dfdyhrec,
     .  dfdyherec,dfdyhedr,dfdyhce,dfdyhece,dfdybrem,
     .  dfdytc(nmax),sumyi,para,everg,w0,den0,onoff(9),
     .   lihlow,lihlte,lihlowh,lihlteh,llihlowh,dfdyh2diss
      parameter(kb=1.380658d-16,sumyi=1.08d0,para=1.0d0,
     .  everg=1.60217733d-12)
c      COMMON /prop/ w,nhd
      common /idata/ den,xinp,ihnu
      common /temp/ tc,gamma,mu
c      common /temp/ tc
      common /prop/ tr,w0,den0
c      common /prop/ tr
      common /onoffd/ onoff
      external fcnj,fcnjt,estate
c
c      call estate(2,tc,ye)
c      y(24)=y(1)*ye(1)
c      y(22)=y(5)*ye(2)


      do 11 i=1,nmax
        dfdx(i)=0.0d0
11    continue
c
      call fcnj(n,y,dfdy)
      do 21 j=1,24
        dfdy(24,j)=0.0d0
21    continue
      print *, 'FCNJT!'
      call fcnjt(n,y,dfdytc)
      print *, 'DONE!'
      do 22 j=1,23
        dfdy(j,24)=dfdytc(j)
c don't need according to Shapiro & Kang ????
c        dfdy(j,j)=dfdy(j,j)-dend/den
22    continue
      tc=y(24)
c
c Total H2 cooling term
c   using Abel 1995 thesis H2 cooling function (ergs/s)
      dfdy(24,1)=-(gamma-1.0d0)*7.243d15/sumyi*para*y(14)*onoff(1)
     .   *4.37d-56*den**(1.5)*den**(-0.36*(dlog10(den))**(0.5))
     .   *(tc**16.05*(tc**dlog10(tc))**(-1.9)
     .   -tr**16.05*(tr**dlog10(tr))**(-1.9))
c H collisional ionization
     .   -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *1.27d-21*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-157809.1d0/tc)
c H collisional excitation
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *7.5d-19/(1.d0+(tc/1.d5)**0.5)*dexp(-118348.d0/tc)
c H2 formation cooling via H- and H2+
     . -(gamma-1.0d0)/kb*den/sumyi*everg*onoff(5)
     . *(3.00d0*1.5d-9*(tc/300.0d0)**-0.1*y(4)
     . + 1.83d0*6.4d-10*y(13))
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
     . -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)*onoff(6)
     . *4.78d-14*(tc/300.0d0)**1.8*exp(-5.31d4/tc)
c
      dfdy(24,2)=
c H Recombination
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *8.7d-27*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
c Bremsstrahlung of H+
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(8)
     . *1.42d-27*tc**0.5*(1.1d0+0.34d0
     . *dexp(-(5.5d0-dlog10(tc))**2/3.0d0))
c
      dfdy(24,4)=
c H2 formation cooling via H- and H2+
     .   -(gamma-1.0d0)/kb*y(1)*den/sumyi*everg*onoff(5)
     . *(3.53d0*1.5d-9*(tc/300.0d0)**-0.1)
c
      dfdy(24,8)=
c He collisional ionization
     .   -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *9.38d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-285335.4d0/tc)
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
     .    -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)*onoff(6)
     .  *1.47d-13*(tc/300.0d0)**1.39*exp(-1.31d5/tc)
c
      dfdy(24,9)=
c He+ collisional ionization
     .   -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *4.95d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-631515.d0/tc)
c He+ radiative recombination
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *1.55d-26*tc**0.3647
c Dielectronic recombination of He
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     .  *1.24d-13*tc**-1.5*dexp(-4.7d5/tc)
     .  *(1.d0+0.3d0*dexp(-9.4d4/tc))
c He+ collisional excitation
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(9)
     . *5.54d-17*tc**-0.397/(1.d0+(tc/1.d5)**0.5)*dexp(-473638.d0/tc)
c Bremsstrahlung of He+
     .  -(gamma-1.0d0)/kb*y(3)*den/sumyi*onoff(8)
     . *1.42d-27*tc**0.5*(1.1d0+0.34d0
     . *dexp(-(5.5d0-dlog10(tc))**2/3.0d0))
c
      dfdy(24,13) =
c H2 formation cooling via H- and H2+
     .   -(gamma-1.0d0)/kb*y(1)*den/sumyi*everg*onoff(5)
     . * 1.83d0*6.4d-10
c
c H2 cooling
c   using Abel 1995 thesis H2 cooling function (ergs/s)
      dfdy(24,14)=-(gamma-1.0d0)*7.243d15/sumyi*para*y(1)*onoff(1)
     .   *4.37d-56*den**(1.5)*den**(-0.36*(dlog10(den))**(0.5))
     .   *(tc**16.05*(tc**dlog10(tc))**(-1.9)
     .   -tr**16.05*(tr**dlog10(tr))**(-1.9))
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
     . -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)*onoff(6)
     .  *2.27d-19*(tc/300.0d0)**3.91*exp(-5.47e4/tc)*2.0d0
c
c HD Cooling
c Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
      lhdlte=-35.6998d0+15.35716d0*dlog10(tc)-5.58513d0*(dlog10(tc))**2
     . +0.8561149d0*(dlog10(tc))**3-1.75538d-2*(dlog10(tc))**4
      hdlte=10**lhdlte/den
c fit to Galli and Palla 1998 low density limit
c      lhdlow=-48.4421d0+30.13495d0*dlog10(tc)-14.4382d0*(dlog10(tc))**2
c     . +3.1536d0*(dlog10(tc))**3-0.249489d0*(dlog10(tc))**4
c      hdlow=10**(lhdlow)
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
      hdlow=(3.0d0*(4.4d-12+3.6d-13*tc**0.77)*dexp(-128.d0/tc)*128.d0
     .  + (5.0d0/3.0d0)*(4.1d-12+2.1e-13*tc**0.92)*dexp(-255.d0/tc)
     .  *255.d0)*kb/1.27d0
c HD heating
c Fit to Lepp and Shull 1984 LTE
      lhdlteh=-35.6998d0+15.35716d0*dlog10(tr)-5.58513d0*(dlog10(tr))**2
     . +0.8561149d0*(dlog10(tr))**3-1.75538d-2*(dlog10(tr))**4
      hdlteh=10**lhdlteh/den
c Galli and Palla 1998 low density limit
c      lhdlowh=-48.4421d0+30.13495d0*dlog10(tr)-14.4382d0*(dlog10(tr))**2
c     . +3.1536d0*(dlog10(tr))**3-0.249489d0*(dlog10(tr))**4
c      hdlowh=10**(lhdlowh)
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
      hdlowh=(3.0d0*(4.4d-12+3.6d-13*tr**0.77)*dexp(-128.d0/tr)*128.d0
     .  + (5.0d0/3.0d0)*(4.1d-12+2.1e-13*tr**0.92)*dexp(-255.d0/tr)
     .  *255.d0)*kb/1.27d0
c Total HD cooling term
      dfdy(24,1)=dfdy(24,1) -(gamma-1.0d0)*7.243d15*y(18)/sumyi*onoff(2)
     .   *(hdlte/(1.0d0+hdlte/hdlow)-hdlteh/(1.0d0+hdlteh/hdlowh))*den
c
c LiH cooling function
c Galli and Palla 1998 low density limit (erg cm3 /s)
      llihlow=-31.47d0+8.817d0*dlog10(tc)-4.144d0*(dlog10(tc))**2
     .   +0.8292d0*(dlog10(tc))**3-0.04996d0*(dlog10(tc))**4
      lihlow=10**(llihlow)
c Fit to LTE of Dalgarno 1984 (erg/s) -> (ergs cm3/s
      lihlte= 2.1d-11*(tc/400.0d0)**0.15*dexp(-1.9d3/tc)/den
       if(tc.le.352.0d0)
     .   lihlte=(3.196d0-2.739d-1*tc+8.276d-3*tc**2)*1.0d-16/den
c heating
      llihlowh=-31.47d0+8.817d0*dlog10(tr)-4.144d0*(dlog10(tr))**2
     .   +0.8292d0*(dlog10(tr))**3-0.04996d0*(dlog10(tr))**4
      lihlowh=10**(llihlowh)
      lihlteh= 2.1d-11*(tr/400.0d0)**0.15*dexp(-1.9d3/tr)/den
       if(tr.le.352.0d0)
     .   lihlteh=(3.196d0-2.739d-1*tr+8.276d-3*tr**2)*1.0d-16/den
       dfdy(24,1)= dfdy(24,1)-(gamma-1.0d0)*7.243d15*y(23)/sumyi
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den*onoff(4)
c HeH+, assume it is LiH/10
       dfdy(24,1)= dfdy(24,1)-(gamma-1.0d0)*7.243d15*y(15)/sumyi
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den/1.0d1*onoff(3)
c ????
c      dfdy(24,1)= -(gamma-1.0d0)*7.243d15*y(18)
c     .   *(hdlte/(1.0d0+hdlte/hdlow)-hdlteh/(1.0d0+hdlteh/hdlowh))*den
c     .   +dfdy241
      dfdy(24,18)= -(gamma-1.0d0)*7.243d15*y(1)/sumyi*onoff(2)
     .   *(hdlte/(1.0d0+hdlte/hdlow)-hdlteh/(1.0d0+hdlteh/hdlowh))*den
c      dfdy(24,18)= -(gamma-1.0d0)*7.243d15
c     .  *2.1d-19*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
       dfdy(24,23)=-(gamma-1.0d0)*7.243d15*y(1)/sumyi*onoff(4)
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den
       dfdy(24,15)=-(gamma-1.0d0)*7.243d15*y(1)/sumyi*onoff(3)
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den/1.0d1
c      dfdy(24,23)= -(gamma-1.0d0)*7.243d15
c     .  *2.1d-12*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
c      dfdy(24,15)= -(gamma-1.0d0)*7.243d15
c     .  *5.0d-14*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
c      dfdy(24,17)= -(gamma-1.0d0)*7.243d15
c     .  *5.7d-14*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
c      dfdy(24,16)= -(gamma-1.0d0)*7.243d15
c     .  *2.4d-8*(((tc/1.0d4)**4.45*dexp(-tc/1.12d3)
c     .   - 1.0d-3*(tc/1.0d4)**2.3*dexp(-tc/1.8d3))
c     .   -((tr/1.0d4)**4.45*dexp(-tr/1.12d3)
c     .   - 1.0d-3*(tr/1.0d4)**2.3*dexp(-tr/1.8d3)))
c
c Other cooling functions
c dfdy(24,3) e-
c  Compton cooling (heating if tr>tc)
       dfdy(24,3)= -(gamma-1.0d0)/kb/sumyi*onoff(7)
     .    *1.017d-37*tr**4*(y(24)-tr)
c  Collisional ionization of H, He, He+
     .   -(gamma-1.0d0)/kb*y(1)*den/sumyi*onoff(9)
     . *1.27d-21*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-157809.1d0/tc)
     .   -(gamma-1.0d0)/kb*y(8)*den/sumyi*onoff(9)
     . *9.38d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-285335.4d0/tc)
     .   -(gamma-1.0d0)/kb*y(9)*den/sumyi*onoff(9)
     . *4.95d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-631515.d0/tc)
c Recombination of H, He, (He+ not yet included)
     .  -(gamma-1.0d0)/kb*y(2)*den/sumyi*onoff(9)
     . *8.7d-27*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
     .  -(gamma-1.0d0)/kb*y(9)*den/sumyi*onoff(9)
     . *1.55d-26*tc**0.3647
c       dydxhe+rec= -(gamma-1.0d0)/kb*y(3)*y(?)*den
c     . *3.48d-26*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
c Dielectronic recombination of He
     .  -(gamma-1.0d0)/kb*y(9)*den/sumyi*onoff(9)
     .  *1.24d-13*tc**-1.5*dexp(-4.7d5/tc)
     .  *(1.d0+0.3d0*dexp(-9.4d4/tc))
c Collisional excitation of H and He+ (He??)
     .  -(gamma-1.0d0)/kb*y(1)*den/sumyi*onoff(9)
     . *7.5d-19/(1.d0+(tc/1.d5)**0.5)*dexp(-118348.d0/tc)
     .  -(gamma-1.0d0)/kb*y(9)*den/sumyi*onoff(9)
     . *5.54d-17*tc**-0.397/(1.d0+(tc/1.d5)**0.5)*dexp(-473638.d0/tc)
c Bremsstrahlung of H+, He+ (add He++)
     .  -(gamma-1.0d0)/kb*(y(2)+y(9))*den/sumyi*onoff(8)
     . *1.42d-27*tc**0.5*(1.1d0+0.34d0
     . *dexp(-(5.5d0-dlog10(tc))**2/3.0d0))
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
     . -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)*onoff(6)
     .  *1.33d-08*(tc/300.0d0)**-0.912*exp(-5.58d4/tc)
c
c using Abel 1996 H2 cooling function
      dfdyh2=  -(gamma-1.0d0)*7.243d15*y(14)/sumyi*para*y(1)
     .   *4.37d-56*den**1.5*den**(-0.36*(dlog10(den))**0.5)
     .   *(16.05d0*tc**15.05*(tc**dlog10(tc))**(-1.9)
     .   - 0.82515d0*tc**15.05*(tc**dlog10(tc))**(-2.9)*dlog10(tc)
     .   *(tc**(dlog10(tc)-1.0d0)))
c
c using Galli and Palla 1998 HD cooling function
      dlte= hdlte
     . *(15.35716d0/tc-5.58513d0*2.0d0/tc*dlog10(tc)
     . +0.8561149d0*3.0d0/tc*(dlog10(tc))**2
     .  -1.75538d-2*4.0d0/tc*(dlog10(tc))**3)
c      dlow= hdlow
c     . *(30.13495d0/tc-14.4382d0*2.0d0/tc+3.1536d0*3.0d0/tc
c     .  -0.249489d-2*4.0d0/tc)
c      ddem=-(hdlte+hdlow)**-2*(dlte+dlow)*hdlte*hdlow
      dlow=(3.0d0*(4.4d-12+3.6d-13*tc**0.77)*dexp(-128.d0/tc)*128.d0
     .  *128.d0/tc**2
     .  + 3.0d0*3.6d-13/tc**0.23*0.77*dexp(-128.d0/tc)*128.d0
     .  + (5.0d0/3.0d0)*(4.1d-12+2.1e-13*tc**0.92)*dexp(-255.d0/tc)
     .  *255.0d0/tc**2
     .  + (5.0d0/3.0d0)*2.1e-13/tc**0.08*dexp(-255.d0/tc)
     .  *255.d0)*kb/1.27d0
      dfdyhd=-(gamma-1.0d0)*7.243d15*y(18)*y(1)*den/sumyi
     . *(hdlow**2*dlte+hdlte**2*dlow)/(hdlte+hdlow)**2
c      dfdyhd= -(gamma-1.0d0)*7.243d15*y(18)
c     .  *2.1d-19*((1/400.0d0)**(0.15)*tc**(-0.85)*0.15*dexp(-1.9d3/tc)
c     .  - (tc/400.0d0)**0.15*dexp(-1.9d3/tc)*1.9d3/tc**2)
c using Galli and Palla 1998 LiH cooling function
      dlte=lihlte*(0.15d0/tc+1.9d3/tc**2)
       if(tc.le.352.0d0) dlte=(-2.739d-1+1.6552d-2*tc)*1.0d-16
      dlow=lihlow*(8.817d0/tc-4.144d0*2.0d0/tc+0.892d0*3.0d0/tc
     .     -4.996d-2*4.0d0/tc)
      dfdylih=-(gamma-1.0d0)*7.243d15*y(23)*y(1)*den/sumyi
     . *(lihlow**2*dlte+lihlte**2*dlow)/(lihlte+lihlow)**2
c HeH+ assume LiH/10
      dfdyhehi=-(gamma-1.0d0)*7.243d15*y(15)*y(1)*den/sumyi
     . *(lihlow**2*dlte+lihlte**2*dlow)/(lihlte+lihlow)**2/1.0d1
c      dfdylih= -(gamma-1.0d0)*7.243d15*y(23)
c     .  *2.1d-12*((1/400.0d0)**(0.15)*tc**(-0.85)*0.15*dexp(-1.9d3/tc)
c     .  - (tc/400.0d0)**0.15*dexp(-1.9d3/tc)*1.9d3/tc**2)
c      dfdyhehi= -(gamma-1.0d0)*7.243d15*y(15)
c     .  *5.0d-14*((1/400.0d0)**(0.15)*tc**(-0.85)*0.15*dexp(-1.9d3/tc)
c     .  - (tc/400.0d0)**0.15*dexp(-1.9d3/tc)*1.9d3/tc**2)
c      dfdyhdi= -(gamma-1.0d0)*7.243d15*y(17)
c     .  *5.7d-14*((1/400.0d0)**(0.15)*tc**(-0.85)*0.15*dexp(-1.9d3/tc)
c     .  - (tc/400.0d0)**0.15*dexp(-1.9d3/tc)*1.9d3/tc**2)
c      dfdyh3i= -(gamma-1.0d0)*7.243d15*y(16)
c     .  *2.4d-8*((1.0d-4)**4.45*4.45*(tc/1.0d4)**3.45*dexp(-tc/1.120d3)
c     .  -1.0d-3*(1.d-4)**2.3*2.3*(tc/1.0d4)**1.3*dexp(-tc/1.8d3)
c     .  -(tr/1.0d4)**4.45*dexp(-tr/1.12d3)/1.12d3
c     .  +1.0d-3*(tr/1.0d4)**2.3*dexp(-tr/1.8d3)/1.8d3)
c Other cooling functions
c  Compton cooling (heating if tr>tc)
       dfdycomp= -(gamma-1.0d0)/kb*y(3)/sumyi
     .    *1.017d-37*tr**4
c  Collisional ionization of H, He, He+
       dfdyhion= -(gamma-1.0d0)/kb*y(3)*y(1)*den/sumyi
     . *1.27d-21*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-157809.1d0/tc)
     . *(157809.1d0/tc**2+0.5d0/tc-0.5d0/((1.0d0+(tc/1.d5)**0.5)
     . *(tc*1.d5)**0.5))
       dfdyheion= -(gamma-1.0d0)/kb*y(3)*y(8)*den/sumyi
     . *9.38d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-285335.4d0/tc)
     . *(285335.4d0/tc**2+0.5d0/tc-0.5d0/((1.0d0+(tc/1.d5)**0.5)
     . *(tc*1.d5)**0.5))
       dfdyhepion= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *4.95d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-631515.d0/tc)
     . *(631515.d0/tc**2+0.5d0/tc-0.5d0/((1.0d0+(tc/1.d5)**0.5)
     . *(tc*1.d5)**0.5))
c Recombination of H, He, (He+ not yet included)
       dfdyhrec= -(gamma-1.0d0)/kb*y(3)*y(2)*den/sumyi
     . *8.7d-27*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
     . *(0.5d0/tc-0.2d0/tc-0.7d0*(tc/1.d6)**0.7
     . /(1.d0+(tc/1.d6)**0.7)/tc)
       dfdyherec= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *1.55d-26*tc**0.3647*0.3647d0/tc
c       dfdyhe+rec= -(gamma-1.0d0)/kb*y(3)*y(?)*den/sumyi
c     . *3.48d-26*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
c     . *(0.5d0/tc-0.2d0/tc-0.7d0*(tc/1.d6)**0.7/(1.d0+(tc/1.d6)**0.7)/tc)
c Dielectronic recocombination of He
       dfdyhedr= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     .  *1.24d-13*tc**-1.5*dexp(-4.7d5/tc)
     .  *((1.d0+0.3d0*dexp(-9.4d4/tc))
     .  *(4.7d5/tc**2-1.5d0/tc) + 0.3d0*dexp(-9.4d4/tc)
     .  *9.4d4/tc**2)
c Collisional excitation of H and He+ (He??)
       dfdyhce= -(gamma-1.0d0)/kb*y(3)*y(1)*den/sumyi
     . *7.5d-19/(1.d0+(tc/1.d5)**0.5)*dexp(-118348.d0/tc)
     . *(118348.d0/tc**2-0.5d0/((1.0d0+(tc/1.d5)**0.5)
     . *(tc*1.d5)**0.5))
       dfdyhece= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *5.54d-17*tc**-0.397/(1.d0+(tc/1.d5)**0.5)*dexp(-473638.d0/tc)
     . *(473638.d0/tc**2-0.5d0/((1.0d0+(tc/1.d5)**0.5)
     . *(tc*1.d5)**0.5)-0.397d0/tc)
c Bremsstrahlung of H+, He+ (add He++)
       dfdybrem= -(gamma-1.0d0)/kb*y(3)*(y(2)+y(9))*den/sumyi
     . *1.42d-27*tc**0.5*((1.1d0+0.34d0
     . *dexp(-(5.5d0-dlog10(tc))**2/3.0d0))
     . *0.5d0/tc+0.34d0*dexp(-(5.5d0-dlog10(tc))**2/3.0d0)
     . *2.0d0*(5.5d0-dlog10(tc))*4.34294482d-1/3.0d0/tc)
c H2 formation cooling via H- and H2+
       dfdyh2form= -(gamma-1.0d0)/kb*y(1)*den/sumyi*everg
     . *3.00d0*1.5d-9*(tc/300.0d0)**-1.1*y(4)*-0.1d0/300.0d0
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
       dfdyh2diss= -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)
     . *(4.78d-14*(tc/300.0d0)**1.8*exp(-5.31d4/tc)*y(1)
     .   *(1.8d0/tc+5.31d4/tc**2)
     .  +1.47d-13*(tc/300.0d0)**1.39*exp(-1.31d5/tc)*y(8)
     .   *(1.39d0/tc+1.31d5/tc**2)
     .  +2.27d-19*(tc/300.0d0)**3.91*exp(-5.47e4/tc)*y(14)
     .   *(3.91d0/tc+5.47d4/tc**2)
     .  +1.33d-08*(tc/300.0d0)**-0.912*exp(-5.58d4/tc)*y(3)
     .   *(-9.12d-1/tc+5.58d4/tc**2))
c
        dfdy(24,24)=(gamma-1.0d0)*(den*mu/w0)**(0.5)*6.15583d-16
     .    + dfdyh2*onoff(1)+dfdyhd*onoff(2)+dfdyhehi*onoff(3)
     .    + dfdylih*onoff(4)+dfdyh2form*onoff(5)+dfdyh2diss*onoff(6)
     .    + dfdycomp*onoff(7)+dfdybrem*onoff(8)
     .    + (dfdyhion+dfdyheion+dfdyhepion+dfdyhrec
     .    + dfdyherec+dfdyhedr+dfdyhce+dfdyhece)*onoff(9)
c       do 10 j=1,n
c          write(6,*) 'dfdy= ',dfdy(1,1),dfdy(2,1),dfdy(4,1)
c 10    continue
      return
      END

      SUBROUTINE derivs(x,y,dydx)
      integer n,ihnu
      DOUBLE PRECISION x,y(*),dydx(*),den,nhd,w(3)
     . ,tc,ye(2),gamma,mu,xinp(1),tr,dydxh2,dydxhd,
     .  dydxlih,dydxhehi,dydxhdi, dydxh3i
     .  ,h2lte,lh2low,h2low,h2lteh,lh2lowh,h2lowh
     .  ,hdlte,lhdlow,hdlow,hdlteh,lhdlowh,hdlowh,
     .  lhdlte,lhdlteh,dydxcomp,kb, sumyi,para,
     .  rkh,cvh,cvl,crh,crl,crlt,tch2,thh2,xlt
     .  ,dydxhion,dydxheion,dydxhepion,dydxhrec
     .  ,dydxherec,dydxhedr,dydxhce,dydxhec,dydxbrem,
     .   dydxh2form,everg,dydxh2diss,w0,den0,onoff(9),
     .   lihlow,lihlte,lihlowh,lihlteh,llihlow,llihlowh
      parameter(kb=1.380658d-16,sumyi=1.08d0,para=1.0d0,
     .   everg=1.60217733d-12)
c      COMMON /prop/ w,den,nhd
      common /idata/ den,xinp,ihnu
      common /ndata/ n
      common /temp/ tc,gamma,mu
      common /prop/ tr,w0,den0
      common /onoffd/ onoff
      external fcn1,estate
c
c      call estate(2,tc,ye)
c      y(24)=y(1)*ye(1)
c      y(22)=y(5)*ye(2)

       call fcn1(n,y,dydx)
c don't need according to Shapiro & Kang ???
c       do 10 i=1,23
c          dydx(i)=dydx(i)-y(i)/den*dend
10     continue
       tc=y(24)
c
c ** Collision temperature equation **
c
c Total H2 cooling/heating (K/s)
c   using H2 cooling function of Abel 1995 thesis (ergs/s)
      dydxh2=   -(gamma-1.0d0)*7.243d15*y(14)/sumyi*para*y(1)
     .   *4.37d-56*den**1.5*den**(-0.36*(dlog10(den))**0.5)
     .   *(tc**16.05*(tc**dlog10(tc))**(-1.9)
     .   -tr**16.05*(tr**dlog10(tr))**(-1.9))
c
c HD cooling
c Fit to Lepp and Shull 1984 LTE (ergs/s) -> hdlte (ergs cm3/s)
      lhdlte=-35.6998d0+15.35716d0*dlog10(tc)-5.58513d0*(dlog10(tc))**2
     . +0.8561149d0*(dlog10(tc))**3-1.75538d-2*(dlog10(tc))**4
      hdlte=10**lhdlte/den
c Fit to Galli and Palla 1998 low density limit
c      lhdlow=-48.4421d0+30.13495d0*dlog10(tc)-14.4382d0*(dlog10(tc))**2
c     . +3.1536d0*(dlog10(tc))**3-0.249489d0*(dlog10(tc))**4
c      hdlow=10**(lhdlow)
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
      hdlow=(3.0d0*(4.4d-12+3.6d-13*tc**0.77)*dexp(-128.d0/tc)*128.d0
     .  + (5.0d0/3.0d0)*(4.1d-12+2.1e-13*tc**0.92)*dexp(-255.d0/tc)
     .  *255.d0)*kb/1.27d0
c heating
c Fit to Lepp and Shull 1984 LTE
      lhdlteh=-35.6998d0+15.35716d0*dlog10(tr)-5.58513d0*(dlog10(tr))**2
     . +0.8561149d0*(dlog10(tr))**3-1.75538d-2*(dlog10(tr))**4
      hdlteh=10**lhdlteh/den
c Galli and Palla 1998 low density limit
c      lhdlowh=-48.4421d0+30.13495d0*dlog10(tr)-14.4382d0*(dlog10(tr))**2
c      . +3.1536d0*(dlog10(tr))**3-0.249489d0*(dlog10(tr))**4
c      hdlowh=10**(lhdlowh)
c Galli and Palla 1998 low density limit (erg cm3 /s)
c  uses HD-He collisional data, so reduced by 1.27
      hdlowh=(3.0d0*(4.4d-12+3.6d-13*tr**0.77)*dexp(-128.d0/tr)*128.d0
     .  + (5.0d0/3.0d0)*(4.1d-12+2.1e-13*tr**0.92)*dexp(-255.d0/tr)
     .  *255.d0)*kb/1.27d0
c
      dydxhd=  -(gamma-1.0d0)*7.243d15*y(18)*y(1)/sumyi
c     . *2.1d-19*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
     .  *(hdlte/(1.0d0+hdlte/hdlow)-hdlteh/(1.0d0+hdlteh/hdlowh))*den
c
c LiH cooling function
c Galli and Palla 1998 low density limit (erg cm3 /s)
      llihlow=-31.47d0+8.817d0*dlog10(tc)-4.144d0*(dlog10(tc))**2
     .   +0.8292d0*(dlog10(tc))**3-0.04996d0*(dlog10(tc))**4
      lihlow=10**(llihlow)
c Fit to LTE of Dalgarno 1984 (erg/s) -> (ergs cm3/s
      lihlte= 2.1d-11*(tc/400.0d0)**0.15*dexp(-1.9d3/tc)/den
       if(tc.le.352.0d0)
     .   lihlte=(3.196d0-2.739d-1*tc+8.276d-3*tc**2)*1.0d-16/den
c heating
      llihlowh=-31.47d0+8.817d0*dlog10(tr)-4.144d0*(dlog10(tr))**2
     .   +0.8292d0*(dlog10(tr))**3-0.04996d0*(dlog10(tr))**4
      lihlowh=10**(llihlowh)
      lihlteh= 2.1d-11*(tr/400.0d0)**0.15*dexp(-1.9d3/tr)/den
       if(tr.le.352.0d0)
     .   lihlteh=(3.196d0-2.739d-1*tr+8.276d-3*tr**2)*1.0d-16/den
      dydxlih=  -(gamma-1.0d0)*7.243d15*y(23)*y(1)/sumyi
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den
c HeH+ cooling function, assume it is LiH/10
      dydxhehi=-(gamma-1.0d0)*7.243d15*y(15)*y(1)/sumyi
     .  *(lihlte/(1.0d0+lihlte/lihlow)-lihlteh/
     .   (1.0d0+lihlteh/lihlowh))*den/1.0d1
c      dydxhehi=  -(gamma-1.0d0)*7.243d15*y(15)
c     . *5.0d-14*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
c      dydxhdi=  -(gamma-1.0d0)*7.243d15*y(17)
c     . *5.7d-14*((tc/400.0d0)**0.15*dexp(-1.9d3/tc)
c     . -(tr/400.0d0)**0.15*dexp(-1.9d3/tr))
c      dydxh3i= -(gamma-1.0d0)*7.243d15*y(16)
c     .  *2.4d-8*(((tc/1.0d4)**4.45*dexp(-tc/1.12d3)
c     .   - 1.0d-3*(tc/1.0d4)**2.3*dexp(-tc/1.8d3))
c     .   -((tr/1.0d4)**4.45*dexp(-tr/1.12d3)
c     .   - 1.0d-3*(tr/1.0d4)**2.3*dexp(-tr/1.8d3)))
c
c Other cooling processes
c  Compton cooling (heating if tr>tc)
       dydxcomp= -(gamma-1.0d0)/kb*y(3)
     .    *1.017d-37*tr**4*(y(24)-tr)/sumyi
c  Collisional ionization of H, He, He+
       dydxhion= -(gamma-1.0d0)/kb*y(3)*y(1)*den/sumyi
     . *1.27d-21*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-157809.1d0/tc)
       dydxheion= -(gamma-1.0d0)/kb*y(3)*y(8)*den/sumyi
     . *9.38d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-285335.4d0/tc)
       dydxhepion= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *4.95d-22*tc**0.5/(1.d0+(tc/1.d5)**0.5)*dexp(-631515.d0/tc)
c Recombination of H, He, (He+ not yet included)
       dydxhrec= -(gamma-1.0d0)/kb*y(3)*y(2)*den/sumyi
     . *8.7d-27*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
       dydxherec= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *1.55d-26*tc**0.3647
c       dydxhe+rec= -(gamma-1.0d0)/kb*y(3)*y(?)*den/sumyi
c     . *3.48d-26*tc**0.5/(1.d0+(tc/1.d6)**0.7)/(tc/1.0d3)**0.2
c Dielectronic recocombination of He
       dydxhedr= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     .  *1.24d-13*tc**-1.5*dexp(-4.7d5/tc)
     .  *(1.d0+0.3d0*dexp(-9.4d4/tc))
c Collisional excitation of H and He+ (He??)
       dydxhce= -(gamma-1.0d0)/kb*y(3)*y(1)*den/sumyi
     . *7.5d-19/(1.d0+(tc/1.d5)**0.5)*dexp(-118348.d0/tc)
       dydxhece= -(gamma-1.0d0)/kb*y(3)*y(9)*den/sumyi
     . *5.54d-17*tc**-0.397/(1.d0+(tc/1.d5)**0.5)*dexp(-473638.d0/tc)
c Bremsstrahlung of H+, He+ (add He++)
       dydxbrem= -(gamma-1.0d0)/kb*y(3)*(y(2)+y(9))*den/sumyi
     . *1.42d-27*tc**0.5*(1.1d0+0.34d0
     . *dexp(-(5.5d0-dlog10(tc))**2/3.0d0))
c H2 formation cooling via H- and H2+
       dydxh2form= -(gamma-1.0d0)/kb*y(1)*den/sumyi*everg
     . *(3.00d0*1.5d-9*(tc/300.0d0)**-0.1*y(4)
     . + 1.83d0*6.4d-10*y(13))
c H2 dissociation cooling via H, He, H2, and e, assuming 4 eV
       dydxh2diss= -(gamma-1.0d0)/kb*den/sumyi*everg*4.0d0*y(14)
     . *(4.78d-14*(tc/300.0d0)**1.8*exp(-5.31d4/tc)*y(1)
     .  +1.47d-13*(tc/300.0d0)**1.39*exp(-1.31d5/tc)*y(8)
     .  +2.27d-19*(tc/300.0d0)**3.91*exp(-5.47e4/tc)*y(14)
     .  +1.33d-08*(tc/300.0d0)**-0.912*exp(-5.58d4/tc)*y(3))
c sum
         dydx(24)=(gamma-1.0d0)*(den*mu/w0)**(0.5)*y(24)*6.15583d-16
     .     +dydxh2*onoff(1)+dydxhd*onoff(2)+dydxhehi*onoff(3)
     .     +dydxlih*onoff(4)+dydxh2form*onoff(5)
     .     +dydxh2diss*onoff(6)+dydxcomp*onoff(7)
     .     +dydxbrem*onoff(8)
     .  +(dydxhion+dydxheion+dydxhepion+dydxhrec+dydxherec+dydxhedr
     .  +dydxhce+dydxhece)*onoff(9)
c       write(*,*) 'dydx=',dydx(1),dydx(2),dydx(4)
c       write(*,*) 'y= ',y(1),y(2),y(4)
c
c      write(*,*) 'check', dydx(10), a(24,tr)
c      print *, dydx(5),dydx(6),'tr=',tr,'tc=',tc
      return
      END
