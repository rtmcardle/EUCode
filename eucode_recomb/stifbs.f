      SUBROUTINE stifbs(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(nv),
     .   y(nv),yscal(nv),SAFE1,SAFE2,
     *REDMAX,REDMIN,TINY,SCALMX
      EXTERNAL derivs
      PARAMETER (NMAX=100,KMAXX=7,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     *REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0)
CU    USES derivs,jacobn,simpr,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nvold,nseq(IMAX)
      DOUBLE PRECISION eps1,epsold,errmax,fact,h,red,
     *scale,work,wrkmin,xest,xnew,
     *a(IMAX),alf(KMAXX,KMAXX),dfdx(NMAX),dfdy(NMAX,NMAX),err(KMAXX),
     *yerr(NMAX),ysav(NMAX),yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,nvold,xnew
      DATA first/.true./,epsold/-1.d0/,nvold/-1/
      DATA nseq /2,6,10,14,22,34,50,70/
      if(eps.ne.epsold.or.nv.ne.nvold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        nvold=nv
        a(1)=nv+a(1)
        do 14 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
14      continue
        do 15 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
15      continue
1       kmax=kopt
      endif
      h=htry
      do 16 i=1,nv
        ysav(i)=y(i)
16    continue
      call jacobn(x,y,dfdx,dfdy,nv,nmax)
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 18 k=1,kmax
        xnew=x+h
c        if(xnew.eq.x)pause 'stepsize underflow in stifbs'
        if(xnew.eq.x) STOP
        call simpr(ysav,dydx,dfdx,dfdy,nmax,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 17 i=1,nv
            errmax=dmax1(errmax,dabs(yerr(i)/yscal(i)))
17        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.d0)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
18    continue
3     red=dmin1(red,REDMIN)
      red=dmax1(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 19 kk=1,km
        fact=dmax1(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
19    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=dmax1(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END

