c Calls recfast for e,H,H+,He,He+
      PROGRAM xstiff
C      driver for routine stiff, d chemistry
      INTEGER KMAXX,NMAX,nm,nz,neq,ihnu,nq,ne,nnz
      PARAMETER (KMAXX=200,NMAX=100,nm=21,neq=27
     .   ,ne=2)
      DOUBLE PRECISION dxsav,eps,hstart,x1,x2,y(neq),xp,yp,
     .   omegab,h100,z,tr,den,w0,dtdz,w(3),xs,tc,tr1,ye(ne),
     .   zm(nm),tm(nm),tm2(nm),zstep0,zstep,zstart,xinp(8),
     .   yrd,yr1,yr2,yhe,dhp,lihp,zstep1,
     .   wdm,ys(5),yr(3),
     .   zinitial,zfinal,tnow,wvac,h0,nnow,zeq
     .   ,OmegaT,OmegaL,OmegaK,OmegaC
     .   ,a,b,deplete
      INTEGER kmax,kount,nbad,nok
      COMMON /path/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
      common /idata/ den,xinp,ihnu
      common /ndata/ nq
      common /prop/ tr
      common /temp/ tc
      common /exci/ yrd
      common /seager/ ys
      common/Cosmo/tnow,hO,nnow,zeq,OmegaT,OmegaL,OmegaK
      common/zLIST/zinitial,zfinal,nnz
      common/star/deplete
      EXTERNAL input,rates,derivs,stifbs,tspln,splint,estate,derivsp,
     . recfast
c
c read in model parameters
c
      open(unit=10, file='model.data')
      read (10,*) omegab, h100, w0
      read (10,*) yhe, dhp, lihp
      read (10,*) zstep0, zstep1, nz, zstart, xs
      close(unit=10)
c
      deplete=1.0d0
c
       eps=1.0d-5
       hstart=100.0d0
       nq=neq
       wdm=w0-omegab
       OmegaC=wdm
       wvac=1.d0-w0
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**(-1.5d0)
       OmegaT=OmegaC+omegab
       OmegaK=1.d0-OmegaT-OmegaL
       tnow=2.728d0
       h0=0.0d0
       nnow=0.0d0
       zeq=0.0d0
       zstep=zstep0
c
c  spline fit LiH+ equilibrium constants k1 and k2
      call kspln
      call kspln1
      call kspln2
c
c  spline fit to matter temperature
      call tspln(nm,zm,tm,tm2)
c
c look-up Dell'Antonio-Rybicki n(H) and n(e)
       open(unit=15, file='hh+eIIIext.data')
112      read(15,*) z,yr1,yr2
         if (z .ne. zstart) goto 112
	 yrd=yr1/yr2
         ye(1)=yrd
	 ye(2)=yrd
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**(-1.5d0)
       OmegaK=1.d0-OmegaT-OmegaL
       if(z.lt.500.0d0) call estate(ne,tc,ye)
c
      kmax=KMAXX
      x1=xs
      den=1.123d-5*(1.0d0-yhe)*omegab*h100**2*(1.0d0+z)**3
c
c read in list of species (species.data) and reactions (rates.data)
      call input
c
c black-body radiation temperature
         tr=2.728d0*(1.0d0+z)
c  use matter temperature for z < 600 for Models 1, 3, and 5
c z< 1000 for models 2 and 4
      if (z .ge. 600.0d0) then
         tc=2.728d0*(1.0d0+z)
      else
         call splint(zm,tm,tm2,nm,z,tr1)
         tc=tr1
      end if
c
c  call recfast to get H,H+,e,He,He+
          zinitial=z
          zfinal=z-zstep
          yr(1)=1.0d0
          yr(2)=1.0d0
          yr(3)=tc
      call recfast(omegab,wdm,h100*100.0d0,yhe,yr)
         ys(1)=(1.0d0-yr(1))
         ys(2)=yr(1)
         ys(3)=yr(1)+yr(2)*yhe/(4.0d0*(1.0d0-yhe))
         ys(4)=(1.0d0-yr(2))*yhe/(4.0d0*(1.0d0-yhe))
         ys(5)=yr(2)*yhe/(4.0d0*(1.0d0-yhe))

c  calculate rates
      call rates(tr,tc,z)
c
c starting abundances H-,D,D+,D-,Li,Li+, Li-,           7
c                     H2+,H2,HeH+,H3+,HD+,HD,HeD+,H2D+, 8
c                     He2+,LiH+,LiH                     3
c                     H*,D*,H,H+,e,He,He+                7
c
       open(unit=14, file='initabundance.data')
       do j=1,neq
           read (14,*) y(j)
       end do
       close(unit=14)
c
c      write(6,*) 'time', 'H','p','e','H-'
c      write(6,*) x1,y
      do it1=1,nz
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**(-1.5d0)
       OmegaK=1.d0-OmegaT-OmegaL
c time to redshift conversion
      if(z.lt.30.0d0) zstep=zstep1
      if(z.lt.20.0d0) zstep=zstep0
      dtdz=3.15d7*0.98d10/h100/(1.0d0+z)**2
     . /dsqrt(1.0d0+w0*z)
      x2=x1+dtdz*zstep
      den=1.123d-5*(1.0d0-yhe)*omegab*h100**2*(1.0d0+z)**3
c black-body radiation temperature
      tr=2.728d0*(1.0d0+z)
c  use matter temperature for z < 600 for models 1, 3, and 5
c z < 1000 for models 2 and 4
      if (z .ge. 600.0d0) then
         tc=2.728d0*(1.0d0+z)
      else
         call splint(zm,tm,tm2,nm,z,tr1)
         tc=tr1
         tr=2.728d0*(1.0d0+z)
      end if
c
c  calculate rates
      call rates(tr,tc,z)
c
c      print *, 'tr=',tr,'tc=',tc,'z=',z
      call odeint(y,neq,x1,x2,eps,hstart,0.d0,nok,nbad,derivs,stifbs)
       write(16,998) z,(y(i),i=1,neq),
     .              ye(1)*y(1),ye(2)*y(2)
       call derivsp(z,y)
      write(35,998) z,tr,tc,den,x2
      x1=x2
      z=z-zstep
      zinitial=z
      zfinal=zinitial-zstep
      call estate(ne,tc,ye)
      if(z.ge.500.0d0) then
111      read(15,*) zr,yr1,yr2
         if (zr .ne. z) goto 111
	 yrd=yr1/yr2
         ye(1)=yrd
         ye(2)=yrd
      endif
          write(17,998) z,
     .    y(1),y(2),y(3),tc,tr
       tc=tr1
      enddo
c
      open(unit=20, file='halfabund.data')
      do j=1,neq
        write (20,*) y(j)
      end do
      close(unit=20)
c

999   write(*,*) 'NORMAL COMPLETION'
998   format(30(E11.4E3,1x))
      STOP
      END
