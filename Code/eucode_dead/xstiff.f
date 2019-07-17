c Calls recfast for e,H,H+,He,He+
      PROGRAM xstiff
C      driver for routine stiff, d chemistry
      INTEGER KMAXX,NMAX,nm,nz,neq,ihnu,nq,ne,nnz
      PARAMETER (KMAXX=200,NMAX=100,nm=21,neq=24
     .   ,ne=2,nstars=var0)
      DOUBLE PRECISION dxsav,eps,hstart,x1,x2,y(neq),xp,yp,
     .   omegab,h100,z,tr,den,w0,dtdz,w(3),xs,tc,tr1,ye(ne),
     .   zm(nm),tm(nm),tm2(nm),zstep0,zstep,zstart,xinp(8),
     .   yrd,yr1,yr2,yhe,dhp,lihp,zstep1,
     .   wdm,ys(5),yr(3), zs,tstar,zstep2,
     .   zinitial,zfinal,tnow,wvac,h0,nnow,zeq
     .   ,OmegaT,OmegaL,OmegaK,OmegaC
     .   ,a,b,weight

      INTEGER kmax,kount,nbad,nok
      CHARACTER*15 filename
      CHARACTER*16 filetwo
      CHARACTER*7 filethree
      LOGICAL file_exists
      COMMON /path/ kmax,kount,dxsav,xp(KMAXX),yp(NMAX,KMAXX)
      common /idata/ den,xinp,ihnu
      common /ndata/ nq
      common /prop/ tr
      common /temp/ tc
      common /exci/ yrd
      common /seager/ ys
      common/Cosmo/tnow,hO,nnow,zeq,OmegaT,OmegaL,OmegaK
      common/zLIST/zinitial,zfinal,nnz
      common/loop/m
      EXTERNAL input,rates,derivs,stifbs,tspln,splint,estate,derivsp,
     . recfast
c
c read in model parameters
c
      open(unit=19, file="point.data")
      read(19,*) weight
      close(unit=12)
c
      if (weight.ne.0.0) then
      do m=1,nstars
c
      write(filetwo, fmt='(a,i3.3,a)')'./sabund/',m-1,'.txt'
      INQUIRE(FILE=filetwo, EXIST=file_exists)
      if (file_exists) then
      write(filename, fmt='(a,i3.3,a)')'./model/',m-1,'.txt'
      open(unit=10, file=filename)
      read (10,*) omegab, h100, w0
      read (10,*) yhe, dhp, lihp
      read (10,*) zstep0, zstep1, nz, zstart, xs
      close(unit=10)
c
       eps=1.0d-3
       hstart=10.0d0
       nq=neq
       wdm=w0-omegab
       OmegaC=wdm
       wvac=1.d0-w0
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**-1.5d0
       OmegaT=OmegaC+omegab
       OmegaK=1.d0-OmegaT-OmegaL
       tnow=2.728d0
       h0=0.0d0
       nnow=0.0d0
       zeq=0.0d0
       zstep=zstep0
c          write(17,998) tnow,hO,nnow,zeq,OmegaT,OmegaL,OmegaK
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
         if (z .gt. zstart) goto 112
	 yrd=yr1/yr2
         ye(1)=yrd
	 ye(2)=yrd
       z=zstart
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**-1.5d0
       OmegaK=1.d0-OmegaT-OmegaL
       if(z.lt.500.0d0) call estate(ne,tc,ye)
c      n(H)=w(1), n(H 2s,2p)=w(3) and n(e)=n(H+)=w(2)
c
c      write(*,*) z, yrd
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
c      if (z .ge. 1000.0d0) then
         tc=2.728d0*(1.0d0+z)
      else
c         call splint(zm,tm,tm2,nm,z,tr1)
c         tc=tr1
          tc=y(24)
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
       open(unit=14, file=filetwo)

       read (14,*) ys(1)
       read (14,*) ys(2)
       read (14,*) ys(3)
       read (14,*) ys(4)
       read (14,*) ys(5)

       read (14,*) y(1)
       read (14,*) y(2)
       read (14,*) y(3)
       read (14,*) y(4)
       read (14,*) y(5)
       read (14,*) y(6)
       read (14,*) y(7)
       read (14,*) y(8)
       read (14,*) y(9)
       read (14,*) y(10)
       read (14,*) y(11)
       read (14,*) y(12)
       read (14,*) y(13)
       read (14,*) y(14)
       read (14,*) y(15)
       read (14,*) y(16)
       read (14,*) y(17)
       read (14,*) y(18)
       read (14,*) y(19)
       read (14,*) y(20)
       read (14,*) y(21)
       read (14,*) y(22)
       close(unit=14)
       else
        goto 222
      endif

c
c      write(6,*) 'time', 'H','p','e','H-'
c      write(6,*) x1,y
c      print *, 'Star ',m
      do it1=1,nz
       OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**-1.5d0
       OmegaK=1.d0-OmegaT-OmegaL
c time to redshift conversion
      if(z.le.30.0d0) zstep=zstep1
      if(z.le.15.0d0) zstep=zstep0
      dtdz=3.15d7*0.98d10/h100/(1.0d0+z)**2
     . /dsqrt(1.0d0+w0*z)
      x2=x1+dtdz*zstep
      den=1.123d-5*(1.0d0-yhe)*omegab*h100**2*(1.0d0+z)**3
c black-body radiation temperature
      tr=2.728d0*(1.0d0+z)
c  use matter temperature for z < 600 for models 1, 3, and 5
c z < 1000 for models 2 and 4
      if (z .ge. 600.0d0) then
c      if (z .ge. 1000.0d0) then
         tc=2.728d0*(1.0d0+z)
      else
c         call splint(zm,tm,tm2,nm,z,tr1)
c         tc=yr(3)
          tc=y(24)
        end if
c      end if
c
c  calculate rates
      call rates(tr,tc,z)
c
c      print *, 'tr=',tr,'tc=',tc,'z=',z
      call odeint(y,neq,x1,x2,eps,hstart,0.d0,nok,nbad,derivs,stifbs)
c     write(*,'(/1x,a,t30,i3)') 'Successful steps:',nok
c     write(*,'(1x,a,t30,i3)') 'Bad steps:',nbad
c      print*,'ye=1',ye(1),ye(2)
       write(filethree, fmt='(i3.3,a)') m-1,'.txt'
       open(unit=16, file=filethree)
       write(16,998) z,(ys(i),i=1,5),(y(i),i=1,neq),
     .              ye(1)*ys(1),ye(2)*y(2)
c       write(50,998) z,(y(i),i=neq-3,neq)
       call derivsp(z,y)
c      print*, tr,tc
      write(35,998) z,tr,tc,den,x2
c      print*, tr,tc
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
c  call recfast to get H,H+,e,He,He+
c          write(17,998) hO,nnow,zeq,OmegaT,zinitial,zfinal
      call recfast(omegab,wdm,h100*100.0d0,yhe,yr)

         ys(1)=(1.0d0-yr(1))
         ys(2)=yr(1)
         ys(3)=yr(1)+yr(2)*yhe/(4.0d0*(1.0d0-yhe))
         ys(4)=(1.0d0-yr(2))*yhe/(4.0d0*(1.0d0-yhe))
         ys(5)=yr(2)*yhe/(4.0d0*(1.0d0-yhe))

c        print*, tr,tc
          write(17,998) z, (ys(i),i=1,5),
     .    yr(1),yr(2),yr(3),tc,tr
c         print*, tr,tc
          write(99,*) OmegaL
       tc=yr(3)
c      write(*,*) zr, ye(1), ye(2)
      enddo
      close(unit=16)
222   enddo
      endif
999   write(*,*) 'NORMAL COMPLETION'
998   format(30(E11.4E3,1x))
      STOP
      END
