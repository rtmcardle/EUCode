c Calls recfast for e,H,H+,He,He+
      PROGRAM xstiff
C      driver for routine stiff, d chemistry
      INTEGER KMAXX,NMAX,nm,nz,neq,ihnu,nq,ne,nnz,nstars,nlines,i,
     .   n,counter
      PARAMETER (KMAXX=200,NMAX=100,nm=21,neq=27
     .   ,ne=2,nstars=4,nlines=1576)
      DOUBLE PRECISION dxsav,eps,hstart,x1,x2,y(neq),xp,yp,
     .   omegab,h100,z,tr,den,w0,dtdz,w(3),xs,tc,tr1,ye(ne),
     .   zm(nm),tm(nm),tm2(nm),zstep0,zstep,zstart,xinp(8),
     .   yrd,yr1,yr2,yhe,dhp,lihp,zstep1,wdm,ys(5),yr(3),
     .   zinitial,zfinal,tnow,wvac,h0,nnow,zeq,
     .   OmegaT,OmegaL,OmegaK,OmegaC,a,pi,zn,weight,deplete,
     .   mass(nstars),lifetime(nstars),teff(nstars),srad(nstars),
     .   Q(nstars),zbirth(nstars),r(nstars),theta(nstars),phi(nstars),
     .   dist(nstars),zed(nlines,nstars),frad(nlines,nstars)
      INTEGER kmax,kount,nbad,nok
      CHARACTER*3 string
      CHARACTER*16 filetwo
      CHARACTER*15 filethree
      CHARACTER*13 fileone
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
      open(unit=12, file="stars.data")
      open(unit=19, file="point.data")
      read(19,*) weight
      do n=1,nstars
        read(12,*) mass(n),lifetime(n),teff(n),srad(n),Q(n),zbirth(n),
     .        r(n),theta(n),phi(n)
        read(19,*) dist(n)
        write(string,'(I3.3)') n-1
        open(unit=13, file="zrd_"//string//".txt")
        do m=1,nlines
          read(13,*) zed(m,n),frad(m,n)
        end do
        close(unit=13)
      end do
      close(unit=12)
      close(unit=19)
c
      if (weight.ne.0.0) then
      do n=1,nstars
c
      open(10, file='model.data')
      read (10,*) omegab, h100, w0
      read (10,*) yhe, dhp, lihp
      read (10,*) zstep0, zstep1, nz, zstart, xs
      close(unit=10)
c
      pi=3.14159265d0
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
      deplete=1.0d0
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
c ADD 'IF < 500 PULL LINE WHERE Z = 500'
      if (zstart .gt. 500) then
        open(unit=15, file='hh+eIIIext.data')
112     read(15,*) z,yr1,yr2
        if (z .ne. zstart) goto 112
        else
          open(unit=18, file='hh+eIIIext.data')
113       read(18,*) zn,yr1,yr2
          if (zn .ne. 500.0) goto 113
          z=zstart
          close(unit=18)
      end if
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
      do l=1,neq
        read (14,*) y(l)
      end do
      close(unit=14)


        if (dist(n).ne.0.0d0) then
c        print *, 'Star ',n
        do it1=1,nz
          m=it1
          OmegaL=wvac+(1.5d-4)*(1/(1.d0+z))**(-1.5d0)
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
            tc=2.728d0*(1.0d0+z)
          else
            call splint(zm,tm,tm2,nm,z,tr1)
            tc=tr1
            if (frad(m,n).ge.dist(n)) then
              tr=teff(n)
              deplete=(srad(n)/dist(n))**2
c              print *, n,' in range. Deplete: ', deplete, ' tr=', tr
            else
              tr=2.728d0*(1.0d0+z)
              if ((z.lt.zbirth(n)).AND.(frad(m,n).eq.0.0d0)) then
c               write(*,*) 'MADE IT TO 1'
                write(filetwo, fmt='(a,i3.3,a)')'./sabund/',n-1,'.txt'
                open(unit=20, file=filetwo)
                do j=1,neq
c                 write (*,*) y(j)
                  write (20,*) y(j)
                end do
                close(unit=20)
                write(filethree, fmt='(a,i3.3,a)')'./model/',n-1,'.txt'
                open(unit=21, file=filethree)
                write (21,*) omegab, h100, w0
                write (21,*) yhe, ehp, lihp
                write (21,*) zstep0, zstep1, nz+1-m, z, xs
                close(unit=21)
                goto 222
              end if
c              print *, n,' not in range. Deplete: ', deplete
            end if
          end if
c
c  calculate rates
        call rates(tr,tc,z)
c
c      print *, 'tr=',tr,'tc=',tc,'z=',z
        call odeint(y,neq,x1,x2,eps,hstart,0.d0,nok,nbad,derivs,stifbs)
        write(fileone, fmt='(a,i3.3,a)')'abund.',n-1,'.txt'
        open(unit=16, file=fileone, access='append')
        write(16,998) z,(y(i),i=1,neq),
     .             ye(1)*y(1),ye(2)*y(2)
        call derivsp(z,y)
        write(35,998) z,tr,tc,den,x2
        x1=x2
        z=z-zstep
        zinitial=z
        zfinal=zinitial-zstep
        call estate(ne,tc,ye)
        if(z.ge.500.0d0) then
111        read(15,*) zr,yr1,yr2
          if (zr .ne. z) goto 111
      	 yrd=yr1/yr2
          ye(1)=yrd
          ye(2)=yrd
          endif
            write(17,998) z,
     .     y(1),y(2),y(3),tc,tr
        tc=tr1
        enddo
c       write(*,*) 'MADE IT TO 3'
222     close(unit=16)
        end if
      enddo
      end if
999   write(*,*) 'NORMAL COMPLETION'
998   format(30(E11.4E3,1x))
      STOP
      END
