c+++++++++++++++++++++++++++++++++++++++++++++++++++++++
c subroutine to read species and rates to be included in
c calculation from species.data and rates.data.
c
c output is n=number of species and
c fills common blocks io,idata,rdata
c
c----------------------------------------
c written by S. Lepp ca 1984
c modified by P. C. Stancil, 6-25-96
c----------------------------------------
c name  character array of gas particle names, could
c       be product or reactant, read from species.data
c
c nminp character array of particle names whose abundances
c       are not to be solved for, in this case photons,
c       H*, and D*
c
c xkp   array of fit parameters for each reaction
c       with the definition: xkp(i,1)=a(cm^3/s),
c       xkp(i,2)=b, xkp(i,3)=c(K), and xkp(i,4)=d
c       where for collision reactions the absolute
c       value of d is a branching ratio and
c
c       if d>0   alpha(T)=a*|d|*(T/300)**b*exp(-c/T)
c       if d<0   alpha(T)=a*|d|*(T/300)**b*exp(-T/c).
c
c       For photo-destruction reactions the value of
c       d corresponds to a look-up table in subrountine
c       ke.f of equilibrium constants:
c
c       if d>0   alpha(T)=a*(T/300)**b*exp(-c/T)*k(|d|)
c       if d<0   alpha(T)=a*(T/300)**b*exp(-T/c)*k(|d|).
c
c rch   temporary storage character array of reactant
c       and product names for a reaction
c
c r     2D integer array containing index number for
c       reactants and products of each reaction.
c       Positive integer for gas particles and -1 for
c       photons. Array filled by call to integer
c       function indx and comparison to particle names
c       from species.data.
c
c ihnu  index number, -1, for photon
c
c idim  number of non-solve particles: photon, H*, and D*,He,He+
c
c nreacs total number of reactions. Determine here from
c        reading rates.data
c
c den, xk, xinp defined in rates.f
c -------------------------------------------------------
      subroutine input
      implicit double precision (a-h,o-z)
c
      integer maxr,r,maxx
c
      parameter (idim = 8)
      parameter (maxr =260)
      parameter (maxx =100)
c
      character*7 name(maxx),nminp(idim),rch(6),symbol(5)
      character*1 flag
c
      dimension r(maxr,6),xk(maxr),xkp(maxr,4),xinp(idim)
      common /io/ name,nminp
      common /idata/ den,xinp,ihnu
      common /rdata/r,xk,xkp,nreacs
      external indx
c
      open(unit=7,file='rates.data',status='old')
      open(unit=8,file='species.data',status='old')
c
c read gas particle names to be included in calculation
c
      do 13 i=1,maxx
  13  name(i)='*****'
      symbol(1)='  +'
      symbol(2)='  ->'
      i = 1
501   read(8,503) name(i)
503   format(a7)
      if(name(i).eq.'end    ') goto 502
      i=i+1
      goto 501
502   name(i)='       '
      n=i-1
c      print*,'n=',n
c
c assign non-solve particle parameters
c
      nminp(1)='photon '
      ihnu=-1
      nminp(2)='H*     '
      nminp(3)='D*     '
      nminp(4)='H      '
      nminp(5)='H+     '
      nminp(6)='e      '
      nminp(7)='He     '
      nminp(8)='He+    '
c
c read reactions and rate coefficient fit parameters
c
      i = 1
    1   read(7,990)iflag, (rch(j), j=1,6),(xkp(i,j),j=1,4)
  990   format(i4,1x,6(a7),4e10.2)
        xkp(i,3)=-xkp(i,3)
        if (iflag.eq.9999) goto 5
c
c assign reactant and product indexes
c
        do 103 j=1,6
          r(i,j)=indx(rch(j),name,n,nminp,idim)
  103   continue
c
  4   i=i+1
      if (i.le.maxr) goto 1
      write(6,*) i
      stop 'too many reacs increase parameter maxr'
    5 nreacs= i-1
c      print*,'nreacs=',nreacs
       close(unit=7)
       close(unit=8)
       return
       end
