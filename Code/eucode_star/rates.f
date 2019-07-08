c 
c 6 june 01 PCS data from recfast added (H,H+,e,He,He+)
c 06-04-01 PCS hacking on li rr rate
c 06-04-01 PCS hacking on li photoionization rate
c adds input of Seager H and He data
c includes fix for liH radiative association rate 15 Jan 1998
c++++++++++++++++++++++++++++++++++++++++++++++
c subroutine to calculate all rate coefficients
c at a specific radiation and matter temperature.
c
c rate coefficients are placed in common block
c rdata
c
c----------------------------------------
c written by S. Lepp ca 1984
c modified by P. C. Stancil, 6-25-96
c----------------------------------------
c tr   radiation temperature (K)
c
c tc   matter or collision temperature (K)
c
c xk   array of radiative coefficients for each
c      reaction
c
c ki   integer corresponding to equilibrium
c      constant in function ke 
c
c den  total number density (1/cm^3). Must be calculated
c      in a driver and supplied in common block idata
c
c xinp photon 'fractional abundance' and H* and D* fractional abundance
c 
c r, nreacs, name, nminp, xkp, idim, 
c ihnu and the form of the rate coefficient fits
c are defined in input.f
c
c------------------------------------------------------     
      subroutine rates(tr,tc,z)
      implicit double precision (a-h,o-z)
      double precision ke,kee,tc,tr,z,yrd,nhd
c
      integer dim,maxr,r,ki,maxx,nreacs
c
c  dim2 is dim + num conservation equations
c
      parameter (idim = 8)
      parameter (maxr =260)
      parameter (maxx =100)
c
      dimension r(maxr,6)
      character*7 name(maxx),nminp(idim),rch(6),symbol(5)
      character*1 flag
c
      double precision xinp(idim),xk(maxr),xkp(maxr,4),xtmp(3),
     .  ys(5),deplete
      common /io/ name,nminp
      common /idata/ den,xinp,ihnu
      common /rdata/r,xk,xkp,nreacs
      common /exci/ yrd
      common /seager/ ys
      common /star/ deplete
      external ke,estate
c
c   determine photon and H* abundances
c  
      xinp(1)=1.0d0/den
      call estate(2,tc,xtmp)
      if(z.ge.500.0d0) then
         xtmp(1)=yrd
	 xtmp(2)=yrd
      endif
      xinp(2)=xtmp(1)
      xinp(3)=xtmp(2)
      xinp(4)=ys(1)
      xinp(5)=ys(2)
      xinp(6)=ys(3)
      xinp(7)=ys(4)
      xinp(8)=ys(5)
c
      do 4 i=1,nreacs
          xk(i)=0.0d0
c
c calculate reaction rates
c
          if ((r(i,1).eq.ihnu).or.(r(i,2).eq.ihnu) )goto 3
c if one of the reactants is a photon, skip down to 3
c
      xk(i)=xkp(i,1)
      if(xkp(i,2).ne.0.0d0) xk(i)=xk(i)*(tc/300.0d0)**xkp(i,2)
      if(xkp(i,3).ne.0.0d0) then
        if(xkp(i,4).gt.0.0d0) then
           xk(i)=xk(i)*exp(xkp(i,3)/tc)
        elseif(xkp(i,4).lt.0.0d0) then 
           xk(i)=xk(i)*exp(tc/xkp(i,3))
        else
           print*,'problem with rate exp term'
        endif
      endif
c fix for LiH radiative assocition rate
      if((i.eq.120).and.(tc.lt.380.0d0)) then
	xk(i)=3.74d-20*(tc/300.0d0)**0.11*dexp(-tc/2.0d3)
	goto 4
      endif
c
      if(xkp(i,4).ne.0.0d0) xk(i)=xk(i)*dabs(xkp(i,4))
c Verner & Ferland (1996) Li RR rate
      if(i.eq.23) xk(i)=1.036d-11*((tc/1.077d2)**0.5
     .  *(1.0d0+(tc/1.077d2)**0.5)**0.612
     .  *(1.0d0+(tc/1.177d7)**0.5)**1.388)**(-1)
      goto  4
c
c calculate photo reaction rates
c

  3   xk(i)=xkp(i,1)*deplete
      if(xkp(i,2).ne.0.0d0) xk(i)=xk(i)*(tr/300.0d0)**xkp(i,2)
      if(xkp(i,3).ne.0.0d0) then
        if(xkp(i,4).gt.0.0d0) then
           xk(i)=xk(i)*exp(xkp(i,3)/tr)
        elseif(xkp(i,4).lt.0.0d0) then
           xk(i)=xk(i)*exp(tr/xkp(i,3))
        else
           print*,'problem with rate exp term'
        endif
      endif
c
c get equilibrium constant at radiation temperature
c 
      if(xkp(i,4).ne.0.0d0) then
         ki=idint(dabs(xkp(i,4)))
         kee=ke(ki,tr)*deplete 
c fix for LiH radiative assocition rate
      if((xkp(i,4).eq.-1.50d1).and.(tr.lt.380.0d0)) then
	xk(i)=3.74d-20*(tr/300.0d0)**0.11*dexp(-tr/2.0d3)*deplete
      endif
c uuses Li PI by detailed balance
c       if((dabs(xkp(i,4)).gt.5.0d0).or.(dabs(xkp(i,4)).eq.4.0d0))
       if(dabs(xkp(i,4)).gt.5.0d0)
     .    xk(i)=xk(i)*kee
cccccc
cc  print particular photo-rate
      if(i.eq.36) write(80,*) tr, xk(i),i
c         write(3,*) 'k called',i,xk(i),ki,kee
      else
         print*, 'no partition function for reaction ',i
      endif
  4   continue
c
c      write(3,999) 'rad temp',' = ',tr,'den',' = ',den
c      write(3,999) (nminp(j),' = ',xinp(j), j=1,idim)
  999 format(3(a5,a3,1pe10.2,2x))
c      write(3,*)
c      write(3,*) 'nreacs =',nreacs
c      do 105 i=1,nreacs
c      do 106 j=1,6
c      rch(j)=' '
c      if (j.gt.3) symbol(j-1)=' '
c      if (r(i,j).gt.0) rch(j)=name(r(i,j))
c      if (r(i,j).lt.0) rch(j)=nminp(-r(i,j))
c      if (j.gt.3 .and. rch(j).ne.' ') symbol(j-1)='  +'
c  106 continue
c      write(3,998) i,rch(1),symbol(1),rch(2),symbol(2),rch(3)
c     1 ,symbol(3),rch(4),symbol(4),rch(5),symbol(5),rch(6),xk(i)
  998 format (i4,1x,5(a7,a5),a5,1pe10.2)
  105 continue
      return
      end

