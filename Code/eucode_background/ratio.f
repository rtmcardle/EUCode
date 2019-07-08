       integer i, j, n
       parameter (n=27)
       double precision z, a(n), dhprim,d2,dh,d2p,dhp,
     .     d3p,d2hp,dh2p
       parameter(dhprim=2.55d-5)
c
       open(unit=10,file='ratios.data')
c
       do i=1,5000
          read(16,*,end=89) z,(a(j),j=1,n)
          d2=a(24)/a(14)/dhprim**2
          dh=a(18)/a(14)/dhprim
          d2p=a(25)/a(13)/dhprim**2
          dhp=a(17)/a(13)/dhprim
          d3p=a(27)/a(16)/dhprim**3
          d2hp=a(26)/a(16)/dhprim**2
          dh2p=a(20)/a(16)/dhprim
          write(10,*) z,d2,dh,d2p,dhp,d3p,d2hp,dh2p 
        end do
 89     continue
c
        stop
        end
