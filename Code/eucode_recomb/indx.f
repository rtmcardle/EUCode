c+++++++++++++++++++++++++++++++++++++++++++++++++++
c integer function to determine an index number for
c reactants and products of each reaction.
c
c returns integer indx
c-----------------------------------------------------
c written by S. Lepp ca 1984
c modified by P. C. Stancil, 6-25-96
c-----------------------------------------------------
c item   particle name to determine index for
c
c table1 contains list of gas particle names previously
c        read from species.data
c
c n1     total number of gas particle types
c
c table2 contains list of non-solve particles, photons
c        in this case
c
c n2     number of non-solve particles
c
c-----------------------------------------------------
      integer function indx(item,table1,n1,table2,n2)
      implicit double precision (a-h,o-z)
      integer n1,n2
      character*7 item,table1(n1),table2(n2)
c
c
      indx = 0
      if (item.eq.'       ') return
c
c assign index for gas particle
c
      do 101 i=1,n1
      if (item.eq.table1(i)) indx = i
  101 continue
      if (indx.ne.0) return
c
c assign index for non-gas particle, photon
c
      do 102 i=1,n2
      if (item.eq.table2(i)) indx = -i
  102 continue
      if (indx.ne.0) return
      write(6,*) 'item =',item
      write(6,*) 'can not find item'
      return
      end
