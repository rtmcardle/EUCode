       double precision x(3),xprime(3),pd(3,3),t
       integer i,j,n
       external input,rates,indx,fcn1,fcnj
c
       t=1000.0d1
c       call init(n,x)
       call input(n,x)
       call rates(t,n,x)
       write(6,*) (x(i), i=1,n)
       call fcn1(n,x,xprime)
       write(6,*) (xprime(i), i=1,n)
       call fcnj(n,x,pd)
       do 10 j=1,n
          write(6,*) (pd(i,j), i=1,n)
 10    continue 
       stop
       end
