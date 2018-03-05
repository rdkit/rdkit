      subroutine driver(n,m,nLevels,len,inData,d,iopt,ia,ib,crit)
      implicit none

c     locals
      integer n,m,nLevels
      real*8 inData(n,m)
      integer iopt
      integer ia(n),ib(n)
      real*8 crit(n)
      integer len
      real*8  d(len)

c     extern functions
      integer ioffset

c     locals
c      real*8 membr(n)
c      integer nn(n)
c      real*8 disnn(n)
c      logical flag(n)
      real*8 membr(10)
      integer nn(10)
      real*8 disnn(10)
      logical flag(10)

      integer i,j,k,idx
      real*8 dist

      do i=1,n
         do j = i+1,n
            idx = ioffset(n,i,j)
            dist = 0.0
            do k=1,m
               dist = dist + (inData(i,k)-inData(j,k))**2
            enddo
            d(idx) = dist
         enddo
      enddo
      CALL HC(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,FLAG,D)

      end
      

      subroutine distdriver(n,len,d,iopt,ia,ib,crit)
      implicit none

c     args
      integer n
      integer len
      real*8  d(len)
      integer iopt
      integer ia(n),ib(n)
      real*8 crit(n)

c     locals
c      real*8 membr(n)
c      integer nn(n)
c      real*8 disnn(n)
c      logical flag(n)
      real*8 membr(10)
      integer nn(10)
      real*8 disnn(10)
      logical flag(10)
      integer i

      CALL HC(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,FLAG,D)
c      write(6,*) "ia: ",(ia(i),i=1,N)
c      write(6,*) "ib: ",(ib(i),i=1,N)
c      write(6,*) "crit: ",(crit(i),i=1,N)
      end
      
