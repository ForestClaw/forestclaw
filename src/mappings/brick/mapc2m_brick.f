      subroutine mapc2m_brick(blockno,xc,yc,xp,yp,zp,mi,mj)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno, mi, mj, i, j

      integer ix(0:15), jy(0:15)
      data ix /0,1,0,1,2,3,2,3,0,1,0,1,2,3,2,3/
      data jy /0,0,1,1,0,0,1,1,2,2,3,3,2,2,3,3/

c      if (mi .ne. 4 .or. mj .ne. 4) then
c         write(6,*) 'mapc2m_brick.f : Must have mi == mj == 4'
c         stop
c      endif

      i = ix(blockno)
      j = jy(blockno)
      xp = i + xc
      yp = j + yc
      zp = 0

      end
