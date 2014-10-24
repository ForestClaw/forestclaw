      subroutine mapc2m_brick(blockno,xc,yc,xp,yp,zp,mi,mj)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno, mi, mj, i, j

c     # Map into [0,mi] x [0,mj]
c     # Blocks are order i first, then j.
c     # blockno = j*mi + i
      j = blockno/mi
      i = blockno - j*mi
      xp = i + xc
      yp = j + yc
      zp = 0

      end
