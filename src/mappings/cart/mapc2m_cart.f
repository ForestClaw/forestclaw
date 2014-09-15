      subroutine mapc2m_cart(blockno,xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno

c     # return a map in [-1,1],[-1,1]
      xp = 2.d0*(xc - 0.5d0)
      yp = 2.d0*(yc - 0.5d0)
      zp = 0


      end
