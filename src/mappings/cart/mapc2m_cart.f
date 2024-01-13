      subroutine fclaw_map_2d_c2m_cart(blockno,xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp
      integer blockno

c     # return a map in [-1,1],[-1,1]
      xp = 2*xc - 1
      yp = 2*yc - 1
      zp = 0


      end
