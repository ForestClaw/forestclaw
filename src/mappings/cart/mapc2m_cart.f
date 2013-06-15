      subroutine mapc2m_cart(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      xp = xc
      yp = xc + yc
      zp = 0


      end
