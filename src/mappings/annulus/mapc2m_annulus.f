      subroutine mapc2m_annulus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp, r, alpha
      double precision pi
      common /compi/ pi

      r = alpha + (1-alpha)*yc
      xp = r*cos(2*pi*xc)
      yp = r*sin(2*pi*xc)
      zp = 0


      end
