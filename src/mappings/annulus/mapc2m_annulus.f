      subroutine mapc2m_annulus(blockno,xc,yc,xp,yp,zp,beta)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp, r, beta
      double precision pi
      common /compi/ pi

      r = beta + (1-beta)*yc
      xp = r*cos(2*pi*xc)
      yp = r*sin(2*pi*xc)
      zp = 0


      end
