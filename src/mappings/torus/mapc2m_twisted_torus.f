      subroutine mapc2m_torus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision alpha, r

      double precision pi

      common /compi/ pi

      r = 1 + alpha*cos(2*pi*(yc + xc))

      xp = r*cos(2*pi*xc)
      yp = r*sin(2*pi*xc)
      zp = alpha*sin(2*pi*(yc + xc))

      end
