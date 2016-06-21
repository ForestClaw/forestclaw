      subroutine mapc2m_twisted_torus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision alpha, r

      double precision pi

      common /compi/ pi

      r = 1.d0 + alpha*cos(2.d0*pi*(yc + xc))

      xp = r*cos(2.d0*pi*xc)
      yp = r*sin(2.d0*pi*xc)
      zp = alpha*sin(2.d0*pi*(yc + xc))

      end
