      subroutine mapc2m_torus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision alpha, r

      double precision pi, pi2

      common /compi/ pi

      pi2 = 2*pi

      r = 1 + alpha*cos(pi2*yc)

      xp = r*cos(pi2*xc)
      yp = r*sin(pi2*xc)
      zp = alpha*sin(pi2*yc)

      end
