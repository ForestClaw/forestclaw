      subroutine mapc2m_twisted_torus(blockno,xc,yc,xp,yp,zp,alpha)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision alpha, r

      double precision pi, pi2
      common /compi/ pi, pi2

      r = 1.d0 + alpha*cos(pi2*(yc + xc))

      xp = r*cos(pi2*xc)
      yp = r*sin(pi2*xc)
      zp = alpha*sin(pi2*(yc + xc))

      end
