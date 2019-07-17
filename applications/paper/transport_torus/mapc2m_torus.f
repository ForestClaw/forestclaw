      subroutine mapc2m_torus(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision r1, R

      r1 = alpha*(1 + beta*sin(pi2*xc))
      R = 1 + r1*cos(pi2*yc)

      xp = R*cos(pi2*xc)
      yp = R*sin(pi2*xc)
      zp = r1*sin(pi2*yc)

      end
