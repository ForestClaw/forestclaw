      subroutine mapc2m_torus(xc,yc,xp,yp,zp,alpha, beta)
      implicit none

      double precision xc,yc,xp,yp,zp
      double precision alpha, beta, r1, R


      double precision pi, pi2
      common /compi/ pi, pi2

      r1 = alpha*(1 + beta*sin(pi2*xc))
      R = 1 + r1*cos(pi2*yc)

      xp = R*cos(pi2*xc)
      yp = R*sin(pi2*xc)
      zp = r1*sin(pi2*yc)

      end

