      subroutine mapc2m_torus(xc,yc,xp,yp,zp)
      implicit none

      double precision xc,yc,xp,yp,zp

      double precision alpha, beta
      common /torus_comm/ alpha, beta

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision r1, R, theta, phi

      call map_comp2torus(xc,yc,theta,phi)      

      r1 = alpha*(1 + beta*sin(theta))
      R = 1 + r1*cos(phi)

      xp = R*cos(theta)
      yp = R*sin(theta)
      zp = r1*sin(phi)

      end
