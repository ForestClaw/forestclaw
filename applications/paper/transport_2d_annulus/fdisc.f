      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision x0, y0, r0, r, ravg, th

      call mapc2m_annulus(xc,yc,xp,yp,zp)

      ravg = (1 + beta)/2.d0
      th = pi2*(0.25 + 1.d0/32.d0)
      x0 = ravg*cos(th)
      y0 = ravg*sin(th)

      r = sqrt((xp - x0)**2 + (yp-y0)**2)

      r0 = init_radius
      fdisc = r - r0

      end
