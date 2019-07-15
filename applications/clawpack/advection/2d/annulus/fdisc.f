      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision beta, theta(2)
      common /annulus_comm/ beta, theta 

      integer example
      common /example_comm/ example

      integer*8 cont, get_context

      double precision th, tp, x0,y0,r,r0, ravg, rd

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      r0 = init_radius
      ravg = (1 + beta)/2.d0
      th = (0.5 + 1.d0/16.d0)*pi
      x0 = ravg*cos(th)
      y0 = ravg*sin(th)
      rd = sqrt(xp**2 + yp**2)

      r = sqrt((xp - x0)**2 + (yp-y0)**2)
      fdisc = r - r0


      end
