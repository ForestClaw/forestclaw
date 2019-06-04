      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision init_radius
      common /initradius_comm/ init_radius

      double precision beta
      common /annulus_comm/ beta      

      double precision x0,y0,r,r0, ravg

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

cc     # Torus or annulus
c      th = atan2(yp,xp)
c      tp = abs(th - pi/2.d0)
c      fdisc = tp - pi/8.d0

      r0 = init_radius
      ravg = (1 + beta)/2.d0
      x0 = -ravg
      y0 = 0

      r = sqrt((xp - x0)**2 + (yp-y0)**2)
      fdisc = r - r0


      end
