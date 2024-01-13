      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context
      double precision th, tp

      double precision pi, pi2
      common /compi/ pi, pi2

      cont = fclaw_map_get_context()

      call fclaw_map_2d_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

c     # Torus or annulus
      th = atan2(yp,xp)
      tp = abs(th)
      fdisc = tp - pi/4.d0

      end
