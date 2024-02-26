      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc
      integer blockno

      double precision pi, pi2
      common /compi/ pi, pi2

      integer*8 cont, fclaw_map_get_context
      double precision xp, yp, zp
      double precision th, tp

      cont = fclaw_map_get_context()

      pi = 4.0*atan(1.0)
      pi2 = 2*pi

      call fclaw_map_2d_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

c     # Torus or annulus
      th = atan2(yp,xp)
      tp = abs(th - pi/2.d0)
      fdisc = tp - pi/8.d0

      end
