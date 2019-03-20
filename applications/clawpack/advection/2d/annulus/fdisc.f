      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp

      double precision pi, pi2
      common /compi/ pi, pi2

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

c     # Torus or annulus
      th = atan2(yp,xp)
      tp = abs(th - pi/2.d0)
      fdisc = tp - pi/8.d0

      end
