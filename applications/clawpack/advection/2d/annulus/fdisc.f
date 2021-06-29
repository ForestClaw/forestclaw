      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc
      integer blockno

      double precision pi
      common /compi/ pi

      integer*8 cont, get_context
      double precision xp, yp, zp
      double precision th, tp

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

c     # Torus or annulus
      th = atan2(yp,xp)
      tp = abs(th - pi/2.d0)
      fdisc = tp - pi/8.d0

      end
