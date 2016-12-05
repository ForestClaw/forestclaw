      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno
      integer*8 cont, get_context
      double precision th, tp, r2

      double precision pi
      common /compi/ pi

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      th = atan2(yp,xp)
      tp = abs(th + pi/2.d0)
      fdisc = tp - pi/8.d0

      end
