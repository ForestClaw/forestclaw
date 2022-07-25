      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      integer blockno
      double precision xc,yc, zc

      integer*8 cont, fclaw_map_get_context

      double precision r, xp, yp, zp

      cont = fclaw_map_get_context()

      call fclaw3d_map_c2m(cont,
     &         blockno,xc,yc,zc,xp,yp,zp)

      r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2 + (zp-0.5)**2)

      fdisc = r-0.25d0
      end
