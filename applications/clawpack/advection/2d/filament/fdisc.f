      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision r

      cont = fclaw_map_get_context()

      call fclaw_map_2d_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)

      r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2)

      fdisc = r-0.25d0
      end
