      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision r

      logical fclaw2d_map_is_used

      cont = fclaw_map_get_context()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)
      else
         xp = xc
         yp = yc
      endif

      r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2)

      fdisc = r-0.25d0
      end
