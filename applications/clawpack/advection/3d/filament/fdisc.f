      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      integer blockno
      double precision xc,yc, zc

      integer*8 cont, fclaw_map_get_context
      integer fclaw2d_map_is_used

      double precision r, xp, yp, zp

      cont = fclaw_map_get_context()


      if (fclaw2d_map_is_used(cont) .ne. 0) then
         call fclaw3d_map_c2m(cont,
     &         blockno,xc,yc,zc,xp,yp,zp)
      else
         xp = xc
         yp = yc
         zp = zc
      endif

      r = sqrt((xp-0.5d0)**2 + (yp-1.d0)**2 + (zp-0.5)**2)

      fdisc = r-0.25d0
      end
