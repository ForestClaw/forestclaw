      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      integer blockno
      double precision xc,yc, zc

      integer*8 cont, fclaw_map_get_context
      integer fclaw2d_map_is_used

      double precision xp, yp, zp

      cont = fclaw_map_get_context()

      if (fclaw2d_map_is_used(cont) .ne. 0) then
         call fclaw3d_map_c2m(cont,
     &         blockno,xc,yc,zc,xp,yp,zp)
      else
         xp = xc
         yp = yc
         zp = zc
      endif

      fdisc = xp - 0.5

      end
