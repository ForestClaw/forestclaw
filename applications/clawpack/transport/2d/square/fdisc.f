      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc
      integer blockno

      integer example
      common /example_comm/ example

      integer*8 cont, fclaw_map_get_context
      logical fclaw2d_map_is_used

      double precision xp, yp, zp

      double precision r, r0


      cont = fclaw_map_get_context()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)
      else
         xp = xc
         yp = yc
      endif

      r0 = 0.2
      r = sqrt((xp-0.5)**2 + (yp-0.5)**2)
      fdisc = r - r0


      end
