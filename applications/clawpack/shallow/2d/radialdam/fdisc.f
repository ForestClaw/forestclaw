      double precision function fdisc(blockno, x,y)
      implicit none

      integer blockno
      double precision x,y

      integer*8 cont, fclaw_map_get_context
      logical fclaw2d_map_is_used

      double precision x0,y0,r0
      common /cdisc/ x0,y0,r0

      double precision r2,xp,yp,zp

      cont = fclaw_map_get_context()

      IF (fclaw2d_map_is_used(cont)) THEN
         CALL fclaw2d_map_c2m(cont,blockno,x,y,xp,yp,zp)
         r2 = (xp-x0)**2 + (yp-y0)**2
      ELSE
         r2 = (x-x0)**2 + (y-y0)**2
      ENDIF

      fdisc = r2 - r0**2
c
      return
      end
