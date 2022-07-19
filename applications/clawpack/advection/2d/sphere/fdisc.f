      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision pi, pi2
      common /compi/ pi, pi2

      cont = fclaw_map_get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

c     #  Make a beach ball.
      if (xp .ne. 0) then
         fdisc = yp/xp
      else
         fdisc = 1
      endif

      end
