      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision phi0, phi, xp, yp, zp, rp

      cont = fclaw_map_get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

!!c     #  Make a beach ball.
!!      if (xp .ne. 0) then
!!         fdisc = yp/xp
!!      else
!!         fdisc = 1
!!      endif


      rp = sqrt(xp**2 + yp**2 + zp**2)

      phi = asin(zp/rp)

      phi0 = pi/4.d0
      fdisc = abs(phi - phi0) - pi/32.0

      end
