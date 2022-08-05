      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      double precision xc,yc, zc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision maxelev 
      common /com_sphere/ maxelev

      double precision phi, R
      double precision rp, phi0

      cont = fclaw_map_get_context()

      call fclaw3d_map_c2m(cont,
     &      blockno,xc,yc,zc,xp,yp,zp)

      rp = sqrt(xp**2 + yp**2 + zp**2)

!!      xcenter = Rcenter*cos(phi)*cos(theta)
!!      ycenter = Rcenter*cos(phi)*sin(theta)
!!      zcenter = Rcenter*sin(phi)
!!
!!      r = sqrt((xp - xcenter)**2 + (yp - ycenter)**2)

      phi = asin(zp/rp)

      phi0 = pi/4.d0
      fdisc = abs(phi - phi0) - pi/32.0

!!c     #  Make a beach ball.
!!      if (xp .ne. 0) then
!!         fdisc = yp/xp
!!      else
!!         fdisc = 1
!!      endif

!!      if (xp .gt. 0) then
!!         fdisc = 1
!!      else
!!            fdisc = -1
!!      endif


      end
