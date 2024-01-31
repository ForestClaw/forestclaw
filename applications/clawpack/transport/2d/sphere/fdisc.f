      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, phi, theta
      integer blockno

      double precision pi, pi2
      common /compi/ pi, pi2

      integer initchoice
      common /initchoice_comm/ initchoice

      integer*8 cont, fclaw_map_get_context

      double precision xp, yp, zp

c      double precision r, r0

      cont = fclaw_map_get_context()

      call fclaw_map_2d_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)

      call map2spherical(xp,yp,zp,theta,phi)

      fdisc = abs(theta-pi) - pi/6

c      r0 = 0.2
c      r = sqrt((xp-0.5)**2 + (yp-0.5)**2)
c      fdisc = r - r0


      end
