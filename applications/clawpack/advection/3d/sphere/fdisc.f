      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      double precision xc,yc, zc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      double precision pi, pi2
      common /compi/ pi, pi2

      double precision phi, R
      double precision rp, phi0

      cont = fclaw_map_get_context()

      call fclaw3d_map_c2m(cont,
     &      blockno,xc,yc,zc,xp,yp,zp)

      rp = sqrt(xp**2 + yp**2 + zp**2)

      phi = asin(zp/rp)

      phi0 = pi/4.d0
      fdisc = abs(phi - phi0) - pi/32.0


      end
