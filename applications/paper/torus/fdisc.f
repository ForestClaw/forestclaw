      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno


      double precision pi, pi2
      common /compi/ pi, pi2

      double precision alpha, beta
      common /torus_comm/ alpha, beta
      
      integer*8 cont, get_context
      double precision th, r0,x0,y0,z0,r

      cont = get_context()

      call fclaw2d_map_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)

      th = pi2*(0.25 + 1.d0/16.d0)
      r0 = 0.15
      x0 = cos(th)
      y0 = sin(th)
      z0 = alpha

c     # Distance from thc
      r = sqrt((xp-x0)**2 + (yp-y0)**2 + (zp-z0)**2)

      fdisc = r - r0 

      end
