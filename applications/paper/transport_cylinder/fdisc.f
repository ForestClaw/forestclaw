      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno

      double precision pi, pi2
      common /compi/ pi, pi2

      integer example
      common /example_comm/ example  

      DOUBLE PRECISION xc0, yc0, r0
      COMMON /cylinder_init_comm/ xc0, yc0, r0

      double precision xp0,yp0,zp0,r

      call mapc2m_cylinder(xc,yc,xp,yp,zp)

      call mapc2m_cylinder(xc0,yc0,xp0,yp0,zp0)

c     # Distance from (xc0,yc0)
c      r = sqrt((xp-xp0)**2 + (yp-yp0)**2 + (zp-zp0)**2)      
c      fdisc = r - r0 

      r = abs(yc - .5d0)
      fdisc = r - 0.25

      end
