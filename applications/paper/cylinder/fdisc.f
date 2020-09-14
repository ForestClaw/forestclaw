      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp, rp
      integer blockno


      double precision pi, pi2
      common /compi/ pi, pi2

      double precision r_cyl, h_cyl
      common /cylinder_comm/ r_cyl, h_cyl

      integer example
      common /example_comm/ example  

      DOUBLE PRECISION xc0, yc0, r0
      COMMON /cylinder_init_comm/ xc0, yc0, r0


      integer*8 cont, get_context
      double precision xp0, yp0, zp0, r

      cont = get_context()

      call fclaw2d_map_c2m(cont, blockno,xc,yc,xp,yp,zp)

      call fclaw2d_map_c2m(cont, blockno,xc0,yc0,xp0,yp0,zp0)

c     # Distance from thc
      r = sqrt((xp-xp0)**2 + (yp-yp0)**2 + (zp-zp0)**2)

      fdisc = r - r0 

      end
