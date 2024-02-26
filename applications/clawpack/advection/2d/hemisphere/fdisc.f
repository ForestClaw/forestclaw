      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, fclaw_map_get_context

      cont = fclaw_map_get_context()

      call fclaw_map_2d_c2m(cont,
     &      blockno,xc,yc,xp,yp,zp)


      fdisc = abs(xp) - 0.25d0
      end
