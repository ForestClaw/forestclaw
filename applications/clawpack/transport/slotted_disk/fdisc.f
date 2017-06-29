      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      double precision q, slotted_disk_sum
      integer*8 cont, get_context
      double precision r,hmax,b,c,qv

      cont = get_context()

      call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)

      call get_td_sdisk_parms(r,hmax,b,c)

c     # Returns 0 or 1.
      q = slotted_disk_sum(xp,yp,zp)

      qv = (q-b)/c
      fdisc = 1 - 2*qv
      end
