      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc
      integer blockno

      double precision r, hmax, b, c
      common /slotteddisk_parms/ r, hmax, b, c

      integer*8 cont, fclaw_map_get_context
      
      double precision q, qv, slotted_disk_sum
      double precision xp, yp, zp

      cont = fclaw_map_get_context()

      call fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)

c      call get_td_sdisk_parms(r,hmax,b,c)

c     # Returns 0 or 1.
      q = slotted_disk_sum(xp,yp,zp)

      qv = (q-b)/c
      fdisc = 1 - 2*qv
      end
