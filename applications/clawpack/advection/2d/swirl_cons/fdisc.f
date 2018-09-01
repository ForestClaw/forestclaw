      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      double precision dxc, xm, ym, w, wx, wy

      double precision r


      double precision pi, r2
      integer example

      common /compi/ pi
      common /comex/ example

      cont = get_context()

      if (fclaw2d_map_is_used(cont)) then
         call fclaw2d_map_c2m(cont,
     &         blockno,xc,yc,xp,yp,zp)
      else
         xp = xc
         yp = yc
      endif
  
      dxc = 0.25d0
      xm = 0.2d0
      ym = 0.5d0
      w = 0.3d0
      wx = 0.1d0
      wy = 0.35d0
c      r2 = (xp-0.5d0)**2 + (yp-0.5d0)**2
c      fdisc = r2 - 0.25**2
      if (abs(xp-xm) .le. wx .and. abs(yp-ym) .le. wy) then      
         fdisc = -1
      else
         fdisc = 1
      endif


      end
