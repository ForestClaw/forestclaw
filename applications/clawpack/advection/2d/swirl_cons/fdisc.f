      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context
      logical fclaw2d_map_is_used

      double precision dxc, xm, ym, w

      double precision r


      double precision pi
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
      xm = 0.5
      ym = 0.5
      w = dxc/2.d0
      if (abs(yp-ym) .le. w .and. abs(xp-xm) .le. w) then
         fdisc = -1
      else
         fdisc = 1
      endif


      end
