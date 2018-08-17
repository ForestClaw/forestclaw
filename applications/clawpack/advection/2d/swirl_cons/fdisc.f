      double precision function  fdisc(blockno,xc,yc)
      implicit none

      double precision xc,yc, xp, yp, zp
      integer blockno
      integer*8 cont, get_context

      double precision r

      logical fclaw2d_map_is_used

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
  

c      if (abs(xp - 0.125d0) .gt. 0.025d0 .or. 
c      if (abs(xp - 0.375d0) .gt. 0.025d0 .or. 
c     &    abs(yp - 0.5d0) .gt. 0.125d0) then
      if (abs(xp - 0.375d0) .le. 0.125d0) then
         fdisc = 1
      else
         fdisc = -1
      endif

      end
