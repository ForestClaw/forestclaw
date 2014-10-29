      subroutine mapc2m_latlong(blockno,xc,yc,xp,yp,zp)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision pi, deg2rad,xc1,yc1

      common /compi/ pi

c     # blockno is not used here;  assume that [xc,yc] is
c     # in a box [long0, long1]x[lat0,lat1], in
c     # degrees.

      deg2rad = pi/180.d0

      xc1 = deg2rad*xc
      yc1 = deg2rad*yc

      xp = cos(yc1)*cos(xc1)
      yp = cos(yc1)*sin(xc1)
      zp = sin(yc1)

      end
