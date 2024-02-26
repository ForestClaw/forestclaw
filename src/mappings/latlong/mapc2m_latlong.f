c     # [xc,yc] are in [-180,180]x[0,360]
      subroutine fclaw_map_2d_c2m_latlong(blockno,xc,yc,xp,yp,zp)
      implicit none

      integer blockno
      double precision xc,yc,xp,yp,zp
      double precision deg2rad,xc1,yc1

      double precision pi, pi2
      common /compi/ pi, pi2

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
