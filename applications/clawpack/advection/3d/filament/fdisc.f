      double precision function  fdisc(blockno,xc,yc,zc)
      implicit none

      integer blockno
      double precision xc,yc, zc

      double precision r

      r = sqrt((xc-0.5d0)**2 + (yc-1.d0)**2 + (zc-0.5)**2)

      fdisc = r-0.25d0
      end
