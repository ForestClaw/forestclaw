      subroutine setaux_manifold(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,
     &      maux,aux,xp,yp,zp,xd,yd,zd,area)
      implicit none

      integer maxmx, maxmy,mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision  aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)


      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

      return
      end
